# FastAPI: Primer Design from RefSeq NP_XXXXXX directly

from fastapi import FastAPI, HTTPException, Query
from fastapi.middleware.cors import CORSMiddleware
import requests, re
from Bio.Seq import Seq
from xml.etree import ElementTree

app = FastAPI()

app.add_middleware(
    CORSMiddleware,
    allow_origins=["*"],
    allow_credentials=True,
    allow_methods=["*"],
    allow_headers=["*"],
)

PREFERRED_CODON = {
    "A": "GCC", "R": "CGT", "N": "AAC", "D": "GAC", "C": "TGC",
    "Q": "CAG", "E": "GAG", "G": "GGC", "H": "CAC", "I": "ATC",
    "L": "CTG", "K": "AAG", "M": "ATG", "F": "TTC", "P": "CCC",
    "S": "AGC", "T": "ACC", "W": "TGG", "Y": "TAC", "V": "GTG",
}

def reverse_complement(seq: str) -> str:
    return str(Seq(seq).reverse_complement())

def gc_content(seq: str) -> float:
    seq = seq.upper()
    return (seq.count("G") + seq.count("C")) / len(seq) * 100

def calc_tm(seq: str) -> float:
    seq = seq.upper()
    return 2 * (seq.count("A") + seq.count("T")) + 4 * (seq.count("G") + seq.count("C"))

def get_mrna_id_from_protein(refseq_protein_id: str) -> str:
    url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/elink.fcgi"
    params = {
        "dbfrom": "protein",
        "db": "nuccore",
        "id": refseq_protein_id,
        "retmode": "xml"
    }
    r = requests.get(url, params=params)
    root = ElementTree.fromstring(r.content)
    ids = root.findall(".//LinkSetDb/Link/Id")
    if not ids:
        raise ValueError("연결된 mRNA (NM_) ID가 없습니다.")
    return ids[0].text

def get_cds(nm_id: str) -> str:
    url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"
    params = {
        "db": "nuccore",
        "id": nm_id,
        "rettype": "fasta",
        "retmode": "text"
    }
    r = requests.get(url, params=params)
    lines = r.text.splitlines()
    return ''.join([l.strip() for l in lines if not l.startswith(">")]).upper()

def get_protein_seq_from_np(refseq_protein_id: str) -> str:
    url = f"https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"
    params = {
        "db": "protein",
        "id": refseq_protein_id,
        "rettype": "fasta",
        "retmode": "text"
    }
    r = requests.get(url, params=params)
    lines = r.text.strip().split("\n")
    return ''.join(lines[1:]).strip()

def design_primer(cds: str, protein_seq: str, mutation: str, flank: int = 15):
    match = re.match(r"([A-Za-z])(\d+)([A-Za-z*])", mutation)
    if not match:
        raise ValueError("Mutation 형식은 S76A처럼 입력해주세요.")
    orig_aa, pos, new_aa = match.groups()
    pos = int(pos)

    trimmed_cds = cds[:len(cds) - len(cds) % 3]
    translated = str(Seq(trimmed_cds).translate())

    if pos > len(translated):
        raise ValueError("단백질 길이를 초과한 위치입니다.")
    if translated[pos - 1] != orig_aa.upper():
        raise ValueError(f"{pos}번째 아미노산은 {translated[pos - 1]}이며, {orig_aa}와 일치하지 않습니다.")

    codon_start = (pos - 1) * 3
    new_codon = PREFERRED_CODON.get(new_aa.upper())
    if not new_codon:
        raise ValueError("지원되지 않는 새로운 아미노산입니다.")

    up = cds[codon_start - flank: codon_start]
    down = cds[codon_start + 3: codon_start + 3 + flank]
    fwd = up + new_codon + down
    rev = reverse_complement(fwd)
    return {
        "forward_primer": fwd,
        "reverse_primer": rev,
        "tm": round(calc_tm(fwd), 1),
        "gc_percent": round(gc_content(fwd), 1)
    }

@app.get("/primer")
def primer_endpoint(
    refseq_protein: str = Query(...),
    mutation: str = Query(...)
):
    try:
        protein_seq = get_protein_seq_from_np(refseq_protein)
        nm_id = get_mrna_id_from_protein(refseq_protein)
        cds = get_cds(nm_id)
        primer = design_primer(cds, protein_seq, mutation)
    except Exception as e:
        import traceback
        traceback.print_exc()
        raise HTTPException(status_code=500, detail=str(e))
    return {
        "refseq_protein": refseq_protein,
        "mutation": mutation.upper(),
        **primer
    }