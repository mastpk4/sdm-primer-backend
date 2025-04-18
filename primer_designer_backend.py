from fastapi import FastAPI, HTTPException, Query
from fastapi.middleware.cors import CORSMiddleware
import re
import requests
from xml.etree import ElementTree
from Bio.Seq import Seq  # Biopython 필요

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
    at = seq.count("A") + seq.count("T")
    gc = seq.count("G") + seq.count("C")
    return 2 * at + 4 * gc

def get_refseq_protein_from_uniprot(uniprot_id: str) -> str:
    url = f"https://rest.uniprot.org/uniprotkb/{uniprot_id}.json"
    r = requests.get(url)
    if r.status_code != 200:
        raise ValueError("UniProt ID not found")
    data = r.json()
    for ref in data.get("uniProtKBCrossReferences", []):
        if ref.get("database") == "RefSeq":
            return ref.get("id")
    raise ValueError("RefSeq 단백질 ID를 찾을 수 없습니다.")

def get_nucleotide_id_from_protein(refseq_protein_id: str) -> str:
    url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/elink.fcgi"
    params = {
        "dbfrom": "protein",
        "db": "nuccore",
        "id": refseq_protein_id,
        "retmode": "xml"
    }
    r = requests.get(url, params=params)
    root = ElementTree.fromstring(r.content)
    linksets = root.findall(".//LinkSetDb/Link/Id")
    if not linksets:
        raise ValueError("연결된 mRNA ID를 찾을 수 없습니다")
    return linksets[0].text

def get_cds_sequence(nucleotide_id: str) -> str:
    url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"
    params = {
        "db": "nuccore",
        "id": nucleotide_id,
        "rettype": "fasta",
        "retmode": "text"
    }
    r = requests.get(url, params=params)
    lines = r.text.splitlines()
    return "".join(line.strip() for line in lines if not line.startswith(">")).upper()

def fetch_cds_from_uniprot(uniprot_id: str) -> str:
    refseq_protein = get_refseq_protein_from_uniprot(uniprot_id)
    nucleotide_id = get_nucleotide_id_from_protein(refseq_protein)
    cds = get_cds_sequence(nucleotide_id)
    return cds

def design_primer(cds: str, mutation: str, flank: int = 15):
    match = re.match(r"([A-Za-z])(\d+)([A-Za-z*])", mutation)
    if not match:
        raise ValueError("Mutation 형식은 예: S76A")
    orig_aa, pos, new_aa = match.groups()
    pos = int(pos)

    translated = str(Seq(cds).translate(to_stop=True))
    if pos > len(translated):
        raise ValueError("변이 위치가 단백질 길이를 초과합니다.")
    if translated[pos - 1] != orig_aa.upper():
        raise ValueError(f"{pos}번째 아미노산은 {translated[pos - 1]}이며, {orig_aa}와 일치하지 않습니다.")

    codon_start = (pos - 1) * 3
    new_codon = PREFERRED_CODON.get(new_aa.upper())
    if not new_codon:
        raise ValueError("지원되지 않는 아미노산입니다.")

    up = cds[codon_start - flank : codon_start]
    down = cds[codon_start + 3 : codon_start + 3 + flank]
    fwd = up + new_codon + down
    rev = reverse_complement(fwd)
    return {
        "forward_primer": fwd,
        "reverse_primer": rev,
        "tm": round(calc_tm(fwd), 1),
        "gc_percent": round(gc_content(fwd), 1),
    }

@app.get("/primer")
def primer_endpoint(
    uniprot_id: str = Query(..., min_length=6, max_length=12),
    mutation: str = Query(..., regex=r"^[A-Za-z]\d+[A-Za-z*]$")
):
    try:
        cds = fetch_cds_from_uniprot(uniprot_id)
        primer_data = design_primer(cds, mutation)
    except Exception as e:
        raise HTTPException(status_code=500, detail=str(e))
    return {
        "uniprot_id": uniprot_id,
        "mutation": mutation.upper(),
        **primer_data
    }
