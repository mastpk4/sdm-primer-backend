from fastapi import FastAPI, HTTPException, Query
from fastapi.middleware.cors import CORSMiddleware
import requests
from Bio import Entrez, SeqIO
from Bio.Seq import Seq
import re

app = FastAPI()

app.add_middleware(
    CORSMiddleware,
    allow_origins=["*"], allow_credentials=True, allow_methods=["*"], allow_headers=["*"],
)

Entrez.email = "your_email@example.com"  # 자신의 이메일로 설정

PREFERRED_CODON = {
    "A": "GCC", "R": "CGT", "N": "AAC", "D": "GAC", "C": "TGC",
    "Q": "CAG", "E": "GAG", "G": "GGC", "H": "CAC", "I": "ATC",
    "L": "CTG", "K": "AAG", "M": "ATG", "F": "TTC", "P": "CCC",
    "S": "AGC", "T": "ACC", "W": "TGG", "Y": "TAC", "V": "GTG",
}

def reverse_complement(seq: str) -> str:
    return str(Seq(seq).reverse_complement())

def gc_content(seq: str) -> float:
    return (seq.count("G") + seq.count("C")) / len(seq) * 100

def calc_tm(seq: str) -> float:
    return 2 * (seq.count("A") + seq.count("T")) + 4 * (seq.count("G") + seq.count("C"))

def get_nm_from_np(np_id: str) -> str:
    res = requests.get("https://eutils.ncbi.nlm.nih.gov/entrez/eutils/elink.fcgi", params={
        "dbfrom": "protein", "db": "nuccore", "id": np_id,
        "linkname": "protein_nuccore_mrna", "retmode": "json"
    })
    data = res.json()
    try:
        return data['linksets'][0]['linksetdbs'][0]['links'][0]
    except:
        raise HTTPException(status_code=500, detail="해당 NP ID에 대한 NM ID를 찾을 수 없습니다.")

def fetch_cds_and_protein(np_id: str):
    nm_id = get_nm_from_np(np_id)
    handle = Entrez.efetch(db="nuccore", id=nm_id, rettype="gb", retmode="text")
    record = SeqIO.read(handle, "gb")
    handle.close()
    for feature in record.features:
        if feature.type == "CDS" and "translation" in feature.qualifiers:
            cds = feature.extract(record.seq)
            protein = feature.qualifiers["translation"][0]
            return str(cds), protein
    raise HTTPException(status_code=500, detail="CDS를 찾을 수 없습니다.")

def design_primer(cds: str, protein: str, mutation: str, flank: int = 15):
    match = re.match(r"([A-Za-z])(\d+)([A-Za-z*])", mutation)
    if not match:
        raise ValueError("변이 형식이 잘못되었습니다 (예: S76A)")
    orig_aa, pos, new_aa = match.groups()
    pos = int(pos)

    if pos > len(protein):
        raise ValueError("변이 위치가 단백질 길이를 초과합니다.")
    if protein[pos - 1] != orig_aa:
        raise ValueError(f"{pos}번째 아미노산은 {protein[pos - 1]}이며, {orig_aa}와 일치하지 않습니다.")

    codon_start = (pos - 1) * 3
    new_codon = PREFERRED_CODON[new_aa.upper()]
    up = cds[codon_start - flank:codon_start]
    down = cds[codon_start + 3:codon_start + 3 + flank]
    fwd = up + new_codon + down
    rev = reverse_complement(fwd)

    # 하이라이트 적용된 서열
    wt_codon = cds[codon_start:codon_start + 3]
    wt_dna_context = cds[codon_start - 30: codon_start] + f"<b>{wt_codon}</b>" + cds[codon_start + 3: codon_start + 33]
    wt_protein_context = (
        protein[max(0, pos - 11): pos - 1]
        + f"<b>{orig_aa}</b>"
        + protein[pos: min(len(protein), pos + 10)]
    )
    fwd_display = up + f"<b>{new_codon}</b>" + down
    rev_display = reverse_complement(up + f"<b>{new_codon}</b>" + down)

    return {
        "forward_primer": fwd_display,
        "reverse_primer": rev_display,
        "tm": round(calc_tm(fwd), 1),
        "gc_percent": round(gc_content(fwd), 1),
        "wt_protein_context": wt_protein_context,
        "wt_dna_context": wt_dna_context,
    }

@app.get("/primer")
def primer_endpoint(
    refseq_protein: str = Query(...),
    mutation: str = Query(..., pattern=r"^[A-Za-z]\d+[A-Za-z*]$")
):
    try:
        cds, protein_seq = fetch_cds_and_protein(refseq_protein)
        primer = design_primer(cds, protein_seq, mutation)
    except Exception as e:
        raise HTTPException(status_code=500, detail=str(e))
    return {
        "refseq_protein": refseq_protein,
        "mutation": mutation.upper(),
        **primer
    }
