from fastapi import FastAPI, HTTPException, Query
from fastapi.middleware.cors import CORSMiddleware
import requests
import re

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
    trans = str.maketrans("ATCGatcg", "TAGCtagc")
    return seq.translate(trans)[::-1]

def gc_content(seq: str) -> float:
    seq = seq.upper()
    return (seq.count("G") + seq.count("C")) / len(seq) * 100

def calc_tm(seq: str) -> float:
    seq = seq.upper()
    at = seq.count("A") + seq.count("T")
    gc = seq.count("G") + seq.count("C")
    return 2 * at + 4 * gc

def fetch_cds_from_uniprot(uniprot_id: str) -> str:
    url = f"https://rest.uniprot.org/uniprotkb/{uniprot_id}.fasta?query=cds"
    r = requests.get(url)
    if r.status_code != 200:
        raise HTTPException(status_code=404, detail="UniProt ID not found")
    fasta = r.text.splitlines()
    seq = "".join(line.strip() for line in fasta if not line.startswith(">"))
    return seq.upper()

def design_primer(cds: str, mutation: str, flank: int = 15):
    match = re.match(r"([A-Za-z])(\d+)([A-Za-z*])", mutation)
    if not match:
        raise ValueError("Mutation format must be like S76A")
    orig_aa, pos, new_aa = match.groups()
    pos = int(pos)
    codon_start = (pos - 1) * 3
    if codon_start + 3 > len(cds):
        raise ValueError("Position out of range")

    new_codon = PREFERRED_CODON[new_aa.upper()]
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
    uniprot_id: str = Query(..., min_length=6, max_length=10),
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
