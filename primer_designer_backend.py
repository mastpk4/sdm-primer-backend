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

CODON_TABLE = {
    'TTT':'F','TTC':'F','TTA':'L','TTG':'L','CTT':'L','CTC':'L','CTA':'L','CTG':'L',
    'ATT':'I','ATC':'I','ATA':'I','ATG':'M','GTT':'V','GTC':'V','GTA':'V','GTG':'V',
    'TCT':'S','TCC':'S','TCA':'S','TCG':'S','CCT':'P','CCC':'P','CCA':'P','CCG':'P',
    'ACT':'T','ACC':'T','ACA':'T','ACG':'T','GCT':'A','GCC':'A','GCA':'A','GCG':'A',
    'TAT':'Y','TAC':'Y','TAA':'*','TAG':'*','CAT':'H','CAC':'H','CAA':'Q','CAG':'Q',
    'AAT':'N','AAC':'N','AAA':'K','AAG':'K','GAT':'D','GAC':'D','GAA':'E','GAG':'E',
    'TGT':'C','TGC':'C','TGA':'*','TGG':'W','CGT':'R','CGC':'R','CGA':'R','CGG':'R',
    'AGT':'S','AGC':'S','AGA':'R','AGG':'R','GGT':'G','GGC':'G','GGA':'G','GGG':'G',
}

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
    url = f"https://rest.uniprot.org/uniprotkb/{uniprot_id}.fasta"
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
    codon = cds[codon_start:codon_start+3]
    if len(codon) != 3:
        raise ValueError("Invalid position")
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