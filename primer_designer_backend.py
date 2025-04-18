from fastapi import FastAPI, HTTPException, Query
from fastapi.middleware.cors import CORSMiddleware
import re

app = FastAPI()

# CORS 허용 (웹페이지에서 접근 가능하도록)
app.add_middleware(
    CORSMiddleware,
    allow_origins=["*"],
    allow_credentials=True,
    allow_methods=["*"],
    allow_headers=["*"],
)

# 아미노산 → 가장 흔한 codon (사람 기준)
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

# ✅ 실제 UniProt에서 CDS를 못 가져오므로, 예시용 DNA 직접 넣기
def fetch_cds_from_uniprot(uniprot_id: str) -> str:
    if uniprot_id == "P17676":
        return (
            "ATGGACTACAAGGACGACGATGACAAGGCCGCCATCGACTTCGCCCCGTACCTGGAGCCGCCATCCTGCCGAGCCAGCAGCAGCGAGCAGCTGGGCCCCTGCTGCCCCAGCAGCAGCAG"
            "CCAGCGCCTGCCTGCTGCTGGAGCTGGTGCAGCGGGACGGCCTGCAGGTGGAGGGGGCCGAGCGCGCAGCGCGAGTGCAGGAGCGGGGCCGCCCTCTGGGGAGGGGGGACCCGAGGA"
            "AGGGGGAGCGGCGGGGGCTCCGGGTGAGGGAGCCGGAGTGGAGGCGGAGGAGGAGGCGCTGGAGCCCTGGGCTGGAGCGGGGGCGGCAGCGCCGCCACCCACAGCGCCCGGGGCGAG"
            "CGGCCCTGCTGCTGCGGGCTGCTGCCGGCGGCGGCGGCGGCGGCAGGGCGGCAGCGGCGGCAGCGGCGGCGGGCGGCGGCGGCGGCGGCGGCAGCGGCGGCAGGGCGGCGGCAGCGG"
            "CGGCGGCGGCGGCAGCGGCGGCAGCGGCGGCGGCGGCGGCAGCGGCGGCGGCAGCGGCGGCAGCGGCGGCGGCGGCGGCAGCGGCGGCAGCGGCGGCAGCGGCGGCAGCGGCGGCAG"
            "CGGCGGCAGCGGCGGCGGCAGCGGCGGCGGCAGCGGCGGCAGCGGCGGCAGCGGCGGCAGCGGCGGCAGCGGCGGCAGCGGCGGCGGCAGCGGCGGCAGCGGCGGCAGCGGCGGCGG"
            "CAGCGGCGGCAGCGGCGGCGGCAGCGGCGGCAGCGGCGGCGGCAGCGGCGGCAGCGGCGGCAGCGGCGGCGGCAGCGGCGGCAGCGGCGGCAGCGGCGGCAGCGGCGGCAGCGGCGG"
            "CAGCGGCGGCAGCGGCGGCAGCGGCGGCAGCGGCGGCAGCGGCGGCGGCAGCGGCGGCGGCAGCGGCGGCAGCGGCGGCGGCAGCGGCGGCAGCGGCGGCAGCGGCGGCAGCGGCGG"
        )
    raise HTTPException(status_code=404, detail="Only P17676 is supported for now.")

def design_primer(cds: str, mutation: str, flank: int = 15):
    match = re.match(r"([A-Za-z])(\d+)([A-Za-z*])", mutation)
    if not match:
        raise ValueError("Mutation 형식은 예: S76A")

    orig_aa, pos, new_aa = match.groups()
    pos = int(pos)
    codon_start = (pos - 1) * 3

    if codon_start + 3 > len(cds):
        raise ValueError("변이 위치가 CDS 길이를 초과합니다")

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
