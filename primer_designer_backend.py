# 🔁 primer_designer_backend.py

from fastapi import FastAPI, HTTPException, Query
from fastapi.middleware.cors import CORSMiddleware
from Bio import Entrez, SeqIO
import re

app = FastAPI()

# CORS 허용
app.add_middleware(
    CORSMiddleware,
    allow_origins=["*"],
    allow_credentials=True,
    allow_methods=["*"],
    allow_headers=["*"],
)

# 가장 흔한 codon 테이블
PREFERRED_CODON = {
    "A": "GCC", "R": "CGT", "N": "AAC", "D": "GAC", "C": "TGC",
    "Q": "CAG", "E": "GAG", "G": "GGC", "H": "CAC", "I": "ATC",
    "L": "CTG", "K": "AAG", "M": "ATG", "F": "TTC", "P": "CCC",
    "S": "AGC", "T": "ACC", "W": "TGG", "Y": "TAC", "V": "GTG",
}

# 보조 함수들
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

# 🧬 NCBI에서 CDS & 단백질 직접 가져오기
def fetch_cds_and_protein(np_id: str):
    Entrez.email = "your_email@example.com"
    handle = Entrez.elink(dbfrom="protein", db="nuccore", id=np_id, linkname="protein_nuccore_refseq")
    record = Entrez.read(handle)
    linked = record[0]["LinkSetDb"]
    if not linked:
        raise HTTPException(status_code=404, detail="NM ID를 찾을 수 없습니다")
    nm_id = linked[0]['Link'][0]['Id']

    # 유전체 DB에서 CDS 추출
    handle = Entrez.efetch(db="nuccore", id=nm_id, rettype="gb", retmode="text")
    gb_record = SeqIO.read(handle, "genbank")
    for feature in gb_record.features:
        if feature.type == "CDS" and "protein_id" in feature.qualifiers:
            if np_id in feature.qualifiers["protein_id"]:
                cds = feature.extract(gb_record.seq)
                protein_seq = feature.qualifiers.get("translation", [""])[0]
                return str(cds), protein_seq
    raise HTTPException(status_code=404, detail="CDS not found")

# 🧪 프라이머 설계
def design_primer(cds: str, protein_seq: str, mutation: str, flank: int = 15):
    match = re.match(r"([A-Za-z])(\d+)([A-Za-z\*])", mutation)
    if not match:
        raise ValueError("Mutation 형식은 예: S76A")

    orig_aa, pos, new_aa = match.groups()
    pos = int(pos)
    aa_idx = pos - 1
    codon_start = aa_idx * 3

    # 오류 처리
    translated = str(SeqIO.Seq.Seq(cds).translate())
    if aa_idx >= len(translated):
        raise ValueError("변이 위치가 CDS 범위를 초과합니다")
    if translated[aa_idx] != orig_aa:
        raise ValueError(f"{pos}번째 아미노산은 {translated[aa_idx]}이며, {orig_aa}와 일치하지 않습니다.")

    # 단백질 컨텍스트 (소문자로 표시)
    left = protein_seq[max(0, aa_idx-10):aa_idx]
    right = protein_seq[aa_idx+1:aa_idx+11]
    protein_context = f"{left}{orig_aa.lower()}{right}"

    # Wild-type DNA 서열
    dna_left = cds[max(0, codon_start-30):codon_start]
    wt_codon = cds[codon_start:codon_start+3]
    dna_right = cds[codon_start+3:codon_start+33]
    wt_dna = f"{dna_left}{wt_codon.lower()}{dna_right}"

    # Mutant codon
    new_codon = PREFERRED_CODON[new_aa.upper()].lower()
    up = cds[codon_start - flank : codon_start]
    down = cds[codon_start + 3 : codon_start + 3 + flank]
    fwd = up + new_codon + down
    rev = reverse_complement(fwd)

    return {
        "protein_context": protein_context,
        "wt_dna": wt_dna,
        "forward_primer": fwd,
        "reverse_primer": rev,
        "tm": round(calc_tm(fwd), 1),
        "gc_percent": round(gc_content(fwd), 1),
    }

@app.get("/primer")
def primer_endpoint(
    refseq_protein: str = Query(..., min_length=5, max_length=20),
    mutation: str = Query(..., regex=r"^[A-Za-z]\d+[A-Za-z\*]$")
):
    try:
        cds, protein_seq = fetch_cds_and_protein(refseq_protein)
        primer = design_primer(cds, protein_seq, mutation)
    except Exception as e:
        raise HTTPException(status_code=500, detail=str(e))
    return {
        "refseq_protein": refseq_protein,
        "mutation": mutation,
        **primer
    }