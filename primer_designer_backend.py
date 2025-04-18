###############################################
# primer_designer_backend.py  ―  FastAPI 서버 #
###############################################
from fastapi import FastAPI, HTTPException, Query
from fastapi.middleware.cors import CORSMiddleware
import re, requests
from Bio import Entrez, SeqIO
from Bio.Seq import Seq

app = FastAPI()

# CORS 허용
app.add_middleware(
    CORSMiddleware,
    allow_origins=["*"],
    allow_credentials=True,
    allow_methods=["*"],
    allow_headers=["*"],
)

Entrez.email = "your_email@example.com"   # ✏️ 본인 이메일로 변경

# 가장 흔한 코돈(휴먼)
PREFERRED_CODON = {
    "A": "GCC", "R": "CGT", "N": "AAC", "D": "GAC", "C": "TGC",
    "Q": "CAG", "E": "GAG", "G": "GGC", "H": "CAC", "I": "ATC",
    "L": "CTG", "K": "AAG", "M": "ATG", "F": "TTC", "P": "CCC",
    "S": "AGC", "T": "ACC", "W": "TGG", "Y": "TAC", "V": "GTG",
    "*": "TAA"
}

# ────────── 유틸 함수 ──────────
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

# NP → NM 변환
def get_nm_from_np(np_id: str) -> str:
    url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/elink.fcgi"
    params = {
        "dbfrom": "protein", "db": "nuccore",
        "id": np_id, "linkname": "protein_nuccore_mrna",
        "retmode": "json"
    }
    data = requests.get(url, params=params).json()
    try:
        return data['linksets'][0]['linksetdbs'][0]['links'][0]
    except Exception:
        raise HTTPException(500, "해당 NP ID에 NM ID를 찾지 못했습니다.")

# NM → CDS + 단백질 서열
def fetch_cds_and_protein(np_id: str):
    nm_id = get_nm_from_np(np_id)
    gb = SeqIO.read(
        Entrez.efetch(db="nuccore", id=nm_id, rettype="gb", retmode="text"),
        "gb")
    for feat in gb.features:
        if feat.type == "CDS" and "translation" in feat.qualifiers:
            cds_seq = feat.extract(gb.seq)
            prot = feat.qualifiers["translation"][0]
            return str(cds_seq), prot
    raise HTTPException(500, "NM 기록에서 CDS를 찾지 못했습니다.")

# 프라이머 설계
def design_primer(cds: str, protein: str, mutation: str, flank: int = 15):
    m = re.match(r"([A-Za-z])(\\d+)([A-Za-z*])", mutation)
    if not m:
        raise ValueError("변이 형식 오류 (예: S76A)")
    wt_aa, pos, mut_aa = m.groups()
    pos = int(pos)

    if pos > len(protein):
        raise ValueError("변이 위치가 단백질 길이를 초과합니다.")
    if protein[pos-1] != wt_aa:
        raise ValueError(f"{pos}번째 아미노산은 {protein[pos-1]}이며, {wt_aa}와 다릅니다.")

    codon_start = (pos-1)*3
    wt_codon   = cds[codon_start:codon_start+3]
    mut_codon  = PREFERRED_CODON[mut_aa.upper()]
    up   = cds[codon_start-flank: codon_start]
    down = cds[codon_start+3  : codon_start+3+flank]

    forward = up + mut_codon + down
    reverse = reverse_complement(forward)

    return {
        "forward_primer": forward,
        "reverse_primer": reverse,
        "tm": round(calc_tm(forward), 1),
        "gc_percent": round(gc_content(forward), 1),
        # 강조용 추가 필드
        "wt_codon": wt_codon,
        "mut_codon": mut_codon,
        "codon_offset": len(up)          # forward primer 안에서 변이 코돈 시작 index
    }

# ────────── 엔드포인트 ──────────
@app.get("/primer")
def primer_endpoint(
    refseq_protein: str = Query(..., min_length=6),
    mutation: str       = Query(..., regex=r"^[A-Za-z]\\d+[A-Za-z*]$")
):
    try:
        cds, prot = fetch_cds_and_protein(refseq_protein)
        result = design_primer(cds, prot, mutation)
        return {"refseq_protein": refseq_protein,
                "mutation": mutation.upper(),
                **result}
    except Exception as e:
        raise HTTPException(500, detail=str(e))
