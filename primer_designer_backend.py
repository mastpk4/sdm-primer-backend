from fastapi import FastAPI, HTTPException, Query
from fastapi.middleware.cors import CORSMiddleware
import requests, re
from Bio import Entrez, SeqIO
from Bio.Seq import Seq

app = FastAPI()

# CORS 설정
app.add_middleware(
    CORSMiddleware,
    allow_origins=["*"],
    allow_credentials=True,
    allow_methods=["*"],
    allow_headers=["*"],
)

Entrez.email = "your_email@example.com"  # 여기에 본인 이메일 입력

# 가장 흔한 codon table (인간 기준)
PREFERRED_CODON = {
    "A": "GCC", "R": "CGT", "N": "AAC", "D": "GAC", "C": "TGC",
    "Q": "CAG", "E": "GAG", "G": "GGC", "H": "CAC", "I": "ATC",
    "L": "CTG", "K": "AAG", "M": "ATG", "F": "TTC", "P": "CCC",
    "S": "AGC", "T": "ACC", "W": "TGG", "Y": "TAC", "V": "GTG",
}

# Reverse complement
def reverse_complement(seq: str) -> str:
    return str(Seq(seq).reverse_complement())

# GC content
def gc_content(seq: str) -> float:
    seq = seq.upper()
    return (seq.count("G") + seq.count("C")) / len(seq) * 100

# Tm 계산 (간단한 공식)
def calc_tm(seq: str) -> float:
    seq = seq.upper()
    return 2 * (seq.count("A") + seq.count("T")) + 4 * (seq.count("G") + seq.count("C"))

# NP ID로부터 NM ID 가져오기
def get_nm_from_np(np_id: str) -> str:
    url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/elink.fcgi"
    params = {
        "dbfrom": "protein", "db": "nuccore",
        "id": np_id, "linkname": "protein_nuccore_mrna",
        "retmode": "json"
    }
    r = requests.get(url, params=params)
    data = r.json()
    try:
        nm_id = data['linksets'][0]['linksetdbs'][0]['links'][0]
        return nm_id
    except:
        raise HTTPException(status_code=500, detail="NM ID를 찾을 수 없습니다.")

# CDS와 단백질 서열 가져오기
def fetch_cds_and_protein(np_id: str):
    nm_id = get_nm_from_np(np_id)
    handle = Entrez.efetch(db="nuccore", id=nm_id, rettype="gb", retmode="text")
    record = SeqIO.read(handle, "gb")
    handle.close()
    for feature in record.features:
        if feature.type == "CDS" and 'translation' in feature.qualifiers:
            cds_seq = feature.extract(record.seq)
            protein_seq = feature.qualifiers['translation'][0]
            return str(cds_seq), protein_seq
    raise HTTPException(status_code=500, detail="CDS를 찾을 수 없습니다.")

# 프라이머 디자인
def design_primer(cds: str, protein_seq: str, mutation: str, flank: int = 15):
    match = re.match(r"([A-Za-z])([0-9]+)([A-Za-z*])", mutation)
    if not match:
        raise ValueError("변이 형식이 잘못되었습니다 (예: S76A)")
    orig_aa, pos, new_aa = match.groups()
    pos = int(pos)

    if pos > len(protein_seq):
        raise ValueError("변이 위치가 단백질 길이를 초과합니다.")
    if protein_seq[pos - 1] != orig_aa:
        raise ValueError(f"{pos}번째 아미노산은 {protein_seq[pos - 1]}이며, {orig_aa}와 일치하지 않습니다.")

    codon_start = (pos - 1) * 3
    if codon_start + 3 > len(cds):
        raise ValueError("CDS에서 코돈 위치가 잘못되었습니다.")

    new_codon = PREFERRED_CODON[new_aa.upper()]
    upstream = cds[codon_start - flank:codon_start]
    downstream = cds[codon_start + 3:codon_start + 3 + flank]
    forward_primer = upstream + new_codon + downstream
    reverse_primer = reverse_complement(forward_primer)

    return {
        "forward_primer": forward_primer,
        "reverse_primer": reverse_primer,
        "tm": round(calc_tm(forward_primer), 1),
        "gc_percent": round(gc_content(forward_primer), 1)
    }

# API 엔드포인트
@app.get("/primer")
def primer_endpoint(
    refseq_protein: str = Query(..., min_length=6),
    mutation: str = Query(...)  # 정규식 제거
):
    try:
        cds, protein_seq = fetch_cds_and_protein(refseq_protein)
        primer_data = design_primer(cds, protein_seq, mutation)
    except Exception as e:
        raise HTTPException(status_code=500, detail=str(e))

    return {
        "refseq_protein": refseq_protein,
        "mutation": mutation.upper(),
        **primer_data
    }
