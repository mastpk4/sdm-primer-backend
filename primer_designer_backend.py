from fastapi import FastAPI, HTTPException, Query
from fastapi.middleware.cors import CORSMiddleware
from Bio import Entrez, SeqIO
from Bio.Seq import Seq
import re
import requests
import os

app = FastAPI()

# CORS 설정
app.add_middleware(
    CORSMiddleware,
    allow_origins=["*"],
    allow_credentials=True,
    allow_methods=["*"],
    allow_headers=["*"],
)

# NCBI 이메일 설정
Entrez.email = "your_email@example.com"

# 아미노산 → 가장 흔한 codon (사람 기준)
PREFERRED_CODON = {
    "A": "GCC", "R": "CGT", "N": "AAC", "D": "GAC", "C": "TGC",
    "Q": "CAG", "E": "GAG", "G": "GGC", "H": "CAC", "I": "ATC",
    "L": "CTG", "K": "AAG", "M": "ATG", "F": "TTC", "P": "CCC",
    "S": "AGC", "T": "ACC", "W": "TGG", "Y": "TAC", "V": "GTG",
    "*": "TAA"
}

def fetch_cds_and_protein(np_id):
    # NP ID를 통해 NCBI에서 연결된 NM ID 획득
    handle = Entrez.efetch(db="protein", id=np_id, rettype="gb", retmode="text")
    record = SeqIO.read(handle, "genbank")
    handle.close()

    # CDS feature와 유전자 정보 파싱
    cds_feature = next((f for f in record.features if f.type == "CDS"), None)
    if not cds_feature:
        raise HTTPException(status_code=404, detail="CDS 정보를 찾을 수 없습니다.")

    # 유전자 서열 (GenBank에서 참조된 nucleotide accession 가져오기)
    protein_seq = str(record.seq.translate(to_stop=True))
    protein_id = cds_feature.qualifiers.get("protein_id", [None])[0]
    db_xref = cds_feature.qualifiers.get("db_xref", [""])[0]
    gene = cds_feature.qualifiers.get("gene", [""])[0]
    coded_by = cds_feature.qualifiers.get("coded_by", [""])[0]

    match = re.search(r"(NM_\d+\.\d+):?(\d+)?\.\.(\d+)?", coded_by)
    if not match:
        raise HTTPException(status_code=404, detail="CDS 위치 정보를 찾을 수 없습니다.")
    nm_id, start, end = match.groups()

    handle = Entrez.efetch(db="nucleotide", id=nm_id, rettype="fasta", retmode="text")
    record = SeqIO.read(handle, "fasta")
    handle.close()
    cds_seq = record.seq[int(start)-1:int(end)]

    return str(cds_seq), protein_seq

def reverse_complement(seq):
    return str(Seq(seq).reverse_complement())

def gc_content(seq):
    seq = seq.upper()
    return (seq.count("G") + seq.count("C")) / len(seq) * 100

def calc_tm(seq):
    seq = seq.upper()
    at = seq.count("A") + seq.count("T")
    gc = seq.count("G") + seq.count("C")
    return 2 * at + 4 * gc

def design_primer(cds, protein, mutation, flank=15):
    match = re.match(r"([A-Za-z])(\d+)([A-Za-z*])", mutation)
    if not match:
        raise ValueError("Mutation 형식은 예: S76A")

    orig_aa, pos, new_aa = match.groups()
    pos = int(pos)

    translated = str(Seq(cds).translate())
    if translated[pos - 1] != orig_aa:
        raise ValueError(f"{pos}\ubc88\uc9f8 \uc544\ubbf8\ub178\uc0ac\ub294 {translated[pos - 1]}\uc774\uba70, {orig_aa}\uc640 \uc77c\uce58\ud558\uc9c0 \uc54a\uc2b5\ub2c8\ub2e4.")

    codon_start = (pos - 1) * 3
    new_codon = PREFERRED_CODON[new_aa.upper()]
    up = cds[codon_start - flank: codon_start]
    down = cds[codon_start + 3: codon_start + 3 + flank]
    fwd = up + new_codon + down
    rev = reverse_complement(fwd)

    # context amino acid 표시용
    context = protein[max(0, pos - 11):pos - 1] + f"[{orig_aa}>{new_aa}]" + protein[pos:pos + 10]
    wt_dna = cds[codon_start - flank: codon_start + 3 + flank]

    return {
        "protein_context": context,
        "wt_dna": str(wt_dna),
        "forward_primer": fwd,
        "reverse_primer": rev,
        "tm": round(calc_tm(fwd), 1),
        "gc_percent": round(gc_content(fwd), 1)
    }

@app.get("/primer")
def primer_endpoint(
    refseq_protein: str = Query(..., min_length=6, max_length=15),
    mutation: str = Query(..., regex=r"^[A-Za-z]\d+[A-Za-z*]$")
):
    try:
        cds, protein = fetch_cds_and_protein(refseq_protein)
        primer = design_primer(cds, protein, mutation)
    except Exception as e:
        raise HTTPException(status_code=500, detail=str(e))
    return {
        "refseq_protein": refseq_protein,
        "mutation": mutation.upper(),
        **primer
    }
