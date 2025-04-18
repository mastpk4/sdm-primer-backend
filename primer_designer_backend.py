# primer_designer_backend.py
from fastapi import FastAPI, HTTPException, Query
from fastapi.middleware.cors import CORSMiddleware
from Bio import Entrez, SeqIO, Seq
import re

app = FastAPI()

# CORS 허용
app.add_middleware(
    CORSMiddleware,
    allow_origins=["*"], allow_credentials=True,
    allow_methods=["*"], allow_headers=["*"],
)

Entrez.email = "your_email@example.com"  # 실제 메일 주소로 변경 추천

PREFERRED_CODON = {
    "A": "GCC", "R": "CGT", "N": "AAC", "D": "GAC", "C": "TGC",
    "Q": "CAG", "E": "GAG", "G": "GGC", "H": "CAC", "I": "ATC",
    "L": "CTG", "K": "AAG", "M": "ATG", "F": "TTC", "P": "CCC",
    "S": "AGC", "T": "ACC", "W": "TGG", "Y": "TAC", "V": "GTG",
    "*": "TAA"
}

def reverse_complement(seq: str) -> str:
    return str(Seq.Seq(seq).reverse_complement())

def calc_tm(seq: str) -> float:
    seq = seq.upper()
    return 2 * (seq.count("A") + seq.count("T")) + 4 * (seq.count("G") + seq.count("C"))

def fetch_protein_sequence(np_id: str) -> str:
    handle = Entrez.efetch(db="protein", id=np_id, rettype="fasta", retmode="text")
    record = SeqIO.read(handle, "fasta")
    return str(record.seq)

def fetch_nucleotide_cds(np_id: str) -> tuple[str, str]:
    # NP → NM 매핑
    summary = Entrez.read(Entrez.esummary(db="protein", id=np_id))
    linked_nuc = summary[0]["Caption"].replace("NP_", "NM_")  # crude mapping
    # Fetch NM record
    handle = Entrez.efetch(db="nucleotide", id=linked_nuc, rettype="gb", retmode="text")
    record = SeqIO.read(handle, "genbank")
    cds_feature = next(f for f in record.features if f.type == "CDS")
    cds_seq = cds_feature.extract(record.seq)
    return str(cds_seq), linked_nuc

def design_primer(cds: str, protein: str, mutation: str, flank: int = 15):
    match = re.match(r"([A-Za-z])(\d+)([A-Za-z*])", mutation)
    if not match:
        raise ValueError("Mutation 형식 오류 (예: S76A)")
    orig_aa, pos, new_aa = match.groups()
    pos = int(pos)
    translated = str(Seq.Seq(cds).translate())

    # CDS에서 번역한 단백질이 NP 시퀀스의 일부인지 확인
    if translated.find(protein[:20]) == -1 and protein.find(translated[:20]) == -1:
        raise ValueError("CDS에서 번역한 단백질이 UniProt과 일치하지 않습니다.")

    if translated[pos - 1] != orig_aa:
        raise ValueError(f"{pos}번째 아미노산은 {translated[pos - 1]}이며, {orig_aa}와 일치하지 않습니다.")

    codon_start = (pos - 1) * 3
    up = cds[max(0, codon_start - flank):codon_start]
    down = cds[codon_start + 3: codon_start + 3 + flank]
    mutant_codon = PREFERRED_CODON[new_aa.upper()]
    fwd = up + mutant_codon.lower() + down
    rev = reverse_complement(fwd)
    wt_dna = up + cds[codon_start:codon_start + 3].lower() + down

    protein_context = translated[max(0, pos - 11): pos - 1] + "[" + orig_aa + "→" + new_aa + "]" + translated[pos: pos + 10]

    return {
        "protein_context": protein_context,
        "wt_dna": wt_dna,
        "forward_primer": fwd,
        "reverse_primer": rev,
        "tm": round(calc_tm(fwd), 1),
        "gc_percent": round((fwd.upper().count("G") + fwd.upper().count("C")) / len(fwd) * 100, 1)
    }

@app.get("/primer")
def primer_endpoint(
    refseq_protein: str = Query(...),
    mutation: str = Query(..., regex=r"^[A-Za-z]\d+[A-Za-z*]$")
):
    try:
        protein = fetch_protein_sequence(refseq_protein)
        cds, nm_id = fetch_nucleotide_cds(refseq_protein)
        primer = design_primer(cds, protein, mutation)
    except Exception as e:
        raise HTTPException(status_code=500, detail=str(e))
    return {
        "refseq_protein": refseq_protein,
        "mutation": mutation.upper(),
        **primer
    }
