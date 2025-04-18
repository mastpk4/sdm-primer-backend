from fastapi import FastAPI, HTTPException, Query
from fastapi.middleware.cors import CORSMiddleware
import requests
from Bio import Entrez, SeqIO
from Bio.Seq import Seq
import re

app = FastAPI()

app.add_middleware(
    CORSMiddleware,
    allow_origins=["*"],
    allow_credentials=True,
    allow_methods=["*"],
    allow_headers=["*"],
)

Entrez.email = "your_email@example.com"

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

def get_nm_from_np(np_id: str) -> str:
    url = f"https://eutils.ncbi.nlm.nih.gov/entrez/eutils/elink.fcgi"
    params = {
        "dbfrom": "protein",
        "db": "nuccore",
        "id": np_id,
        "linkname": "protein_nuccore_mrna",
        "retmode": "json"
    }
    r = requests.get(url, params=params)
    data = r.json()
    try:
        return data['linksets'][0]['linksetdbs'][0]['links'][0]
    except:
        raise HTTPException(status_code=500, detail="NM ID not found for the given NP ID.")

def fetch_cds_and_protein(np_id: str):
    nm_id = get_nm_from_np(np_id)
    handle_prot = Entrez.efetch(db="protein", id=np_id, rettype="gb", retmode="text")
    prot_record = SeqIO.read(handle_prot, "gb")
    handle_prot.close()
    protein_name = prot_record.description

    handle = Entrez.efetch(db="nuccore", id=nm_id, rettype="gb", retmode="text")
    record = SeqIO.read(handle, "gb")
    handle.close()

    for feature in record.features:
        if feature.type == "CDS" and 'translation' in feature.qualifiers:
            cds_seq = feature.extract(record.seq)
            protein_seq = feature.qualifiers['translation'][0]
            return str(cds_seq), protein_seq, protein_name

    raise HTTPException(status_code=500, detail="CDS not found.")

def highlight_target(seq, codon, flank):
    return seq[:flank] + f"<b style='color:red'>{codon}</b>" + seq[flank+3:]

def design_primer(cds: str, protein_seq: str, mutation: str, flank: int = 15):
    match = re.match(r"([A-Za-z])([0-9]+)([A-Za-z*])", mutation)
    if not match:
        raise ValueError("Mutation format invalid (e.g., S76A)")

    orig_aa, pos, new_aa = match.groups()
    pos = int(pos)

    if pos > len(protein_seq):
        raise ValueError("Mutation position exceeds protein length.")

    if protein_seq[pos - 1] != orig_aa:
        raise ValueError(f"At position {pos}, expected {orig_aa} but found {protein_seq[pos - 1]}.")

    codon_start = (pos - 1) * 3
    new_codon = PREFERRED_CODON[new_aa.upper()]

    if codon_start + 3 > len(cds):
        raise ValueError("CDS codon location invalid.")

    up = cds[codon_start - flank:codon_start]
    down = cds[codon_start + 3:codon_start + 3 + flank]
    fwd = up + new_codon + down
    rev = reverse_complement(fwd)

    wt_codon = cds[codon_start:codon_start+3]
    wt_protein_context = protein_seq[pos-11:pos+10].replace(orig_aa, f"<b style='color:red'>{orig_aa}</b>", 1)
    wt_dna_context = cds[codon_start-30:codon_start+33].replace(wt_codon, f"<b style='color:red'>{wt_codon}</b>", 1)
    fwd_highlight = highlight_target(fwd, new_codon, flank)
    rev_highlight = highlight_target(rev, reverse_complement(new_codon), flank)

    return {
        "wt_protein_context": wt_protein_context,
        "wt_dna_context": wt_dna_context,
        "forward_primer": fwd_highlight,
        "reverse_primer": rev_highlight,
        "tm": round(calc_tm(fwd), 1),
        "gc_percent": round(gc_content(fwd), 1)
    }

@app.get("/primer")
def primer_endpoint(
    refseq_protein: str = Query(..., min_length=6),
    mutation: str = Query(..., pattern=r"^[A-Za-z]\d+[A-Za-z*]$")
):
    try:
        cds, protein_seq, protein_name = fetch_cds_and_protein(refseq_protein)
        primer = design_primer(cds, protein_seq, mutation)
    except Exception as e:
        raise HTTPException(status_code=500, detail=str(e))

    return {
        "refseq_protein": refseq_protein,
        "protein_name": protein_name,
        "mutation": mutation.upper(),
        **primer
    }
