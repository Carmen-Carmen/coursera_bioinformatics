import random
from enum import Enum

DNA_BASES = [
    "A", 
    "T", 
    "C", 
    "G", 
]

RNA_BASES = [
    "A", 
    "U", 
    "C", 
    "G", 
]

def generate_random_strand(length, type="DNA"):
    nucleotides = []
    if type == "DNA":
        nucleotides = DNA_BASES
    elif type == "RNA":
        nucleotides = RNA_BASES

    if len(nucleotides) == 0:
        return ""

    random_strand = "".join(
        random.choice(nucleotides) for _ in range(length)
    )

    return random_strand

class complement_direction(Enum):
    DNA_to_DNA = "D2D"
    DNA_to_RNA = "D2R"
    RNA_to_DNA = "R2D"
    RNA_to_RNA = "R2R"

DNA_PAIR_DNA = {
    "A": "T", 
    "T": "A", 
    "C": "G", 
    "G": "C", 
}

RNA_PAIR_RNA = {
    "A": "U", 
    "U": "A", 
    "C": "G", 
    "G": "C", 
}

DNA_PAIR_RNA = {
    "A": "U", 
    "T": "A", 
    "C": "G", 
    "G": "C", 
}

RNA_PAIR_DNA = {
    "A": "T", 
    "U": "A", 
    "C": "G", 
    "G": "C", 
}

def get_reverse_complement_strand(strand, direction=complement_direction.DNA_to_DNA):
    pairing = {}
    if direction == complement_direction.DNA_to_DNA:
        pairing = DNA_PAIR_DNA
    elif direction == complement_direction.DNA_to_RNA:
        pairing = DNA_PAIR_RNA
    elif direction == complement_direction.RNA_to_DNA:
        pairing = RNA_PAIR_DNA
    elif direction == complement_direction.RNA_to_RNA:
        pairing = RNA_PAIR_RNA
    else:
        raise ValueError("Invalid pairing direction")
    
    # complement
    complementary_strand = ""
    for nt in strand:
        complementary_strand += pairing[nt]
    
    # reverse
    reverse_complement_strand = "".join(
        reversed(complementary_strand)
    )

    return reverse_complement_strand
