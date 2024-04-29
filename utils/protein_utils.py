# util functions for problems concerning proteins, peptides and amino acids

RNA_codon_to_AA = {
    "AAA": "K", "AAC": "N", "AAG": "K", "AAU": "N", "ACA": "T", "ACC": "T", "ACG": "T", "ACU": "T",
    "AGA": "R", "AGC": "S", "AGG": "R", "AGU": "S", "AUA": "I", "AUC": "I", "AUG": "M", "AUU": "I",
    "CAA": "Q", "CAC": "H", "CAG": "Q", "CAU": "H", "CCA": "P", "CCC": "P", "CCG": "P", "CCU": "P",
    "CGA": "R", "CGC": "R", "CGG": "R", "CGU": "R", "CUA": "L", "CUC": "L", "CUG": "L", "CUU": "L",
    "GAA": "E", "GAC": "D", "GAG": "E", "GAU": "D", "GCA": "A", "GCC": "A", "GCG": "A", "GCU": "A",
    "GGA": "G", "GGC": "G", "GGG": "G", "GGU": "G", "GUA": "V", "GUC": "V", "GUG": "V", "GUU": "V",
    "UAA": "*", "UAC": "Y", "UAG": "*", "UAU": "Y", "UCA": "S", "UCC": "S", "UCG": "S", "UCU": "S",
    "UGA": "*", "UGC": "C", "UGG": "W", "UGU": "C", "UUA": "L", "UUC": "F", "UUG": "L", "UUU": "F",
}

AA_3_letter_to_abbr = {
    "Gly": "G", "Ala": "A", "Val": "V", "Leu": "L", 
    "Gly": "G", "Ala": "A", "Val": "V", "Leu": "L", 
    "Ile": "I", "Pro": "P", "Phe": "F", "Tyr": "Y", 
    "Trp": "W", "Ser": "S", "Thr": "T", "Cys": "C", 
    "Met": "M", "Asn": "N", "Gln": "Q", "Asp": "D", 
    "Glu": "E", "Lys": "K", "Arg": "R", "His": "H"
}

AA_abbr_to_3_letter = {
    "G": "Gly", "A": "Ala", "V": "Val", "L": "Leu", 
    "I": "Ile", "P": "Pro", "F": "Phe", "Y": "Tyr", 
    "W": "Trp", "S": "Ser", "T": "Thr", "C": "Cys", 
    "M": "Met", "N": "Asn", "Q": "Gln", "D": "Asp", 
    "E": "Glu", "K": "Lys", "R": "Arg", "H": "His",
}

AA_to_RNA_codon = {
    "*": ['UAA', 'UAG', 'UGA'], 
    "A": ['GCA', 'GCC', 'GCG', 'GCU'], 
    "C": ['UGC', 'UGU'], 
    "D": ['GAC', 'GAU'], 
    "E": ['GAA', 'GAG'], 
    "F": ['UUC', 'UUU'], 
    "G": ['GGA', 'GGC', 'GGG', 'GGU'], 
    "H": ['CAC', 'CAU'], 
    "I": ['AUA', 'AUC', 'AUU'], 
    "K": ['AAA', 'AAG'], 
    "L": ['CUA', 'CUC', 'CUG', 'CUU', 'UUA', 'UUG'], 
    "M": ['AUG'], 
    "N": ['AAC', 'AAU'], 
    "P": ['CCA', 'CCC', 'CCG', 'CCU'], 
    "Q": ['CAA', 'CAG'], 
    "R": ['AGA', 'AGG', 'CGA', 'CGC', 'CGG', 'CGU'], 
    "S": ['AGC', 'AGU', 'UCA', 'UCC', 'UCG', 'UCU'], 
    "T": ['ACA', 'ACC', 'ACG', 'ACU'], 
    "V": ['GUA', 'GUC', 'GUG', 'GUU'], 
    "W": ['UGG'], 
    "Y": ['UAC', 'UAU'],
}

AA_weight = {
    # "G": 57, 
    # "A": 71, 
    # "S": 87, 
    # "P": 97, 
    # "V": 99, 
    # "T": 101, 
    # "C": 103, 
    # "I": 113, 
    # "L": 113, 
    # "N": 114, 
    # "D": 115, 
    # "K": 128, 
    # "Q": 128, 
    # "E": 129, 
    # "M": 131, 
    # "H": 137, 
    # "F": 147, 
    # "R": 156, 
    # "Y": 163, 
    # "W": 186, 
    # below are fake AAs: 
    "X": 4, "Z": 5
}

weight_AA = {
    # 57: 'G', 71: 'A', 87: 'S', 97: 'P',
    # 99: 'V', 101: 'T', 103: 'C', 113:'I/L',
    # 114: 'N', 115: 'D', 128: 'K/Q', 129: 'E',
    # 131: 'M', 137: 'H', 147: 'F', 156: 'R', 
    # 163: 'Y', 186: 'W', 
    # below are fake AAs: 
    4: "X", 5: "Z"
}

# ALL_AA_WEIGHTS = list(weight_AA.keys()) + [113, 128]
ALL_AA_WEIGHTS = list(weight_AA.keys())

def translate_RNA_to_peptide(RNA_strand):
    peptide_string = ""
    i = 0
    while i < len(RNA_strand):
        codon = RNA_strand[i : i + 3]
        amino_acid = RNA_codon_to_AA[codon]
        if amino_acid == "*":
            break
        peptide_string += amino_acid

        i += 3

    return peptide_string

def get_AA_abbr(AA_3_letter):
    return AA_3_letter_to_abbr[AA_3_letter]

def get_AA_3_letter(AA_abbr):
    return AA_abbr_to_3_letter[AA_abbr]

def get_weight_by_AA(AA):
    return AA_weight[AA]

def get_AA_by_weight(weight: int): 
    AA = weight_AA[weight] 
    if "/" in AA: 
        return AA.split("/")[1]
    else:
        return AA

# all AAs in the given peptide are expressed in 1 single letter
def get_weight_by_peptide(peptide: str): 
    weight = 0
    for AA in peptide: 
        weight += get_weight_by_AA(AA)
    
    return weight
