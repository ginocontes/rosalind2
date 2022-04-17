from collections import defaultdict

def dictionize_nucleotides(nucleotides_string):
    dd = defaultdict(int)
    for nucleotide in nucleotides_string:
        dd[nucleotide] += 1
    return dd


def transcribe_rna(dna_string: str):
    return dna_string.replace('T', 'U')


def reverse_complement(dna_string: str):
    complements = {"A": "T", "C": "G", "T": "A", "G": "C"}
    res = ""
    for nt in reversed(dna_string):
        res += complements[nt]
    return res

def rabbits(n, k):
    f1 = 0
    f2 = 1
    for i in range(1,n):
        temp = f1
        f1 = f2
        f2 = f2 + k*temp
    return f2


def parse_fasta(lines: list[str]):
    res = defaultdict(str)
    for line in lines:
        if line.startswith('>'):
            last_label = line[1:].rstrip()
        else:
            res[last_label] += line.rstrip()
    return res

def parse_fasta_file(path):
    with open(path, "r") as f:
        lines = f.readlines()

    return parse_fasta(lines)


def max_gc_content(path):
    strands = parse_fasta_file(path)
    max_label = ""
    max_content = 0
    for k,v in strands.items():
        gc_content = calculate_GC_Content(v)
        if gc_content > max_content:
            max_content = gc_content
            max_label = k
    return max_label, max_content

def calculate_GC_Content(strand):
    n = len(strand)
    c = 0
    for el in strand:
        if el in 'GC':
            c += 1
    return c/n



def hamming_distance(s1, s2):
    n = min(len(s1), len(s2))
    hd = 0
    for i in range(n):
        if s1[i] != s2[i]:
            hd += 1
    return hd


def mendel_laws(homozygous_dominant, heterozygous, homozygous_recessives):
    """
    :return the probability that two randomly seclected from the population will produce an
    individual possessing a dominant allele(and thus displaying the dominant phenotype).
    :param homozygous_dominant:
    :param heterozygous:
    :param homozygous_recessives:
    :return:
    """
    n = homozygous_dominant + homozygous_recessives + heterozygous
    return 1 - ((((heterozygous * homozygous_recessives)/2) + ((heterozygous*heterozygous-1)/2)) /(n*(n-1))) #out of all the possibile combinations


def parse_codon_table(filename):
    res = dict()
    with open(filename, "r") as f:
        lines = f.readlines()
        for line in lines:
            line = line.rstrip()
            line = list(filter(lambda x: x!="", line.split(" ")))
            for i in range(0, len(line), 2):
                res[line[i]] = line[i+1]
    return res




def rna_translation(s, filename):
    codon_table = parse_codon_table(filename)
    protein_string = ""
    for i in range(0, len(s), 3):
        codon = s[i:i+3]
        amino_acid = codon_table[codon]
        if amino_acid == "Stop":
            return protein_string
        protein_string += amino_acid

    return protein_string