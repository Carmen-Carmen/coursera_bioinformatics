
def compute_Nxx_by_contig_lengths(Nxx, contig_lengths):
    total_length = sum(contig_lengths)
    target_proportion = Nxx / 100

    contig_lengths.sort()
    contig_lengths = reversed(contig_lengths)

    added_len = 0
    for length in contig_lengths:
        added_len += length
        temp_proportion = added_len / total_length

        if temp_proportion >= target_proportion:
            return length

def compute_NGxx_by_contig_lengths(NGxx, contig_lengths, genome_length):
    target_proportion = NGxx / 100

    contig_lengths.sort()
    contig_lengths = reversed(contig_lengths)

    added_len = 0
    for length in contig_lengths:
        added_len += length
        temp_proportion = added_len / genome_length

        if temp_proportion >= target_proportion:
            return length
        
if __name__ == "__main__":
    contig_lengths = [20, 20, 30, 30, 60, 60, 80, 100, 200]
    # computing N50 and N75
    print("N50:", compute_Nxx_by_contig_lengths(50, contig_lengths))
    print("N75:", compute_Nxx_by_contig_lengths(75, contig_lengths))

    # computing NG50
    genome_length = 1000
    print("NG50:", compute_NGxx_by_contig_lengths(50, contig_lengths, genome_length))

    # computing NGA50
    genome_length = 1000
    # the contig in our dataset of length 100 had a misassembly breakpoint in the middle of it
    # 100 --> 50, 50
    contig_lengths_broken = [20, 20, 30, 30, 60, 60, 80, 50, 50, 200]
    print("NGA50:", compute_NGxx_by_contig_lengths(50, contig_lengths_broken, genome_length))