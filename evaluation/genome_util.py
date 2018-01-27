import random
from collections import defaultdict


def read_genome_vmatch(genome_file):
    def process_seq(seq):
        seq = seq.upper()
        for base in 'RYKMSWBDHV':
            seq = seq.replace(base, 'N')
        return seq

    genome_seq = {}  # seq_count : (chrom_name, chrom_seq)
    chrom_seq = ""
    seq_count = 0  # Vmatch is zero-offset
    with open(genome_file) as ref_genome:
        chrom_name = next(ref_genome).rstrip().split(' ')[0][1:]
        for line in ref_genome:
            if line[0] == ">":
                if chrom_seq:
                    genome_seq[seq_count] = [chrom_name, list(chrom_seq)]
                    chrom_seq = ""
                    chrom_name = line.rstrip().split(' ')[0][1:]
                    seq_count += 1
            else:
                chrom_seq += process_seq(line.rstrip())
        genome_seq[seq_count] = [chrom_name, list(chrom_seq)]

    return genome_seq


def write_genome_vmatch(genome_seq, output_file):
    num_chrom = len(genome_seq)
    with open(output_file, 'w') as genome_file:
        for seq_num, seq in genome_seq.items():
            chrom_name = seq[0]
            chrom_seq = seq[1]
            # Sequence name
            genome_file.write('>' + chrom_name + '\n')
            # Writing sequence to file, 70 characters per line
            line_width = 70
            length = len(chrom_seq)
            for i in range(length // line_width):
                genome_file.write(''.join(chrom_seq[i * line_width: (i + 1) * line_width]))
                genome_file.write('\n')
            # Writing the last remainder part of genome
            if length % line_width != 0:
                genome_file.write(''.join(chrom_seq[-(length % line_width):]))
            # If there are more chromosomes to be written
            if num_chrom > 1:
                genome_file.write('\n')
                num_chrom -= 1


def reverse_complement(seq):
    complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'}
    return ''.join([complement[base] for base in seq[::-1]])


def read_vmatch_repeats(vmatch_repeats_file, min_len=150, d=0):
    """
    Returns a list of tuples where each tuple contains one repeat pair characteristics
    [(l, n1, r1, n2, r2), ...]
    """
    # Vmatch:
    # matches are reported in the following way
    # l(S) n(S) r(S) t l(S) n(S) r(S) d e s i
    # where:
    # l = length
    # n = sequence number
    # r = relative position
    # t = type (D=direct, P=palindromic)
    # d = distance value (negative=hamming distance, 0=exact, positive=edit distance)
    # e = E-value
    # s = score value (negative=hamming score, positive=edit score)
    # i = percent identity
    repeats = []
    with open(vmatch_repeats_file) as repeats_file:
        header = next(repeats_file)
        for line in repeats_file:
            features = line.rstrip().split()
            (rep_len, seq_num_1, pos_1, seq_num_2, pos_2, dist) = tuple(map(int, (features[0], features[1], features[2],
                                                                                  features[5], features[6],
                                                                                  features[7])))
            rep_type = features[3]
            if rep_len >= min_len and dist == d:
                repeats.append((rep_len, seq_num_1, pos_1, seq_num_2, pos_2, dist, rep_type))
    return repeats


def check_loc_in_repeats(genome_loc, all_repeats):
    """
    Checks if a mutation position is located in other repeats
    If this mutation is located within a in (a, b) repeat pair, returns b
    mut_loc: (seq_num, pos) tuple
    Returns a list as: [(seq_num, alt_pos)]
    """
    seq_num = genome_loc[0]
    pos = genome_loc[1]
    alt_locations = set()
    for repeat in all_repeats:
        (rep_len, seq_num_1, pos_1, seq_num_2, pos_2, dist, rep_type) = repeat
        if seq_num == seq_num_1 and pos_1 <= pos < pos_1 + rep_len:
            if rep_type == 'D':
                alt_pos = pos_2 + pos - pos_1
            else:
                # if it's palindromic repeat: alt_pos = end_of_repeat_pos - diff
                alt_pos = pos_2 + (rep_len - 1) - (pos - pos_1)
            alt_locations.add((seq_num_2, alt_pos, rep_type))
        elif seq_num == seq_num_2 and pos_2 <= pos < pos_2 + rep_len:
            if rep_type == 'D':
                alt_pos = pos_1 + pos - pos_2
            else:
                alt_pos = pos_1 + (rep_len - 1) - (pos - pos_2)
            alt_locations.add((seq_num_1, alt_pos, rep_type))
        else:
            pass
    return list(alt_locations)


def mutate_genome_repeats_beginning(ref_genome_file, main_repeats_file, all_repeats_file, output_file, mutations_file,
                                    read_len=150):
    """
    Mutating genome repeats obtained from vmatch -supermax repeat finding software
    """
    ref_genome = read_genome_vmatch(ref_genome_file)
    # Super-maximal repeats
    supermax_repeats = read_vmatch_repeats(main_repeats_file)
    # Tandem repeats
    tandem_repeats = read_vmatch_repeats(all_repeats_file)
    random.seed(12)
    nucleotides = {'A', 'C', 'G', 'T'}
    num_mutations = 0
    seen_repeats = set()
    with open(mutations_file, 'w') as mut_file:
        mut_file.write("#Chr\tPos\tRef\tAlt\tDist\tLen\n")
        for repeat in supermax_repeats:
            (rep_len, seq_num_1, pos_1, seq_num_2, pos_2, dist, rep_type) = repeat
            # If we have not mutated one of these repeats before
            if (seq_num_1, pos_1) not in seen_repeats and (seq_num_2, pos_2) not in seen_repeats:
                # Position of mutation from the start of repeat element
                mut_pos = random.choice(range(read_len - 50, read_len - 9))
                # If it's a short range mutation, check that it's not located in tandem repeats
                # Checking the position right before repeat starting point
                loc1_in_repeat = check_loc_in_repeats((seq_num_1, pos_1 - 1), tandem_repeats)
                loc2_in_repeat = check_loc_in_repeats((seq_num_2, pos_2 - 1), tandem_repeats)
                # Which repeat instance to mutate
                if not loc1_in_repeat and not loc2_in_repeat:
                    rep_loc = random.choice([(seq_num_1, pos_1), (seq_num_2, pos_2)])
                elif not loc1_in_repeat:
                    rep_loc = (seq_num_1, pos_1)
                elif not loc2_in_repeat:
                    rep_loc = (seq_num_2, pos_2)
                else:
                    continue
                # Absolute position of mutation on the chromosome
                chr_pos = rep_loc[1] + mut_pos
                # The original reference base
                # seq_num -> chrom_seq -> base
                ref_base = ref_genome[rep_loc[0]][1][chr_pos]
                if ref_base == 'N':  # undetermined base in reference genome
                    continue
                # Mutating the base
                possible_snps = nucleotides - set(ref_base)
                new_base = random.choice(sorted(list(possible_snps)))
                ref_genome[rep_loc[0]][1][chr_pos] = new_base
                # Write mutation to file
                # Position is incremented since position in VCF files starts from one (not zero)
                mut_file.write("{}\t{}\t{}\t{}\t{}\t{}\n".format(ref_genome[rep_loc[0]][0], chr_pos + 1, ref_base,
                                                                 new_base, mut_pos, rep_len))
                num_mutations += 1
                # Mark these repeats as seen
                seen_repeats.add((seq_num_1, pos_1))
                seen_repeats.add((seq_num_2, pos_2))
        mut_file.write("#Number of mutations: {}".format(num_mutations))
    # Write mutated genome to file
    write_genome_vmatch(ref_genome, output_file)
    return True


def mutate_genome_repeats_middle(ref_genome_file, supermax_repeats_file, all_repeats_file, output_file, mutations_file,
                                 read_len=150):
    """
    Mutating genome repeats obtained from vmatch -supermax repeat finding software
    """
    ref_genome = read_genome_vmatch(ref_genome_file)
    # Super-maximal repeats
    supermax_repeats = read_vmatch_repeats(supermax_repeats_file)
    # Tandem repeats
    all_repeats = read_vmatch_repeats(all_repeats_file)
    random.seed(12)
    nucleotides = {'A', 'C', 'G', 'T'}
    num_mutations = 0
    seen_repeats = set()
    with open(mutations_file, 'w') as mut_file:
        mut_file.write("#Chr\tPos\tRef\tAlt\tDist\tLen\tOtherLocs\n")
        for repeat in supermax_repeats:
            (rep_len, seq_num_1, pos_1, seq_num_2, pos_2, dist, rep_type) = repeat
            # If this repeat is at least twice as longer as a read length
            # we can find a position that is out of reach of any unique reads
            if rep_len > 2 * read_len:
                # If we have not mutated one of these repeats before
                if (seq_num_1, pos_1) not in seen_repeats and (seq_num_2, pos_2) not in seen_repeats:

                    # Position of mutation in the middle of the repeat element
                    mut_pos = read_len + random.choice(range(rep_len - 2 * read_len))
                    # Which repeat instance to mutate
                    rep_loc = random.choice([(seq_num_1, pos_1), (seq_num_2, pos_2)])
                    # Absolute position of mutation on the chromosome
                    chr_pos = rep_loc[1] + mut_pos

                    # Find the corresponding repeat locations
                    other_repeats = check_loc_in_repeats((rep_loc[0], chr_pos), all_repeats)

                    # The original reference base
                    ref_base = ref_genome[rep_loc[0]][1][chr_pos]
                    if ref_base == 'N':  # undetermined base in reference genome
                        continue
                    # Mutating the base
                    possible_snps = nucleotides - set(ref_base)
                    new_base = random.choice(sorted(list(possible_snps)))
                    ref_genome[rep_loc[0]][1][chr_pos] = new_base

                    # Finding the acceptable false positives
                    acceptable_fps = []
                    for other_repeat_loc in other_repeats:
                        (seq_num_other, pos_other, rep_type_other) = other_repeat_loc
                        if rep_type_other == 'D':
                            alt_base = new_base
                        else:
                            alt_base = reverse_complement(new_base)
                        chrom_name_alt = ref_genome[seq_num_other][0]
                        acceptable_fps.append((chrom_name_alt, pos_other, alt_base))

                    # Write mutation to file
                    # Position is incremented since position in VCF files starts from one (not zero)
                    mut_file.write("{}\t{}\t{}\t{}\t{}\t{}\t{}\n".format(ref_genome[rep_loc[0]][0], chr_pos + 1,
                                                                         ref_base, new_base, mut_pos, rep_len,
                                                                         acceptable_fps))
                    num_mutations += 1
                    # Mark these repeats as seen
                    seen_repeats.add((seq_num_1, pos_1))
                    seen_repeats.add((seq_num_2, pos_2))
        mut_file.write("#Number of mutations: {}".format(num_mutations))
    # Write mutated genome to file
    write_genome_vmatch(ref_genome, output_file)
    return True


def back_mutate_genome_repeats(ref_genome_file, supermax_repeats_file, all_repeats_file, output_file, mutations_file,
                       read_len=150):
    """
    Mutating genome repeats obtained from vmatch -supermax repeat finding software
    """
    ref_genome = read_genome_vmatch(ref_genome_file)
    # Super-maximal repeats with hamming distance 1
    supermax_repeats = read_vmatch_repeats(supermax_repeats_file, read_len, d=-1)
    # All repeats
    all_repeats = read_vmatch_repeats(all_repeats_file)
    random.seed(12)
    nucleotides = {'A', 'C', 'G', 'T'}
    num_mutations = 0
    seen_repeats = set()
    with open(mutations_file, 'w') as mut_file:
        mut_file.write("#Chr\tPos\tRef\tAlt\tDist\tLen\tOtherLocs\n")
        for repeat in supermax_repeats:
            (rep_len, seq_num_1, pos_1, seq_num_2, pos_2, dist, rep_type) = repeat

            # We should find the point of difference
            seq_1 = ref_genome[seq_num_1][1][pos_1:pos_1+rep_len]
            seq_2 = ref_genome[seq_num_2][1][pos_2:pos_2+rep_len]
            mut_seq_num = -1
            mut_chr_pos = -1
            for i in range(rep_len):
                if seq_1[i] != seq_2[i]:
                    # Check whether the mutation point is not located itself in an exact repeat element
                    other_repeats_1 = check_loc_in_repeats((seq_num_1, pos_1 + i), all_repeats)
                    other_repeats_2 = check_loc_in_repeats((seq_num_2, pos_2 + i), all_repeats)
                    if not other_repeats_1:
                        mut_seq_num = seq_num_1
                        mut_chr_pos = pos_1 + i
                    elif not other_repeats_2:
                        mut_seq_num = seq_num_2
                        mut_chr_pos = pos_2 + i
                    break

            # If we have found a position that is not located inside an exact repeat and
            # we have not mutated one of these repeats before
            if mut_seq_num >= 0 and (seq_num_1, pos_1) not in seen_repeats and (seq_num_2, pos_2) not in seen_repeats:
                # its distance from start or end of repeat is longer than read_len - 50 (eg > 100 bp for a 150 bp read)
                # and (i in range(read_len - 50, rep_len - (read_len - 50))):
                # The original reference base
                ref_base = ref_genome[mut_seq_num][1][mut_chr_pos]
                if ref_base == 'N':  # undetermined base in reference genome
                    continue
                # Mutating the base
                possible_snps = nucleotides - set(ref_base)
                new_base = random.choice(sorted(list(possible_snps)))
                ref_genome[mut_seq_num][1][mut_chr_pos] = new_base

                # Write mutation to file
                # Order: back-mutated nucleotide, true nucleotide
                # Position is incremented since position in VCF files starts from one (not zero)
                mut_file.write("{}\t{}\t{}\t{}\t{}\t{}\n".format(ref_genome[mut_seq_num][0], mut_chr_pos + 1,
                                                                     new_base, ref_base, i, rep_len))
                num_mutations += 1
                # Mark these repeats as seen
                seen_repeats.add((seq_num_1, pos_1))
                seen_repeats.add((seq_num_2, pos_2))
        mut_file.write("#Number of mutations: {}".format(num_mutations))
    # Write mutated genome to file
    write_genome_vmatch(ref_genome, output_file)
    return True
