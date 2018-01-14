import random
from collections import defaultdict


def read_genome_vmatch(genome_file):
    genome_seq = {}  # seq_count : (chrom_name, chrom_seq)
    chrom_seq = ""
    seq_count = 0  # Vmatch is zero-offset
    with open(genome_file) as ref_genome:
        chrom_name = next(ref_genome).split(' ')[0][1:]
        for line in ref_genome:
            if line[0] == ">":
                if chrom_seq:
                    genome_seq[seq_count] = [chrom_name, list(chrom_seq)]
                    chrom_seq = ""
                    chrom_name = line.split(' ')[0][1:]
                    seq_count += 1
            else:
                chrom_seq += line.rstrip()
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


def read_vmatch_repeats(vmatch_repeats_file):
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
            (rep_len, seq_num_1, pos_1, seq_num_2, pos_2) = tuple(map(int, (features[0], features[1], features[2],
                                                                            features[5], features[6])))
            repeats.append((rep_len, seq_num_1, pos_1, seq_num_2, pos_2))
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
    alt_locations = []
    for repeat in all_repeats:
        (rep_len, seq_num_1, pos_1, seq_num_2, pos_2) = repeat
        if seq_num == seq_num_1 and pos_1 <= pos < pos_1 + rep_len:
            alt_locations.append((seq_num_2, pos_2 + pos - pos_1))
        elif seq_num == seq_num_2 and pos_2 <= pos < pos_2 + rep_len:
            alt_locations.append((seq_num_1, pos_1 + pos - pos_2))
        else:
            pass
    return alt_locations


def mutate_genome_repeats(ref_genome_file, main_repeats_file, all_repeats_file, output_file, mutations_file,
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
            (rep_len, seq_num_1, pos_1, seq_num_2, pos_2) = repeat
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
                # If it's a short range mutation, check that it's not located in tandem repeats
                # alt_repeats = check_loc_in_repeats((rep_loc[0], rep_loc[1]-1), tandem_repeats)
                # if alt_repeats:
                #     continue
                # The original reference base
                ref_base = ref_genome[rep_loc[0]][1][chr_pos]
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
        mut_file.write("#Chr\tPos\tRef\tAlt\tDist\tLen\nOtherLocs\n")
        for repeat in supermax_repeats:
            (rep_len, seq_num_1, pos_1, seq_num_2, pos_2) = repeat
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
                    # If it's a short range mutation, check that it's not located in tandem repeats
                    alt_repeats = check_loc_in_repeats((rep_loc[0], chr_pos), all_repeats)
                    if alt_repeats:
                        continue
                    # The original reference base
                    ref_base = ref_genome[rep_loc[0]][1][chr_pos]
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


def read_genome(genome_file):
    genome_header = ""
    genome_seq = ""
    with open(genome_file) as ref_genome:
        for line in ref_genome:
            # Skip header lines
            if line[0] == ">":
                genome_header += line
            else:
                genome_seq += line.rstrip()
    return genome_header, genome_seq


def write_genome(genome_header, genome_seq, genome_file):
    with open(genome_file, 'w') as new_genome:
        new_genome.write(genome_header)
        # Writing the new genome sequence to file, 70 characters per line
        line_width = 70
        length = len(genome_seq)
        for i in range(length // line_width):
            new_genome.write(genome_seq[i * line_width: (i + 1) * line_width])
            new_genome.write("\n")
        # Writing the last remainder part of genome
        if length % line_width != 0:
            new_genome.write(genome_seq[-(length % line_width):])


def reverse_complement(seq):
    complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'}
    return ''.join([complement[base] for base in seq[::-1]])


def genome_slice(genome_file, locations):
    genome_header, genome_seq = read_genome(genome_file)
    for loc in locations:
        # start_pos, length, direction
        if loc[2] == 'P':
            print(reverse_complement(genome_seq[loc[0]: loc[0] + loc[1]]))
        else:
            print(genome_seq[loc[0]: loc[0] + loc[1]])


def extract_genome(ref_genome_file, start_pos, length, output_file, mutate=False, repeats_file_name=None):
    """
    Gets a Fasta file, extracts a part of genome, and writes it back to a new Fasta file as a new reference genome.
    """
    genome_seq = ""
    with open(ref_genome_file) as ref_genome:
        with open(output_file, 'w') as new_genome:
            for line in ref_genome:
                # Skip header lines
                if line[0] == ">":
                    # header_line = line.rstrip().split("|")
                    # new_genome.write("{} {}-{}\n".format("|".join(header_line), start_pos, start_pos + length - 1))
                    # new_genome.write("{}|{}_{}|\n".format("|".join(header_line[:4]),start_pos, start_pos + length - 1))
                    new_genome.write(line)
                else:
                    genome_seq += line.rstrip()

            # Extracted genome sequence
            new_genome_seq = genome_seq[start_pos - 1: start_pos + length - 1]

            if mutate:
                new_genome_seq = mutate_genome(repeats_file_name, new_genome_seq, start_pos, length,
                                               "{}-mutations.txt".format(output_file[:-4]))

            # Writing the new genome sequence to file, 70 characters per line
            line_width = 70
            for i in range(length // line_width):
                new_genome.write(new_genome_seq[i * line_width: (i + 1) * line_width])
                new_genome.write("\n")
            # Writing the last remainder part of genome
            if length % line_width != 0:
                new_genome.write(new_genome_seq[-(length % line_width):])
    return True


def mutate_genome(repeats_file_name, genome_sequence, start_pos, genome_length, mutation_locations_file):
    random.seed(12)

    # [[length, start_1, start_2, ...], ...]
    repeats_list = []
    # [range(start, start + length), ...]
    repeats_ranges = set()
    with open(repeats_file_name) as repeats_file:
        for line in repeats_file:
            fields = line.split()
            # For repeats on reverse strand, remove 'r' from end position
            if fields[1][-1] == 'r':
                fields[1] = fields[1][:-1]
            # Deducting start_pos to make absolute positions for the new genome extract
            (start, end, length) = (int(fields[0]) - start_pos, int(fields[1]) - start_pos, int(fields[2]))
            # Filtering repeats that are smaller than read length
            # if length > 150:
            # Some repeat coordinates may exceed the extract of the genome we are looking for
            if start + length < genome_length and end + length < genome_length:
                if repeats_list and repeats_list[-1][0] == length and repeats_list[-1][1] == start:
                    repeats_list[-1].append(end)
                else:
                    repeats_list.append([length, start, end])

                # Adding this range to ranges of repeat positions
                repeats_ranges.add(range(start, start + length))
                repeats_ranges.add(range(end, end + length))

    # Saving repeat ranges to file
    ranges_list = []
    for rng in repeats_ranges:
        ranges_list.append((list(rng)[0], list(rng)[-1], list(rng)[-1] - list(rng)[0] + 1))
    ranges_list.sort()
    with open("{}-ranges.txt".format(repeats_file_name[:-4]), "w") as repeats_ranges_file:
        for rng in ranges_list:
            repeats_ranges_file.write("{}\t{}\t{}\n".format(rng[0], rng[1], rng[2]))

    # Merging repeat ranges and writing them to file
    merged_ranges = []
    for begin, end, rng_len in ranges_list:
        if merged_ranges and merged_ranges[-1][1] >= begin - 1:
            merged_ranges[-1][1] = max(merged_ranges[-1][1], end)
            # Update repeat length
            merged_ranges[-1][2] = merged_ranges[-1][1] - merged_ranges[-1][0] + 1
        else:
            merged_ranges.append([begin, end, end - begin + 1])

    with open("{}-ranges-merged.txt".format(repeats_file_name[:-4]), "w") as repeats_ranges_file:
        for rng in merged_ranges:
            repeats_ranges_file.write("{}\t{}\t{}\n".format(rng[0], rng[1], rng[2]))

    long_repeats = []
    for rng in merged_ranges:
        if rng[2] > 150:
            long_repeats.append(rng)

    # Mutating genome
    new_genome_seq = list(genome_sequence)
    mutations_list = []

    # 25% of regions
    # num_snps = len(merged_ranges) // 4
    num_snps = len(long_repeats) // 2
    # selected_regions = random.sample(merged_ranges, num_snps)
    selected_regions = random.sample(long_repeats, num_snps)
    nucleotides = set(['A', 'C', 'G', 'T'])
    for region in selected_regions:
        mapping_location = region[0]
        base_position = random.choice(range(100, 140))
        # base_position = 140
        final_position = mapping_location + base_position
        print(region, mapping_location, base_position, final_position)
        # Mutating the base
        possible_snps = nucleotides - set(genome_sequence[final_position])
        new_genome_seq[final_position] = random.choice(sorted(list(possible_snps)))
        # Saving this mutation characteristics
        mutations_list.append((final_position + 1, genome_sequence[final_position], new_genome_seq[final_position],
                               base_position, region))

    # Writing mutation locations to file
    mutations_list.sort()
    with open(mutation_locations_file, "w") as mutations_file:
        mutations_file.write("Number of mutations: {}\n".format(num_snps))
        for item in mutations_list:
            mutations_file.write("{}\t{}\t{}\t{}\t{}\n".format(item[0], item[1], item[2], item[3], item[4]))

    # Returning mutated genome sequence
    return "".join(new_genome_seq)


def toy_genome(ref_genome_file, output_file, mutate=False):
    """
    Gets a Fasta file, extracts a part of genome, and writes it back to a new Fasta file as a new reference genome.
    """
    random.seed(12)
    genome_seq = ""

    with open(output_file, 'w') as new_genome:
        # Reading FASTA genome file
        genome_header, genome_seq = read_genome(ref_genome_file)
        new_genome.write(genome_header)

        # Extracted genome sequence
        new_genome_seq = genome_seq[0:5000] + genome_seq[2000:3000] + genome_seq[6000:8000] \
                         + genome_seq[2000:3000] + genome_seq[9000:11000]

        new_genome_seq = list(new_genome_seq)

        # Mutating the genome
        if mutate:
            nucleotides = set(['A', 'C', 'G', 'T'])
            pos = 2129
            # pos = 2500
            print(new_genome_seq[pos])
            possible_snps = nucleotides - set(new_genome_seq[pos])
            new_genome_seq[pos] = random.choice(sorted(list(possible_snps)))
            print(new_genome_seq[pos])

        new_genome_seq = "".join(new_genome_seq)
        # Writing the new genome sequence to file, 70 characters per line
        length = 11000
        line_width = 70
        for i in range(length // line_width):
            new_genome.write(new_genome_seq[i * line_width: (i + 1) * line_width])
            new_genome.write("\n")
        # Writing the last remainder part of genome
        if length % line_width != 0:
            new_genome.write(new_genome_seq[-(length % line_width):])
    return True


def merge_repeat_ranges(repeats_file_name):
    # [[length, start_1, start_2, ...], ...]
    repeats_list = []
    # [range(start, start + length), ...]
    repeats_ranges = set()
    with open(repeats_file_name) as repeats_file:
        # Skip header lines
        next(repeats_file)
        next(repeats_file)
        for line in repeats_file:
            fields = line.split()
            # For repeats on reverse strand, remove 'r' from end position
            # if fields[1][-1] == 'r':
            #     fields[1] = fields[1][:-1]
            # Deducting start_pos to make absolute positions for the new genome extract
            # Mummer
            # (start, end, length) = (int(fields[0]), int(fields[1]), int(fields[2]))
            # Reputer
            if fields[5] == '-2':
                continue
            (start, end, length) = (int(fields[1]), int(fields[4]), int(fields[0]))
            # Filtering repeats that are smaller than read length
            # if length > 150:
            if repeats_list and repeats_list[-1][0] == length and repeats_list[-1][1] == start:
                repeats_list[-1].append(end)
            else:
                repeats_list.append([length, start, end])

            # Adding this range to ranges of repeat positions
            repeats_ranges.add(range(start, start + length))
            repeats_ranges.add(range(end, end + length))

    # Saving repeat ranges to file
    ranges_list = []
    for rng in repeats_ranges:
        ranges_list.append((list(rng)[0], list(rng)[-1], list(rng)[-1] - list(rng)[0] + 1))
    ranges_list.sort()
    with open("{}-ranges.txt".format(repeats_file_name[:-4]), "w") as repeats_ranges_file:
        for rng in ranges_list:
            repeats_ranges_file.write("{}\t{}\t{}\n".format(rng[0], rng[1], rng[2]))

    # Merging repeat ranges and writing them to file
    merged_ranges = []
    for begin, end, rng_len in ranges_list:
        if merged_ranges and merged_ranges[-1][1] >= begin - 1:
            merged_ranges[-1][1] = max(merged_ranges[-1][1], end)
            # Update repeat length
            merged_ranges[-1][2] = merged_ranges[-1][1] - merged_ranges[-1][0] + 1
        else:
            merged_ranges.append([begin, end, end - begin + 1])

    with open("{}-ranges-merged.txt".format(repeats_file_name[:-4]), "w") as repeats_ranges_file:
        for rng in merged_ranges:
            repeats_ranges_file.write("{}\t{}\t{}\n".format(rng[0], rng[1], rng[2]))

    return merged_ranges


# def find_repeats_snp(seq1_repeats_file_name, seq2_repeats_file_name, seq3_repeats_file_name, seq4_repeats_file_name,
#                      seq5_repeats_file_name, snps_file_name):
def find_repeats_snp(seq1_repeats_file_name, seq2_repeats_file_name, snps_file_name):
    seq1_merged_ranges = merge_repeat_ranges(seq1_repeats_file_name)
    seq2_merged_ranges = merge_repeat_ranges(seq2_repeats_file_name)
    # seq3_merged_ranges = merge_repeat_ranges(seq3_repeats_file_name)
    # seq4_merged_ranges = merge_repeat_ranges(seq4_repeats_file_name)
    # seq5_merged_ranges = merge_repeat_ranges(seq5_repeats_file_name)

    read_len = 300

    with open(snps_file_name, 'r') as snps_file:
        # Skip header line
        next(snps_file)
        for line in snps_file:
            fields = line.split()
            snp = fields[0]
            seq1_pos = int(fields[3])
            seq2_pos = int(fields[6])
            # seq3_pos = int(fields[9])
            # seq4_pos = int(fields[12])
            # seq5_pos = int(fields[15])

            snp_in_repeat = False
            msg = ""

            for start, end, length in seq1_merged_ranges:
                if seq1_pos in range(start, end + 1):
                    if seq1_pos - start < read_len or end - seq1_pos < read_len:
                        msg += "{} {} {} {} ".format(start, length, seq1_pos - start, end - seq1_pos)
                        snp_in_repeat = True
                        break
            if not snp_in_repeat:
                msg += "{} {} {} {} ".format('N', 'N', 'N', 'N')

            for start, end, length in seq2_merged_ranges:
                if seq2_pos in range(start, end + 1):
                    if seq2_pos - start < read_len or end - seq2_pos < read_len:
                        msg += "{} {} {} {} ".format(start, length, seq2_pos - start, end - seq2_pos)
                        snp_in_repeat = True
                        break
            if not snp_in_repeat:
                msg += "{} {} {} {} ".format('N', 'N', 'N', 'N')

            if snp_in_repeat:
                print("{} {}".format(snp, msg))


def k_mismatch_repeats(repeats_file_name, k=1, min_len=200):
    # [[length, start_1, (start_2, D), ...], ...]
    repeats_list = []
    with open(repeats_file_name) as repeats_file:
        # Skip header lines
        next(repeats_file)
        next(repeats_file)
        for line in repeats_file:
            fields = line.split()
            # Reputer repeats with distance k and on Forward strand
            # if int(fields[5]) == -k and fields[2] == 'F' and int(fields[0]) >= min_len:
            (start_1, start_2, length, direction, dist) = (int(fields[1]), int(fields[4]), int(fields[0]), fields[2],
                                                           int(fields[5]))

            if dist == -k and length >= min_len:
                if repeats_list:
                    (pre_len, pre_start_1, pre_start_2, pre_direction) = (repeats_list[-1][0], repeats_list[-1][1],
                                                                          repeats_list[-1][-1][0],
                                                                          repeats_list[-1][-1][1])
                    if pre_len == length:
                        # start_1 = prev start_1
                        if start_1 == pre_start_1:
                            # Add start_2
                            repeats_list[-1].append((start_2, direction))
                        # start_2 = prev start_2
                        elif start_2 == pre_start_2:
                            repeats_list[-1][1] = pre_start_2
                            repeats_list[-1][-1] = (pre_start_1, pre_direction)
                            repeats_list[-1].append((start_1, direction))
                        # start_2 = prev start_1
                        elif start_2 == repeats_list[-1][1]:
                            # Add start_1
                            repeats_list[-1].append((start_1, direction))
                        # It's a new group of repeats with the same length
                        else:
                            repeats_list.append([length, start_1, (start_2, direction)])
                    else:
                        repeats_list.append([length, start_1, (start_2, direction)])
                else:
                    repeats_list.append([length, start_1, (start_2, direction)])

    return repeats_list


def repeat_ranges(repeats_file_name, k=0, rep_len=1):
    """
    Returns a set which includes the ranges of all k-mismatch repeats
    :param repeats_file_name:
    :param k: 0 for exact repeats, k for k-mismatch repeats
    :param rep_len: minimum repeat length considered
    :return:
    """
    repeats_set = set()
    with open(repeats_file_name) as repeats_file:
        # Skip header lines
        next(repeats_file)
        next(repeats_file)
        for line in repeats_file:
            fields = line.split()
            # Reputer repeats with distance k and on Forward strand
            # if int(fields[5]) == -k and fields[2] == 'F' and int(fields[0]) >= min_len:
            (start_1, start_2, length, direction, dist) = (int(fields[1]), int(fields[4]), int(fields[0]), fields[2],
                                                           int(fields[5]))
            if dist == k and length >= rep_len:
                repeats_set.add(range(start_1, start_1 + length))
                repeats_set.add(range(start_2, start_2 + length))
    return repeats_set


def back_mutate_genome(ref_genome_file, repeats_file_name, output_file, read_len=150):
    """

    :param ref_genome_file:
    :param repeats_file_name:
    :param output_file:
    :return:
    """

    def calc_other_rep_locs(repeat_locs_lst, rep_len):
        """
        A helper function that finds the SNP position in other locations of this repeat element
        and considers reverse complement repeats as well.
        """
        other_repeat_locs = []
        for rep in repeat_locs_lst:
            rep_pos = rep[0]
            rep_dir = rep[1]
            if rep_dir == 'F':
                other_repeat_locs.append((rep_pos + i + 1, seq1[i]))
            else:
                other_repeat_locs.append((rep_pos + rep_len - i, reverse_complement(seq1[i])))

        return other_repeat_locs

    # Read reference genome
    genome_header, genome_seq = read_genome(ref_genome_file)
    # Extract the consolidated list of k-mismatch repeats
    repeats_list = k_mismatch_repeats(repeats_file_name, k=1, min_len=read_len)
    # Extract exact repeat ranges
    exact_repeats = repeat_ranges(repeats_file_name, k=0, rep_len=read_len)
    # Find the position of difference and modify the reference genome
    mutation_pos = set()
    mutations_dict = defaultdict(list)
    mutations_list = []
    new_genome_seq = list(genome_seq)
    for repeat in repeats_list:
        (length, start_1, start_2, direction) = (repeat[0], repeat[1], repeat[2][0], repeat[2][1])

        seq1 = genome_seq[start_1: start_1 + length]
        seq2 = genome_seq[start_2: start_2 + length]
        if direction == 'P':
            seq2 = reverse_complement(seq2)

        for i in range(length):
            if seq1[i] != seq2[i]:  # and (i < 150):

                # check whether the mutation point is not located itself in an exact repeat element
                in_exact_repeat = False
                for rng in exact_repeats:
                    if (start_1 + i) in rng:
                        in_exact_repeat = True
                        break
                if not in_exact_repeat:
                    if (start_1 + i) not in mutation_pos:
                        # The order: back-mutated nucleotide, true nucleotide
                        # mutations_list.append([start_1 + i + 1, seq2[i], seq1[i], i, repeat])
                        other_reps = calc_other_rep_locs(repeat[2:], length)
                        mutations_dict[start_1 + i + 1] = [seq2[i], seq1[i], i, other_reps, [length]]
                        # Mutating reference genome
                        new_genome_seq[start_1 + i] = seq2[i]
                        mutation_pos.add(start_1 + i)
                    # Overlapping repeats
                    else:
                        # print(start_1 + i + 1)
                        mutations_dict[start_1 + i + 1][4].append(length)
                        other_reps = calc_other_rep_locs(repeat[2:], length)
                        for rep in other_reps:
                            mutations_dict[start_1 + i + 1][3].append(rep)
                            # print(mutations_dict[start_1 + i + 1])

    # Writing mutated genome sequence to fasta file
    write_genome(genome_header, "".join(new_genome_seq), output_file)
    for key, value in mutations_dict.items():
        mutations_list.append([key] + value)

    # print(mutations_dict)
    # print(mutations_list)
    # Writing mutation locations to file
    mutations_list.sort()
    with open("{}-mutations.txt".format(output_file[:-4]), "w") as mutations_file:
        mutations_file.write("Number of mutations: {}\n".format(len(mutations_list)))
        for item in mutations_list:
            mutations_file.write(
                "{}\t{}\t{}\t{}\t{}\t{}\n".format(item[0], item[1], item[2], item[3], item[4], item[5]))

# toy_genome("./data/genomes/Mycobacterium_tuberculosis_H37Rv_uid57777/NC_000962.fna",
#            "./data/genomes/toy-genome-mutated.fna", mutate=True)

# toy_genome("./data/genomes/Mycobacterium_tuberculosis_H37Rv_uid57777/NC_000962.fna",
#            "./data/genomes/toy-genome-mutated-middle.fna", mutate=True)

# genome_slice("./data/genomes/Mycobacterium_tuberculosis_H37Rv_uid57777/NC_000962.fna",
#              [(1341821, 909, 'F'), (2982945, 909, 'P')])

# extract_genome("./data/genomes/Mycobacterium_tuberculosis_H37Rv_uid57777/NC_000962.fna", 1, 4411532,
#                "./data/genomes/mtb-whole-genome-mutated-100-140-half.fna", mutate=True,
#                repeats_file_name="./data/genomes/mtb-repeats-sorted.txt")

# extract_genome("./data/genomes/Mycobacterium_tuberculosis_H37Rv_uid57777/NC_000962.fna", 3948418, 230,
#                "/home/mohammad/test/first.fna")
#
# extract_genome("./data/genomes/Mycobacterium_tuberculosis_H37Rv_uid57777/NC_000962.fna", 3949618, 230,
#                "/home/mohammad/test/second.fna")

# merge_repeat_ranges("/home/mohammad/pneumoniae-genomes/h1repeats.txt")

# find_repeats_snp("/home/mohammad/pneumoniae/repeats/repeats-KPNIH1.txt",
#                  "/home/mohammad/pneumoniae/repeats/repeats-KPNIH10.txt",
#                  "/home/mohammad/pneumoniae/repeats/repeats-Kpn223.txt",
#                  "/home/mohammad/pneumoniae/repeats/repeats-Kpn555.txt",
#                  "/home/mohammad/pneumoniae/repeats/repeats-NJST258_1.txt",
#                  "/home/mohammad/pneumoniae/mauve-snps.txt")

# find_repeats_snp("/home/mohammad/pneumoniae/repeats-reputer/kp-kpninh1-repeats-reputer-100-ham2.txt",
#                  "/home/mohammad/pneumoniae/repeats-reputer/kp-kpninh10-repeats-reputer-100-ham2.txt",
#                  "/home/mohammad/pneumoniae/snps-mauve/snps-mauve-KPNH1-KPNH10.txt")

# rps = k_mismatch_repeats("/home/mohammad/pneumoniae/repeats-reputer/kp-kpninh1-repeats-reputer-100-ham2.txt")
# print(len(rps), '\n', rps)

# rps = k_mismatch_repeats("/home/mohammad/pneumoniae/repeats-reputer/mtb-repeats-reputer-100-ham2-filtered.txt")
# print(len(rps), '\n', rps)

# back_mutate_genome("/home/mohammad/pneumoniae/genomes/Klebsiella_pneumoniae_KPNIH1.fna",
#                   "/home/mohammad/pneumoniae/repeats-reputer/kp-kpninh1-repeats-reputer-100-ham2.txt",
#                   "/home/mohammad/pneumoniae/genomes/Klebsiella_pneumoniae_KPNIH1-back-mutated-full.fna")

# back_mutate_genome("./data/genomes/Klebsiella_pneumoniae_KPNIH1/Klebsiella_pneumoniae_KPNIH1.fna",
#                    "./data/genomes/repeats-reputer/kp-kpninh1-repeats-reputer-100-ham2.txt",
#                    "./data/genomes/Klebsiella_pneumoniae_KPNIH1-back-mutated-full.fna")

# back_mutate_genome("./data/genomes/Mycobacterium_tuberculosis_H37Rv_uid57777/NC_000962.fna",
#                    "./data/genomes/repeats-reputer/mtb-repeats-reputer-100-ham2.txt",
#                    "./data/genomes/MTB-H37Rv-back-mutated-full.fna", read_len=250)

# extract_genome("./data/genomes/Orientia_tsutsugamushi_Ikeda_uid58869/NC_010793.fna", 1, 2008987,
#                "./data/genomes/ot-whole-genome-mutated-70-140.fna", mutate=True,
#                repeats_file_name="./data/genomes/ot-repeats-sorted.txt")

# def merge_ranges(ranges_file):
#     ranges = []
#     with open(ranges_file) as infile:
#         for line in infile:
#             fields = line.split()
#             ranges.append((int(fields[0]), int(fields[1])))
#     a = sorted(ranges)
#     b = []
#     for begin, end in a:
#         if b and b[-1][1] >= begin - 1:
#             b[-1][1] = max(b[-1][1], end)
#         else:
#             b.append([begin, end])
#
#     with open("{}-merged.txt".format(ranges_file[:-4]), "w") as outfile:
#         for s, e in b:
#             outfile.write("{}\t{}\t{}\n".format(s, e, e - s + 1))
#
# merge_ranges("./data/genomes/mtb-repeats-sorted-ranges.txt")


# print(repeat_ranges("/home/mohammad/pneumoniae/repeats-reputer/kp-kpninh1-repeats-reputer-100-ham2.txt", k=0, rep_len=1))
