import random
from collections import defaultdict


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
            print(reverse_complement(genome_seq[loc[0] : loc[0]+loc[1]]))
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
            # pos = 2129
            pos = 2500
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

            # for start, end, length in seq3_merged_ranges:
            #     if seq3_pos in range(start, end + 1):
            #         if seq3_pos - start < read_len or end - seq3_pos < read_len:
            #             msg += "{} {} {} {} ".format(start, length, seq3_pos - start, end - seq3_pos)
            #             snp_in_repeat = True
            # if not snp_in_repeat:
            #     msg += "{} {} {} {} ".format('N', 'N', 'N', 'N')
            #
            # for start, end, length in seq4_merged_ranges:
            #     if seq4_pos in range(start, end + 1):
            #         if seq4_pos - start < read_len or end - seq4_pos < read_len:
            #             msg += "{} {} {} {} ".format(start, length, seq4_pos - start, end - seq4_pos)
            #             snp_in_repeat = True
            #             break
            # if not snp_in_repeat:
            #     msg += "{} {} {} {} ".format('N', 'N', 'N', 'N')
            #
            # for start, end, length in seq5_merged_ranges:
            #     if seq5_pos in range(start, end + 1):
            #         if seq5_pos - start < read_len or end - seq5_pos < read_len:
            #             msg += "{} {} {} {} ".format(start, length, seq5_pos - start, end - seq5_pos)
            #             snp_in_repeat = True
            #             break
            # if not snp_in_repeat:
            #     msg += "{} {} {} {} ".format('N', 'N', 'N', 'N')

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


def back_mutate_genome(ref_genome_file, repeats_file_name, output_file):
    # Read reference genome
    genome_header, genome_seq = read_genome(ref_genome_file)
    # Extract the consolidated list of k-mismatch repeats
    repeats_list = k_mismatch_repeats(repeats_file_name, k=1, min_len=200)
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
            if seq1[i] != seq2[i]:   #and (i < 150):
                if (start_1 + i) not in mutation_pos:
                    # The order: back-mutated nucleotide, true nucleotide
                    # mutations_list.append([start_1 + i + 1, seq2[i], seq1[i], i, repeat])
                    mutations_dict[start_1 + i + 1] = [seq2[i], seq1[i], i, repeat[2:], [length]]
                    # Mutating reference genome
                    new_genome_seq[start_1 + i] = seq2[i]
                    mutation_pos.add(start_1 + i)
                else:
                    print(start_1 + i + 1)
                    mutations_dict[start_1 + i + 1][4].append(length)
                    for r in repeat[2:]:
                        mutations_dict[start_1 + i + 1][3].append(r)
                    print(mutations_dict[start_1 + i + 1])

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
            mutations_file.write("{}\t{}\t{}\t{}\t{}\t{}\n".format(item[0], item[1], item[2], item[3], item[4], item[5]))



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

back_mutate_genome("/home/mohammad/pneumoniae/genomes/Klebsiella_pneumoniae_KPNIH1.fna",
                   "/home/mohammad/pneumoniae/repeats-reputer/kp-kpninh1-repeats-reputer-100-ham2.txt",
                   "/home/mohammad/pneumoniae/genomes/Klebsiella_pneumoniae_KPNIH1-back-mutated-full.fna")

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

