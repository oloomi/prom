import random


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
                    #header_line = line.rstrip().split("|")
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
            if length > 150:
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


    # Mutating genome
    new_genome_seq = list(genome_sequence)
    mutations_list = []

    # 25% of regions
    num_snps = len(merged_ranges) // 4
    selected_regions = random.sample(merged_ranges, num_snps)
    nucleotides = set(['A', 'C', 'G', 'T'])
    for region in selected_regions:
        mapping_location = region[0]
        base_position = random.choice(range(70, 140))
        final_position = mapping_location + base_position
        print(region, mapping_location, base_position, final_position)
        # Mutating the base
        possible_snps = nucleotides - set(genome_sequence[final_position])
        new_genome_seq[final_position] = random.choice(sorted(list(possible_snps)))
        # Saving this mutation characteristics
        mutations_list.append((final_position + 1, genome_sequence[final_position], new_genome_seq[final_position],
                               region))


    #     # 10 SNPs in repeated regions
    # print("Number of repeat elements larger than 150 bp:", len(repeats_list))
    # # num_snps = 10
    # num_snps = len(repeats_list) // 5
    # selected_regions = random.sample(repeats_list, num_snps)
    # nucleotides = set(['A', 'C', 'G', 'T'])
    # for region in selected_regions:
    #     mapping_location = random.choice(region[1:])
    #     # base_position = random.choice(range(region[0]))
    #     base_position = random.choice(range(70, 140))
    #     final_position = mapping_location + base_position
    #     print(region, mapping_location, base_position, final_position)
    #     # Mutating the base
    #     possible_snps = nucleotides - set(genome_sequence[final_position])
    #     new_genome_seq[final_position] = random.choice(sorted(list(possible_snps)))
    #     # Saving this mutation characteristics
    #     mutations_list.append((final_position + 1, genome_sequence[final_position], new_genome_seq[final_position],
    #                            region))

    # 10 SNPs in normal regions
    # num_normal_mutations = 10
    # normal_mutations = []
    # i = 0
    # while i < num_normal_mutations:
    #     final_position = random.choice(range(genome_length))
    #     in_repeats = False
    #     for repeat_range in repeats_ranges:
    #         if final_position in repeat_range:
    #             in_repeats = True
    #             break
    #     if not in_repeats and final_position not in normal_mutations:
    #         # Mutating the base
    #         possible_snps = nucleotides - set(genome_sequence[final_position])
    #         new_genome_seq[final_position] = random.choice(sorted(list(possible_snps)))
    #         # Saving the mutation characteristics
    #         mutations_list.append(
    #             (final_position + 1, genome_sequence[final_position], new_genome_seq[final_position]))
    #         normal_mutations.append(final_position)
    #         i += 1

    # Writing mutation locations to file
    mutations_list.sort()
    with open(mutation_locations_file, "w") as mutations_file:
        mutations_file.write("Number of mutations: {}\n".format(num_snps))
        for item in mutations_list:
            mutations_file.write("{}\t{}\t{}\t{}\n".format(item[0], item[1], item[2], item[3]))
            # if len(item) > 3:
            #     mutations_file.write("{}\t{}\t{}\t{}\t{}\n".format(item[0], item[1], item[2], item[3][0], item[3][1:]))
            # else:
            #     mutations_file.write("{}\t{}\t{}\n".format(item[0], item[1], item[2]))

    # Returning mutated genome sequence
    return "".join(new_genome_seq)


# extract_genome("./data/genomes/Mycobacterium_tuberculosis_H37Rv_uid57777/NC_000962.fna", 3930000, 20000,
#                "./data/genomes/mtb-genome-extract.fna")

# extract_genome("./data/genomes/Mycobacterium_tuberculosis_H37Rv_uid57777/NC_000962.fna", 1, 4411532,
#                "./data/genomes/mtb-whole-genome-mutated-70-140.fna", mutate=True,
#                repeats_file_name="./data/genomes/mtb-repeats-sorted.txt")

extract_genome("./data/genomes/Orientia_tsutsugamushi_Ikeda_uid58869/NC_010793.fna", 1, 2008987,
               "./data/genomes/ot-whole-genome-mutated-70-140.fna", mutate=True,
               repeats_file_name="./data/genomes/ot-repeats-sorted.txt")

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