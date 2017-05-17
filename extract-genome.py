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

    # print("Repeats:", len(repeats_list))
    # for lst in repeats_list:
    #     print(lst)
    # print("_____________")

    # Saving repeat ranges to file
    ranges_list = []
    for rng in repeats_ranges:
        ranges_list.append((list(rng)[0], list(rng)[-1], list(rng)[-1] - list(rng)[0] + 1))
    ranges_list.sort()
    with open("{}-ranges.txt".format(repeats_file_name[:-4]), "w") as repeats_ranges_file:
        for rng in ranges_list:
            repeats_ranges_file.write("{}\t{}\t{}\n".format(rng[0], rng[1], rng[2]))

    # Mutating genome
    new_genome_seq = list(genome_sequence)
    mutations_list = []

    # 10 SNPs in repeated regions
    print("Number of repeat elements larger than 150 bp:", len(repeats_list))
    # num_snps = 10
    num_snps = len(repeats_list) // 5
    selected_regions = random.sample(repeats_list, num_snps)
    nucleotides = set(['A', 'C', 'G', 'T'])
    for region in selected_regions:
        mapping_location = random.choice(region[1:])
        base_position = random.choice(range(region[0]))
        final_position = mapping_location + base_position
        print(region, mapping_location, base_position, final_position)
        # Mutating the base
        possible_snps = nucleotides - set(genome_sequence[final_position])
        new_genome_seq[final_position] = random.choice(sorted(list(possible_snps)))
        # Saving this mutation characteristics
        mutations_list.append((final_position + 1, genome_sequence[final_position], new_genome_seq[final_position],
                               region))

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
            if len(item) > 3:
                mutations_file.write("{}\t{}\t{}\t{}\n".format(item[0], item[1], item[2], item[3]))
            else:
                mutations_file.write("{}\t{}\t{}\n".format(item[0], item[1], item[2]))

    # Returning mutated genome sequence
    return "".join(new_genome_seq)


# extract_genome("./data/genomes/Mycobacterium_tuberculosis_H37Rv_uid57777/NC_000962.fna", 3930000, 20000,
#                "./data/genomes/mtb-genome-extract.fna")

extract_genome("./data/genomes/Mycobacterium_tuberculosis_H37Rv_uid57777/NC_000962.fna", 1, 4411532,
               "./data/genomes/mtb-whole-genome-mutated.fna", mutate=True,
               repeats_file_name="./data/genomes/mtb-repeats-sorted.txt")


