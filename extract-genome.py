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
                    header_line = line.rstrip().split("|")
                    # Update the genome length field in the header line
                    header_line[1] = str(length)
                    new_genome.write("{} {}:{}\n".format("|".join(header_line), start_pos, start_pos + length - 1))
                else:
                    genome_seq += line.rstrip()

            # Extracted genome sequence
            new_genome_seq = genome_seq[start_pos - 1: start_pos + length - 1]

            if mutate:
                new_genome_seq = mutate_genome(repeats_file_name, new_genome_seq, start_pos, length,
                                                "mutations-{}.txt".format(output_file[:-4]))

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
            # Deducting start_pos to make absolute positions for the new genome extract
            (start, end, length) = (int(fields[0]) - start_pos, int(fields[1]) - start_pos, int(fields[2]))
            if repeats_list and repeats_list[-1][0] == length and repeats_list[-1][1] == start:
                repeats_list[-1].append(end)
            else:
                repeats_list.append([length, start, end])

            # Adding this range to ranges of repeat positions
            repeats_ranges.add(range(start, start + length))

    # print("Repeats:", len(repeats_list))
    # for lst in repeats_list:
    #     print(lst)
    # print("_____________")

    new_genome_seq = list(genome_sequence)

    with open(mutation_locations_file, "w") as mutations_file:
        # 10 SNPs in repeated regions
        selected_regions = random.sample(repeats_list, 10)
        nucleotides = set(['A', 'C', 'G', 'T'])
        for region in selected_regions:
            mapping_location = random.choice(region[1:])
            base_position = random.choice(range(region[0]))
            final_position = mapping_location + base_position
            # Mutating the base
            possible_snps = nucleotides - set(genome_sequence[final_position])
            new_genome_seq[final_position] = random.choice(list(possible_snps))
            # Saving this mutation characteristics to file
            mutations_file.write("{}\t{}\t{}\t{}\n".format(final_position + 1, genome_sequence[final_position],
                                                     new_genome_seq[final_position], region))

        # 10 SNPs in normal regions
        num_normal_mutations = 10
        normal_mutations = []
        i = 0
        while i < num_normal_mutations:
            final_position = random.choice(range(genome_length))
            in_repeats = False
            for repeat_range in repeats_ranges:
                if final_position in repeat_range:
                    in_repeats == True
                    break
            if not in_repeats and final_position not in normal_mutations:
                # Mutating the base
                possible_snps = nucleotides - set(genome_sequence[final_position])
                new_genome_seq[final_position] = random.choice(list(possible_snps))
                mutations_file.write("{}\t{}\t{}\n".format(final_position + 1, genome_sequence[final_position],
                                                         new_genome_seq[final_position]))
                normal_mutations.append(final_position)
                i += 1

    return "".join(new_genome_seq)

# /mnt/e/Codes/data/genomes/Mycobacterium_tuberculosis_H37Rv_uid57777
# extract_genome("../data/genomes/Mycobacterium_tuberculosis_H37Rv_uid57777/NC_000962.fna", 5, 141, "test-mtb-genome.txt")
extract_genome("../data/genomes/Mycobacterium_tuberculosis_H37Rv_uid57777/NC_000962.fna", 3930000, 20000,
               "mtb-genome-extract-mutated.fna", mutate=True, repeats_file_name="./read-mapping/mtb-selected-region-repeats.txt")
# mutate_repeats("./read-mapping/mtb-selected-region-repeats.txt", "", 3930000, 20000, "")