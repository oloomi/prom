import random


def extract_genome(ref_genome_file, start_pos, length, output_file, mutate=False, genome_repeats=None,
                   repeats_file_name=None):
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
                new_genome_seq = mutate_repeats(repeats_file_name, new_genome_seq, start_pos, "mutations-{}.txt".format(output_file[:-3]))

            # Writing the new genome sequence to file, 70 characters per line
            line_width = 70
            for i in range(length // line_width):
                new_genome.write(new_genome_seq[i * line_width: (i + 1) * line_width])
                new_genome.write("\n")
            # Writing the last remainder part of genome
            if length % line_width != 0:
                new_genome.write(new_genome_seq[-(length % line_width):])
    return True


def mutate_repeats(repeats_file_name, genome_sequence, start_pos, mutation_locations_file):
    random.seed(12)

    # [[length, start_1, start_2, ...], ...]
    repeats_list = []
    with open(repeats_file_name) as repeats_file:
        for line in repeats_file:
            fields = line.split()
            # Deducting start_pos to make absolute positions for the new genome extract
            (start, end, length) = (int(fields[0]) - start_pos, int(fields[1]) - start_pos, int(fields[2]))
            if repeats_list and repeats_list[-1][0] == length and repeats_list[-1][1] == start:
                repeats_list[-1].append(end)
            else:
                repeats_list.append([length, start, end])

    print("Repeats:", len(repeats_list))
    for lst in repeats_list:
        print(lst)
    print("_____________")

    selected_regions = random.sample(repeats_list, 20)
    for region in selected_regions:
        mapping_location = random.choice(region[1:])
        base_position = random.choice(range(region[0]))
        final_position = mapping_location + base_position
        print(region)
        print(mapping_location, base_position, final_position)


    return True

# /mnt/e/Codes/data/genomes/Mycobacterium_tuberculosis_H37Rv_uid57777
# extract_genome("../data/genomes/Mycobacterium_tuberculosis_H37Rv_uid57777/NC_000962.fna", 5, 141, "test-mtb-genome.txt")
# extract_genome("../data/genomes/Mycobacterium_tuberculosis_H37Rv_uid57777/NC_000962.fna", 3930000, 20000, "mtb-genome-extract.txt")
mutate_repeats("./read-mapping/mtb-selected-region-repeats.txt", "", 3930000, "")