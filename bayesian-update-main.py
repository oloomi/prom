def extract_genome(ref_genome_file, start_pos, length, output_file):
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

            # Writing the new genome sequence to file, 70 characters per line
            new_genome_seq = genome_seq[start_pos - 1: start_pos + length - 1]
            line_width = 70
            for i in range(length // line_width):
                new_genome.write(new_genome_seq[i * line_width: (i + 1) * line_width])
                new_genome.write("\n")
            # Writing the last remainder part of genome
            if length % line_width != 0:
                new_genome.write(new_genome_seq[-(length % line_width):])
    return True


def bayesian_update():
    return True


# /mnt/e/Codes/data/genomes/Mycobacterium_tuberculosis_H37Rv_uid57777
# extract_genome("../data/genomes/Mycobacterium_tuberculosis_H37Rv_uid57777/NC_000962.fna", 5, 141, "test-mtb-genome.txt")
extract_genome("../data/genomes/Mycobacterium_tuberculosis_H37Rv_uid57777/NC_000962.fna", 3930000, 20000, "mtb-genome-extract.txt")