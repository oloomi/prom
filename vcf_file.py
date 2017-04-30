
def read_vcf_file(vcf_file_name):
    """
    Returns a list of (position, variant) from a VCF file
    """
    variants = []
    with open(vcf_file_name) as vcf_file:
        for line in vcf_file:
            # Skip header lines
            if line[0] == "#":
                continue
            # 0.CHROM   1.POS   2.ID    3.REF   4.ALT	5.QUAL	6.FILTER	7.INFO  8.FORMAT
            fields = line.rstrip().split("\t")
            pos = fields[1]
            alt = fields[4]
            variants.append((pos, alt))

    return variants


def read_benchmark_variants(benchmark_variants_file):
    """
    Returns a list of (position, variant) from the benchmark variations file
    """
    variants = []
    with open(benchmark_variants_file) as benchmark_variants:
        for line in benchmark_variants:
            fields = line.rstrip().split("\t")
            pos = fields[0]
            alt = fields[2]
            variants.append((pos, alt))

    return variants


def compare_variants(benchmark_variants_file, vcf_file_name, output_file_name):
    """
    Compares two benchmark variants with the variants found by variant calling
    :return: Number of false positive, true positive, variants
    """
    benchmark_variants = set(read_benchmark_variants(benchmark_variants_file))
    called_variants = set(read_vcf_file(vcf_file_name))

    true_positives = benchmark_variants & called_variants
    false_positives = called_variants - benchmark_variants
    false_negatives = benchmark_variants - called_variants
    # true_negatives: rest of the genome

    with open(output_file_name, 'w') as output_file:
        output_file.write("True Positives\t{}\n".format(len(true_positives)))
        output_file.write("{}\n".format(true_positives))
        output_file.write("False Positives\t{}\n".format(len(false_positives)))
        output_file.write("{}\n".format(false_positives))
        output_file.write("False Negatives\t{}\n".format(len(false_negatives)))
        output_file.write("{}\n".format(false_negatives))

    # return (true_positives, false_positives, false_negatives)
