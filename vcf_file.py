
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
            # Skip header line
            if line[0].isalpha():
                continue
            fields = line.rstrip().split("\t")
            pos = fields[0]
            alt = fields[2]
            variants.append((pos, alt))

    return variants


def compare_variants(benchmark_variants_file, vcf_files_list, output_file_name):
    """
    Compares two benchmark variants with the variants found by variant calling
    :return: Number of false positive, true positive, variants
    """
    output_file = open(output_file_name, 'w')
    output_file.write("Method\tTruePositive\tFalsePositive\tFalseNegative\tF-Score\n")

    # Reading benchmark variants
    benchmark_variants = set(read_benchmark_variants(benchmark_variants_file))

    for vcf_file_name in vcf_files_list:
        # If the extension is missing
        if vcf_file_name[-4:].lower() != ".vcf":
            vcf_file_name += ".vcf"
        called_variants = set()
        # Reading variants for this mapping
        called_variants = set(read_vcf_file(vcf_file_name))
        # Measures
        true_positives = benchmark_variants & called_variants
        tp = len(true_positives)
        false_positives = called_variants - benchmark_variants
        fp = len(false_positives)
        false_negatives = benchmark_variants - called_variants
        fn = len(false_negatives)
        # true_negatives: rest of the genome

        precision = tp / (tp + fp)
        recall = tp / (tp + fn)
        if precision + recall > 0:
            f1_score = 2 * precision * recall / (precision + recall)
        else:
            f1_score = 0

        output_file.write("{}\t{}\t{}\t{}\t{:.2f}\n".format(vcf_file_name.split("/")[-1][:-4], tp, fp, fn, f1_score))
        if "corrected-mtb" in vcf_file_name:
            output_file.write("\nFalse negatives:\n{}\n".format(false_negatives))
            output_file.write("False positives:\n{}\n\n".format(false_positives))

    # return (true_positives, false_positives, false_negatives)
