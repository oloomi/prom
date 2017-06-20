
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
            qual = float(fields[5])
            if qual > 0:
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


def compare_variants(benchmark_variants_file, vcf_files_list):
    """
    Compares two benchmark variants with the variants found by variant calling
    :return: Number of false positive, true positive, variants
    """
    output = ""
    output += "Method\tTruePositive\tFalsePositive\tFalseNegative\tF-Score\n"

    # Reading benchmark variants
    benchmark_variants = set(read_benchmark_variants(benchmark_variants_file))

    for vcf_file in vcf_files_list:
        method_name = vcf_file[0]
        vcf_file_name = vcf_file[1]
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

        output += "{}\t{}\t{}\t{}\t{:.2f}\n".format(method_name, tp, fp, fn, f1_score)
        if "remu" in vcf_file_name and True:
            output += "\nFalse negatives:\n{}\n".format(false_negatives)
            output += "False positives:\n{}\n\n".format(false_positives)

    return output
    # return (true_positives, false_positives, false_negatives)
