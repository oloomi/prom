import ast


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
            chrom = fields[0]
            pos = int(fields[1])
            alt = fields[4]
            qual = float(fields[5])
            if qual > 0:
                variants.append((chrom, pos, alt))

    return variants


def read_benchmark_variants(benchmark_variants_file, read_len):
    """
    Returns a list of (position, variant) from the benchmark variations file
    """
    variants = []
    acceptable_fp_variants = []

    with open(benchmark_variants_file) as benchmark_variants:
        for line in benchmark_variants:
            # Skip comment lines
            if line[0] == '#':
                continue
            fields = line.rstrip().split("\t")
            chrom = fields[0]
            pos = int(fields[1])
            alt = fields[3]
            dist = int(fields[4])
            rep_len = int(fields[5])
            # Adding the benchmark variant
            variants.append((chrom, pos, alt))

            # If acceptable false positive locations is available
            if len(fields) > 6:
                other_locs = ast.literal_eval(fields[6])
                for pos_alt in other_locs:
                    acceptable_fp_variants.append(pos_alt)

    return variants, acceptable_fp_variants


def compare_variants(benchmark_variants_file, vcf_files_list, original_variants=None):
    """
    Compares two benchmark variants with the variants found by variant calling
    :return: Number of false positive, true positive, variants
    """
    output = ""
    output += "Method\tTruePositive\tFalsePositive\tFalseNegative\tF-Score\tRecall\tAcceptableFP\n"

    # Reading benchmark variants
    variants, acceptable_fp_variants = read_benchmark_variants(benchmark_variants_file, 150)
    benchmark_variants = set(variants)
    acceptable_fps = set(acceptable_fp_variants)

    # If mapping the reads to the original genome, not the back-mutated one
    if original_variants:
        original_variants = set(original_variants)

    for vcf_file in vcf_files_list:
        method_name = vcf_file[0]
        vcf_file_name = vcf_file[1]
        # If the extension is missing
        if vcf_file_name[-4:].lower() != ".vcf":
            vcf_file_name += ".vcf"
        called_variants = set()
        # Reading variants for this mapping
        called_variants = set(read_vcf_file(vcf_file_name))
        if not called_variants:
            print("No variants called in {} !".format(vcf_file))
            output += "{}\t{}\t{}\t{}\t{:.2f}\t{:.2f}\t{}\n".format(method_name, 0, 0, len(benchmark_variants), 0, 0, 0)
            continue
        # Measures
        true_positives = benchmark_variants & called_variants
        tp = len(true_positives)
        false_positives = called_variants - benchmark_variants
        fp = len(false_positives)
        false_negatives = benchmark_variants - called_variants
        fn = len(false_negatives)
        accept_fp = false_positives & acceptable_fps
        acceptable_fp = len(accept_fp)     # acceptable false positives

        if original_variants:
            actual_fp = false_positives - original_variants
            act_fp = len(actual_fp)
        else:
            act_fp = 'NA'

        # true_negatives: rest of the genome

        precision = tp / (tp + fp)
        recall = tp / (tp + fn)
        if precision + recall > 0:
            f1_score = 2 * precision * recall / (precision + recall)
        else:
            f1_score = 0

        output += "{}\t{}\t{}\t{}\t{:.2f}\t{:.2f}\t{}\n".format(method_name, tp, fp, fn, f1_score, recall, acceptable_fp)
        if "REMU" in vcf_file_name and False:
            output += "\nFalse negatives:\n{}\n".format(sorted(list(false_negatives)))
            output += "False positives:\n{}\n\n".format(sorted(list(false_positives)))

    return output
    # return (true_positives, false_positives, false_negatives)
