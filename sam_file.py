import sys
from collections import defaultdict
from genome_util import read_genome
from calc_likelihood import initial_counts

sam_col = {'qname': 0, 'pos': 3, 'cigar': 5, 'seq': 9, 'qual': 10}


def find_mdz_index(sam_fields):
    """
     Finding the index of MD:Z tag, as it can be different for different read mapping softwares
    """
    mdz_index = None
    for i, tag in enumerate(sam_fields[11:]):
        if 'MD:Z' in tag.upper():
            mdz_index = i
    return mdz_index


def add_multiread(sam_fields, multireads_dict, read_counts, mdz_index):
    """
     Checks the alignment of a multiread and add it to dictionary of multireads
    """
    # * means no alignment for a read
    if sam_fields[sam_col['cigar']] != "*":
        if sam_fields[sam_col['cigar']] == '{}M'.format(len(sam_fields[sam_col['seq']])):
            # Correct sequencing errors in the read
            read_seq, num_edits = correct_seq_errors(sam_fields, mdz_index)
            sam_fields[sam_col['seq']] = read_seq
            sam_fields[sam_col['pos']] = int(sam_fields[sam_col['pos']]) - 1
            # Store all alignments of a read, plus their edit distance
            multireads_dict[sam_fields[sam_col['qname']]].append(sam_fields + [num_edits])
        else:
            # Unsupported CIGARs
            read_counts['unsupported'].add(sam_fields[sam_col['qname']])
    else:
        read_counts['unmapped'].add(sam_fields[sam_col['qname']])
    return True


def process_unique_read(sam_fields, read_counts, base_counts, genome_seq):
    """
     Checks the alignment of a multiread and add it to dictionary of multireads
    """
    # * means no alignment for a read
    if sam_fields[sam_col['cigar']] != "*":
        if sam_fields[sam_col['cigar']] == '{}M'.format(len(sam_fields[sam_col['seq']])):
            sam_fields[sam_col['pos']] = int(sam_fields[sam_col['pos']]) - 1
            # Update pseudo-counts according to this unqiue mapping
            initial_counts(base_counts, sam_fields, genome_seq)
            read_counts['unique'] += 1
        else:
            # Unsupported CIGARs
            read_counts['unsupported'].add(sam_fields[sam_col['qname']])
    else:
        read_counts['unmapped'].add(sam_fields[sam_col['qname']])
    return True


def correct_seq_errors(sam_fields, mdz_index):
    (pos, md_z, read_seq, base_quals) = (sam_fields[sam_col['pos']], sam_fields[mdz_index][5:],
                                         sam_fields[sam_col['seq']], sam_fields[sam_col['qual']])
    num_edits = md_z.count('A') + md_z.count('C') + md_z.count('G') + md_z.count('T') + md_z.count('N')
    read_seq = list(read_seq)
    for i, base in enumerate(read_seq):
        # Phred-scale quality score
        if ord(base_quals[i]) - 33 < 20:
            read_seq[i] = 'N'
            # num_edits += 1

    return ''.join(read_seq), num_edits


def filter_alignments(mappings, threshold):
    """
    Filters mapping locations of a multi-read that have more edit operations than the best-match + threshold
    :param mappings: Initial list of multimappings
    :param threshold: an integer
    :return: A list of filtered mappings
    """
    edit_ops = [m[2] for m in mappings]

    # Finding best match
    min_ops = min(edit_ops)

    filtered_mappings = []
    for mapping in mappings:
        num_edits = mapping[2]
        # If this alignment is not too different from the best match
        if num_edits <= min_ops + threshold:
            filtered_mappings.append(mapping)

    return filtered_mappings


def read_all_mappings(sam_file_name, genome_seq, output_file, base_counts):
    # A dictionary like: {read_id: [list of mappings]}
    read_alignments_dict = defaultdict(list)

    read_counts = {'unsupported': set(), 'unmapped': set(), 'unique': 0}

    # Reading the SAM file and creating a dictionary of read_id : alignment
    # The header lines and the unique alignments will be directly written to the output file
    with open(sam_file_name) as sam_file:
        with open(output_file, 'w') as out_file:

            curr_line = next(sam_file)
            # Writing header lines to output SAM file
            while curr_line[0] == '@':
                out_file.write(curr_line)
                curr_line = next(sam_file)

            multimap = False
            unique_count = 0
            multireads_dict = defaultdict(list)
            prev_fields = curr_line.split('\t')
            mdz_index = find_mdz_index(prev_fields)
            if mdz_index is None:
                print('No MD:Z tag found in SAM fields!')
                sys.exit(1)

            # Read and preprocess each alignment
            for curr_line in sam_file:
                curr_fields = curr_line.split('\t')
                # Multimappings found for a previously seen multiread
                if curr_fields[sam_col['qname']] == prev_fields[sam_col['qname']]:
                    add_multiread(prev_fields, multireads_dict, read_counts, mdz_index)
                    multimap = True
                # Previous location is the last location for the seen multiread
                elif multimap == True:
                    add_multiread(prev_fields, multireads_dict, read_counts, mdz_index)
                    multimap = False
                else:
                    process_unique_read(prev_fields, read_counts, base_counts, genome_seq)
                    # add prev as unique
                prev_fields = curr_fields

            # Last line in the file
            if prev_fields[sam_col['qname']] in multireads_dict:
                add_multiread(prev_fields, multireads_dict, read_counts, mdz_index)
            else:
                process_unique_read(prev_fields, read_counts, base_counts, genome_seq)

    read_len = len(prev_fields[sam_col['seq']])

    return multireads_dict, read_len, read_counts


def read_sam_file(sam_file_name, genome_seq):
    # A dictionary like: {read_id: [list of mappings]}
    read_alignments_dict = defaultdict(list)

    total_reads = 0
    not_mapped_reads = set()  # Number of reads that are not mapped to any location
    invalid_cigars = set()  # Number of read mapping locations that contain indels, etc.

    # Reading the SAM file and creating a dictionary of read_id : alignment
    with open(sam_file_name) as sam_file:
        for line in sam_file:
            # Skip header lines
            if line[0] == "@":
                continue

            fields = line.rstrip().split("\t")
            read_id = fields[0]  # QNAME: Query template NAME
            cigar = fields[5]  # CIGAR string (ie. alignment)
            pos = int(fields[3])  # 1-based leftmost mapping POSition
            md_z = fields[-2][5:]  # Alignment eg. MD:Z:118C31 -> 118C31
            read_seq = fields[9]  # Read sequence
            base_quals = fields[10]  # Base quality scores
            # * means no alignment for a read
            if cigar != "*":
                if cigar == '{}M'.format(len(read_seq)):
                    # Correct sequencing errors in the read
                    read_seq, num_edits = correct_seq_errors((pos, md_z, read_seq, base_quals), genome_seq)
                    # Store all alignments of a read
                    read_alignments_dict[read_id].append((pos, cigar, num_edits, read_seq))
                else:
                    # print("Invalid CIGAR:", cigar)
                    invalid_cigars.add(read_id)
            else:
                not_mapped_reads.add(read_id)

    print("Number of reads with non {}M CIGAR:".format(len(read_seq)), len(invalid_cigars))
    print("Number of reads not mapped:", len(not_mapped_reads))
    print("Number of reads in use:", len(read_alignments_dict))

    return read_alignments_dict, len(read_seq)


def write_sam_file(multi_reads_correct_mapping, input_sam_file_name, output_sam_file_name):
    """
    Writes the mappings including best selected mapping for multi-reads to a SAM file
    """

    with open(input_sam_file_name) as in_sam_file:
        with open(output_sam_file_name, 'w') as out_sam_file:
            for line in in_sam_file:
                # Write header lines directly from input file to output file
                if line[0] == "@":
                    out_sam_file.write(line)
                    continue

                fields = line.rstrip().split("\t")
                read_id = fields[0]  # QNAME: Query template NAME
                cigar = fields[5]  # CIGAR string (ie. alignment)
                pos = int(fields[3])  # 1-based leftmost mapping POSition
                read_seq = fields[9]
                # * means no alignment for a read
                if cigar != "*":
                    if cigar == '{}M'.format(len(read_seq)):
                        if read_id in multi_reads_correct_mapping:
                            # If this mapping is the one at the selected mapping location
                            if multi_reads_correct_mapping[read_id] == pos:
                                out_sam_file.write(line)
                            else:
                                continue
                        else:
                            out_sam_file.write(line)
                    else:
                        continue  # ignore mappings with indels at the moment
                else:
                    out_sam_file.write(line)
    return True


def find_unique_reads(sam_file_name):
    # A dictionary like: {read_id: [list of mappings]}
    read_alignments_dict = defaultdict(list)

    # Reading the SAM file and creating a dictionary of read_id : alignment
    with open(sam_file_name) as sam_file:
        for line in sam_file:
            # Skip header lines
            if line[0] == "@":
                continue

            fields = line.rstrip().split("\t")
            read_id = fields[0]  # QNAME: Query template NAME
            cigar = fields[5]  # CIGAR string (ie. alignment)
            # * means no alignment for a read
            if cigar != "*":
                # Store all alignments of a read
                read_alignments_dict[read_id].append(cigar)

    with open(sam_file_name) as in_sam_file:
        with open("{}-unique.sam".format(sam_file_name[:-4]), 'w') as out_sam_file:
            for line in in_sam_file:
                # Write header lines directly from input file to output file
                if line[0] == "@":
                    out_sam_file.write(line)
                    continue

                fields = line.rstrip().split("\t")
                read_id = fields[0]  # QNAME: Query template NAME
                if read_id in read_alignments_dict and len(read_alignments_dict[read_id]) == 1:
                    # if read_id in read_alignments_dict and len(read_alignments_dict[read_id]) == 1 and \
                    #                 read_alignments_dict[read_id][0] == "150M":
                    out_sam_file.write(line)


def unique_reads_write_sam(sam_file_name, unique_reads, initially_resolved_multireads):
    with open(sam_file_name) as in_sam_file:
        with open("{}-unique.sam".format(sam_file_name[:-4]), 'w') as out_sam_file:
            for line in in_sam_file:
                # Write header lines directly from input file to output file
                if line[0] == "@":
                    out_sam_file.write(line)
                    continue

                fields = line.rstrip().split("\t")
                read_id = fields[0]  # QNAME: Query template NAME
                cigar = fields[5]  # CIGAR string (ie. alignment)
                pos = int(fields[3])  # 1-based leftmost mapping POSition
                read_seq = fields[9]
                # * means no alignment for a read
                if cigar == '{}M'.format(len(read_seq)):
                        if read_id in unique_reads:
                            out_sam_file.write(line)
                        elif read_id in initially_resolved_multireads and pos == initially_resolved_multireads[read_id]:
                            out_sam_file.write(line)
