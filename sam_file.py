from collections import defaultdict


def read_sam_file(sam_file_name):
    # A dictionary like: {read_id: [list of mappings]}
    read_alignments_dict = defaultdict(list)

    total_reads = 0
    not_mapped_reads = 0  # Number of reads that are not mapped to any location
    invalid_cigars = 0  # Number of read mapping locations that contain indels, etc.

    # Reading the SAM file and creating a dictionary of read_id : alignment
    with open(sam_file_name) as sam_file:
        for line in sam_file:
            # Skip header lines
            if line[0] == "@":
                continue

            total_reads += 1
            fields = line.rstrip().split("\t")
            read_id = fields[0]  # QNAME: Query template NAME
            cigar = fields[5]  # CIGAR string (ie. alignment)
            pos = int(fields[3])  # 1-based leftmost mapping POSition
            md_z = fields[-2][5:]  # Alignment eg. MD:Z:118C31 -> 118C31
            read_seq = fields[9]  # Read sequence
            # * means no alignment for a read
            if cigar != "*":
                if cigar == '150M':
                    # Store all alignments of a read
                    read_alignments_dict[read_id].append((pos, cigar, md_z, read_seq))
                else:
                    # print("Invalid CIGAR:", cigar)
                    invalid_cigars += 1
            else:
                not_mapped_reads += 1

    print("Total number of reads:", total_reads)
    print("Number of reads with non 150M CIGAR:", invalid_cigars)
    print("Number of reads not mapped:", not_mapped_reads)
    print("Number of reads in use:", len(read_alignments_dict))

    return read_alignments_dict


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
                # * means no alignment for a read
                if cigar != "*":
                    if cigar == '150M':
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
                if read_id in read_alignments_dict and len(read_alignments_dict[read_id]) == 1 and \
                                read_alignments_dict[read_id][0] == "150M":
                    out_sam_file.write(line)
