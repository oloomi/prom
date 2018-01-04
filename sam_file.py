import sys
import os
from collections import defaultdict
from settings import *
from calc_likelihood import initial_counts, update_counts


def read_genome(genome_file):
    genome_header = ""
    genome_seq = ""
    with open(genome_file) as ref_genome:
        for line in ref_genome:
            # Skip header lines
            if line[0] == ">":
                genome_header += line
            else:
                genome_seq += line.rstrip()
    return genome_header, genome_seq


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
            num_edits = find_edit_ops(sam_fields, mdz_index)
            # sam_fields[sam_col['seq']] = read_seq
            sam_fields[sam_col['pos']] = int(sam_fields[sam_col['pos']]) - 1
            # Store all alignments of a read, plus their edit distance
            # NOTE: edit distance is appended to alignment features
            multireads_dict[sam_fields[sam_col['qname']]].append(sam_fields + [num_edits])
            read_counts['multi'] += 1
        else:
            # Unsupported CIGARs
            read_counts['unsupported'].add(sam_fields[sam_col['qname']])
    else:
        read_counts['unmapped'].add(sam_fields[sam_col['qname']])
    return True


def process_unique_read(sam_fields, read_counts, base_counts, genome_seq, outfile):
    """
     Checks the alignment of a multiread and add it to dictionary of multireads
    """
    # * means no alignment for a read
    if sam_fields[sam_col['cigar']] != "*":
        if sam_fields[sam_col['cigar']] == '{}M'.format(len(sam_fields[sam_col['seq']])):
            sam_fields[sam_col['pos']] = int(sam_fields[sam_col['pos']]) - 1
            # Update pseudo-counts according to this unqiue mapping
            initial_counts(base_counts, sam_fields, genome_seq)
            # update_counts(base_counts, sam_fields, 24, genome_seq)
            read_counts['unique'] += 1
        else:
            sam_fields[sam_col['pos']] = int(sam_fields[sam_col['pos']]) - 1
            # Unsupported CIGARs
            read_counts['unsupported'].add(sam_fields[sam_col['qname']])
        # Write alignment to output file
        outfile.write("\t".join(sam_fields[0:sam_col['pos']] + [str(sam_fields[sam_col['pos']] + 1)] +
                                sam_fields[sam_col['pos']+1:]))
    else:
        read_counts['unmapped'].add(sam_fields[sam_col['qname']])
    return True


def find_edit_ops(sam_fields, mdz_index):
    md_z = sam_fields[mdz_index][5:]
    num_edits = md_z.count('A') + md_z.count('C') + md_z.count('G') + md_z.count('T') + md_z.count('N')
    # read_seq = list(read_seq)
    # for i, base in enumerate(read_seq):
    #     # Phred-scale quality score
    #     if ord(base_quals[i]) - 33 < 20:
    #         read_seq[i] = 'N'
    #         # num_edits += 1

    return num_edits


def filter_alignments(mappings, threshold):
    """
    Filters mapping locations of a multi-read that have more edit operations than the best-match + threshold
    :param mappings: Initial list of multimappings
    :param threshold: an integer
    :return: A list of filtered mappings
    """
    edit_ops = [m[-1] for m in mappings]

    # Finding best match
    min_ops = min(edit_ops)

    filtered_mappings = []
    for mapping in mappings:
        num_edits = mapping[-1]
        # Remove the edit distance and reserve two values for each mapping location:
        # sum of probabilities, and current run/iteration probability
        mapping[-1] = 0
        mapping.append(0)
        # If this alignment is not too different from the best match
        if num_edits <= min_ops + threshold:
            filtered_mappings.append(mapping)

    return filtered_mappings


def read_all_mappings(sam_file_name, genome_seq, output_file, base_counts):
    read_counts = {'unsupported': set(), 'unmapped': set(), 'unique': 0, 'multi': 0}

    # Reading the SAM file and creating a dictionary of read_id : alignment
    # The header lines and the unique alignments will be directly written to the output file
    with open(sam_file_name) as sam_file:
        with open(output_file, 'w') as out_file:
            file_size = os.path.getsize(sam_file_name)

            curr_line = next(sam_file)
            progress = len(curr_line)
            # Writing header lines to output SAM file
            while curr_line[0] == '@':
                out_file.write(curr_line)
                curr_line = next(sam_file)
                progress += len(curr_line)

            multimap = False
            unique_count = 0
            multireads_dict = defaultdict(list)
            prev_fields = curr_line.split('\t')
            mdz_index = find_mdz_index(prev_fields)
            if mdz_index is None:
                print('No MD:Z tag found in SAM fields!')
                sys.exit(1)

            progress_percentage = round(progress / file_size * 100)

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
                # It's a uniquely mapped read!
                else:
                    process_unique_read(prev_fields, read_counts, base_counts, genome_seq, out_file)

                prev_fields = curr_fields

                # Progress percentage
                progress += len(curr_line)
                new_prog_perc = round(progress / file_size * 100)
                if new_prog_perc != progress_percentage:
                    sys.stdout.write('\r')
                    sys.stdout.write(" {}%".format(new_prog_perc))
                    progress_percentage = new_prog_perc

            # Last line in the file
            if prev_fields[sam_col['qname']] in multireads_dict:
                add_multiread(prev_fields, multireads_dict, read_counts, mdz_index)
            else:
                process_unique_read(prev_fields, read_counts, base_counts, genome_seq, out_file)

            # For new line character after 100% progress
            print('')

    read_len = len(prev_fields[sam_col['seq']])

    return multireads_dict, read_len, read_counts
