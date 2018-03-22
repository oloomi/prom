import os
import sys


def bwa_complete_secondary(read_len, sam_file_name, output_file):
    """
    Inserts the sequence and quality scores for secondary alignments by looking up primary alignments
    """
    sam_col = {'qname': 0, 'flag': 1, 'rname': 2, 'pos': 3, 'mapq': 4, 'cigar': 5, 'seq': 9, 'qual': 10}
    # The header lines and the unique alignments will be directly written to the output file
    with open(sam_file_name) as sam_file:
        with open(output_file, 'w') as out_file:

            curr_line = next(sam_file)
            # Writing header lines to output SAM file
            while curr_line[0] == '@':
                out_file.write(curr_line)
                curr_line = next(sam_file)

            multimap = False
            prev_fields = curr_line.split('\t')

            # Read and preprocess each alignment
            for curr_line in sam_file:
                curr_fields = curr_line.split('\t')
                # Multimappings found for a previously seen multiread
                if curr_fields[sam_col['qname']] == prev_fields[sam_col['qname']]:
                    if curr_fields[sam_col['cigar']] != '{}M'.format(read_len):
                        continue
                    # BWA does not store the sequence and its quality score for secondary alignments
                    if curr_fields[sam_col['seq']] == '*' and curr_fields[sam_col['cigar']] != '*':
                        curr_fields[sam_col['seq']] = prev_fields[sam_col['seq']]
                        curr_fields[sam_col['qual']] = prev_fields[sam_col['qual']]

                out_file.write('\t'.join(prev_fields))
                prev_fields = curr_fields

            # Last line in the file
            out_file.write('\t'.join(prev_fields))


if __name__ == "__main__":
    if len(sys.argv) >= 4:
        bwa_complete_secondary(sys.argv[1], sys.argv[2], sys.argv[3])
    else:
        print("Read length and input and output file path required!")
