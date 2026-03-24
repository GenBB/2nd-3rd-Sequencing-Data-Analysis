import sys
import pysam

def main():
    inbam = sys.argv[1] #  name-sorted BAM!!!!!!
    intag = sys.argv[2]

    insert_size_threshold = 1000
    read_type_counts_df, total_reads, total_mapped_bases = analyze_bam_name_sorted(inbam, insert_size_threshold)

    print(intag + "\t" + intag.split("_")[0] + "\t" + "Interchromosomal\t" + str(read_type_counts_df["Interchromosomal"]) + "\t" + str(total_reads) + "\t" + str(total_mapped_bases))
    print(intag + "\t" + intag.split("_")[0] + "\t" + "Inverted\t" + str(read_type_counts_df["Inverted"]) + "\t" + str(total_reads) + "\t" + str(total_mapped_bases))
    print(intag + "\t" + intag.split("_")[0] + "\t" + "Outward\t" + str(read_type_counts_df["Outward"]) + "\t" + str(total_reads) + "\t" + str(total_mapped_bases))
    print(intag + "\t" + intag.split("_")[0] + "\t" + "Large_Insert\t" + str(read_type_counts_df["Large_Insert"]) + "\t" + str(total_reads) + "\t" + str(total_mapped_bases))
    print(intag + "\t" + intag.split("_")[0] + "\t" + "Unclassified\t" + str(read_type_counts_df["Unclassified"]) + "\t" + str(total_reads) + "\t" + str(total_mapped_bases))
    print(intag + "\t" + intag.split("_")[0] + "\t" + "Normal\t" + str(read_type_counts_df["Normal"]) + "\t" + str(total_reads) + "\t" + str(total_mapped_bases))

def classify_read_pair_read12(read1, read2, insert_size_threshold):
    if read1.reference_name != read2.reference_name:
        return "Interchromosomal"

    strand_read1 = '-' if read1.is_reverse else '+'
    strand_read2 = '-' if read2.is_reverse else '+'
    read1_pos = read1.reference_start
    read2_pos = read2.reference_start
    insert_size = abs(read1.template_length)

    if strand_read1 == strand_read2:
        return "Inverted"
    elif (
        (strand_read1 == '-' and strand_read2 == '+' and read1_pos < read2_pos) or
        (strand_read1 == '+' and strand_read2 == '-' and read2_pos < read1_pos)
    ):
        return "Outward"
    elif insert_size > insert_size_threshold:
        return "Large_Insert"
    else:
        return "Normal"

def analyze_bam_name_sorted(inbam, insert_size_threshold):
    total_reads = 0
    total_mapped_bases = 0

    read_type_counts_df = {
        "Interchromosomal": 0,
        "Inverted": 0,
        "Outward": 0,
        "Large_Insert": 0,
        "Normal": 0,
        "Unclassified": 0
    }

    with pysam.AlignmentFile(inbam, "rb") as bamfile:
        buffered_read = None
        current_read_name = None

        for read in bamfile:
            if read.is_secondary or read.is_supplementary:
                continue

            if current_read_name == read.query_name:
                
                read1 = buffered_read if buffered_read.is_read1 else read
                read2 = read if buffered_read.is_read1 else buffered_read

                total_reads += 2
                total_mapped_bases += read1.query_alignment_length
                total_mapped_bases += read2.query_alignment_length

                result = classify_read_pair_read12(read1, read2, insert_size_threshold)
                if result:
                    read_type_counts_df[result] += 1
                else:
                    read_type_counts_df["Unclassified"] += 1

                
                buffered_read = None
                current_read_name = None
            else:
                
                buffered_read = read
                current_read_name = read.query_name

    return read_type_counts_df, total_reads, total_mapped_bases

if __name__ == "__main__":
    main()
