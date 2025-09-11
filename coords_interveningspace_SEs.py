import os


def parse_fasta_header(header):
    """Parses the FASTA header to extract chromosome and coordinates."""
    parts = header.split(':')
    chromosome = parts[0]
    start, end = map(int, parts[1].split('-'))
    return chromosome, start, end


def find_intergenic_regions(headers, filenames):
    """Finds the intergenic regions between consecutive FASTA sequences and includes filenames."""
    headers = sorted(zip(headers, filenames), key=lambda x: (x[0][0], x[0][1]))
    intergenic_regions = []

    for i in range(len(headers) - 1):
        (chr1, start1, end1), file1 = headers[i]
        (chr2, start2, end2), file2 = headers[i + 1]

        if chr1 != chr2:
            print(f"Skipping: {chr1} != {chr2}")
            continue  # Skip if not on the same chromosome

        if end1 >= start2:
            print(f"Skipping overlapping or unordered regions: {headers[i][0]} and {headers[i + 1][0]}")
            continue  # Skip if sequences overlap or are not in sequential order

        intergenic_start = end1 + 1
        intergenic_end = start2 - 1
        intergenic_regions.append((chr1, intergenic_start, intergenic_end, f"{file1}_{file2}"))

    return intergenic_regions


def write_bed_file(intergenic_regions, output_file, subdir_name):
    """Writes the intergenic regions to a BED format file with subdirectory name and filenames."""
    with open(output_file, 'w') as f:
        for chromosome, start, end, filename_pair in intergenic_regions:
            f.write(f"{chromosome}\t{start}\t{end}\t{subdir_name}\t{filename_pair}\n")
    print(f"Written to {output_file}")


def read_fasta_headers_from_directory(directory):
    """Reads headers and filenames from all FASTA files in a directory."""
    headers = []
    filenames = []
    for filename in os.listdir(directory):
        if filename.endswith(".txt"):
            file_path = os.path.join(directory, filename)
            with open(file_path, 'r') as f:
                header = f.readline().strip().lstrip('>')
                headers.append(parse_fasta_header(header))
                filenames.append(os.path.splitext(filename)[0])
    return headers, filenames


def process_subdirectories(parent_directory):
    """Processes each subdirectory in the parent directory."""
    for subdir in os.listdir(parent_directory):
        subdir_path = os.path.join(parent_directory, subdir)
        if os.path.isdir(subdir_path):
            print(f"Processing subdirectory: {subdir_path}")
            headers, filenames = read_fasta_headers_from_directory(subdir_path)
            if len(headers) > 1:
                print(f"Headers found: {headers}")
                intergenic_regions = find_intergenic_regions(headers, filenames)
                print(f"Intergenic regions: {intergenic_regions}")
                output_file = os.path.join(subdir_path, "intervening_seq.txt")
                write_bed_file(intergenic_regions, output_file, subdir)
            else:
                print(f"Not enough files in {subdir_path} to find intergenic regions.")


def count_single_file_directories(parent_directory):
    """Counts the number of directories with only one file in the parent directory."""
    single_file_count = 0
    for subdir in os.listdir(parent_directory):
        subdir_path = os.path.join(parent_directory, subdir)
        if os.path.isdir(subdir_path):
            files = [f for f in os.listdir(subdir_path) if f.endswith(".txt")]
            if len(files) == 1:
                single_file_count += 1
    print(f"Number of directories with only one file: {single_file_count}")
    return single_file_count

def concatenate_intervening_seq_files(parent_directory, final_output_file):
    """Concatenates all intervening_seq.txt files from subdirectories into a single BED file."""
    with open(final_output_file, 'w') as outfile:
        for subdir in os.listdir(parent_directory):
            subdir_path = os.path.join(parent_directory, subdir)
            intervening_file = os.path.join(subdir_path, "intervening_seq.txt")
            if os.path.isfile(intervening_file):
                with open(intervening_file, 'r') as infile:
                    for line in infile:
                        outfile.write(line)
    print(f"All intervening sequences concatenated into {final_output_file}")


#parent_directory = "/Users/jillianness/Desktop/sorting_cannavo_data/TE_splitting/ALL_shadowsets_dm3"  # Replace with the actual path to your parent directory
parent_directory = "/Users/jillianness/Desktop/mouse_analysis_031925/Mouse_enhancers/Fasta_files/Shadow_fasta_files/FlilesOutput"

process_subdirectories(parent_directory)
final_bed_path = os.path.join(parent_directory, "all_intervening_sequences.bed")
concatenate_intervening_seq_files(parent_directory, final_bed_path)