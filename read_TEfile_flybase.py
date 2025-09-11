import re


def extract_location_and_name(header):
    # Regular expression to extract the location
    loc_pattern = r"loc=([\w]+):(\d+)\.\.(\d+)"
    loc_match = re.search(loc_pattern, header)

    # Regular expression to extract the name before {}
    name_pattern = r"name=([^{}]+)\{"
    name_match = re.search(name_pattern, header)

    if loc_match and name_match:
        chromosome = loc_match.group(1)
        start = int(loc_match.group(2))
        end = int(loc_match.group(3))
        name = name_match.group(1)
        return {
            'chromosome': chromosome,
            'start': start,
            'end': end,
            'name': name.strip()  # Stripping any extra whitespace
        }
    else:
        return None


def process_headers_from_file(file_path):
    with open(file_path, 'r') as file:
        # Read the whole file content
        content = file.read()

        # Split the content by '>' (which starts each header) and remove any empty strings
        headers = [header for header in content.split('>') if header.strip()]

        results = []

        # Process each header part
        for header in headers:
            header_part = header.split(";")[0]  # Header stops at first ';'
            result = extract_location_and_name(header)
            if result:
                results.append(result)

        return results

def export_to_bed(result_list, output_bed_path):
    with open(output_bed_path, 'w') as bed_file:
        for result in result_list:
            # Write each result as a BED file line: chromosome, start, end, name
            bed_file.write(f"chr{result['chromosome']}\t{result['start']}\t{result['end']}\t{result['name']}\n")


# Usage example
file_path = "/Users/jillianness/Desktop/sorting_cannavo_data/final_merged_libraries/dmel-all-transposon-r6.59.fasta"  # Replace with the path to your input file

output_bed_path = '/Users/jillianness/Desktop/sorting_cannavo_data/final_merged_libraries/FB_TE.bed'
result_list = process_headers_from_file(file_path)
export_to_bed(result_list, output_bed_path)
print(f"Data successfully exported to {output_bed_path}")

# Print the results
for result in result_list:
    print(f"Chromosome: {result['chromosome']}, Start: {result['start']}, End: {result['end']}, Name: {result['name']}")
