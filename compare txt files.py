with open('/Users/jillianness/Desktop/single_enhancer_cannavo/cannavo_shadow_dataset4thcol.txt', 'r') as file1:
    lines_file1 = set(file1.readlines())

with open('/Users/jillianness/Desktop/single_enhancer_cannavo/cannavo_all4thcol.txt', 'r') as file2:
    lines_file2 = set(file2.readlines())

matching_lines = lines_file1.intersection(lines_file2)
non_matching_lines_file1 = lines_file1 - matching_lines
non_matching_lines_file2 = lines_file2 - matching_lines

print("Number of matching lines:", len(matching_lines))

print('file1 unique', len(non_matching_lines_file1))
print('file2 unique', len(non_matching_lines_file2))

print(len(matching_lines))
'''print("\nLines from file1 that don't match:")
for line in non_matching_lines_file1:
    print(line.strip())

print("\nLines from file2 that don't match:")
for line in non_matching_lines_file2:
    print(line.strip())'''

with open('/Users/jillianness/Desktop/single_enhancer_cannavo/non_matching_lines_file2.txt', 'w') as output_file:
    for line in non_matching_lines_file2:
        output_file.write(line)