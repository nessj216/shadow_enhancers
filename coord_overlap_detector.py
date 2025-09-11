import os

# Function to parse coordinates from a header line
def parse_coordinates(header):
    try:
        # Assuming the format is ">chr3R:start-end"
        parts = header.split(b':')
        start_end = parts[1].split(b'-')
        start = int(start_end[0])
        end = int(start_end[1])
        return start, end
    except (IndexError, ValueError):
        return None, None

# Function to check if two sets of coordinates overlap
def coordinates_overlap(coords1, coords2):
    start1, end1 = coords1
    start2, end2 = coords2
    return end1 >= start2 and end2 >= start1

# Directory containing the files to check
directory_path = "//Users/jillianness/Desktop/fly_blast/CannavoFlyDataSet/CurbioOutput/FBgn0038056"

# List all files in the directory
files = os.listdir(directory_path)

# Dictionary to store file paths and their coordinates
file_coordinates = {}

# Parse coordinates for each file
for filename in files:
    file_path = os.path.join(directory_path, filename)
    if os.path.isfile(file_path):  # Check if it's a file, not a directory
        with open(file_path, 'rb') as file:
            header = file.readline().strip()
            start, end = parse_coordinates(header)
            if start is not None and end is not None:
                file_coordinates[file_path] = (start, end)

# Check for overlaps
for file1, coords1 in file_coordinates.items():
    for file2, coords2 in file_coordinates.items():
        if file1 != file2 and coordinates_overlap(coords1, coords2):
            print(f"Overlap detected between {file1} and {file2}")
    print(coords1, coords2)