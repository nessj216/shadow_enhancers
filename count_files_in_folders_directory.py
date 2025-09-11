import os

def count_txt_files_in_folder(folder_path):
    """Return the number of .txt files in the given folder (non-recursive)."""
    count = 0
    for entry in os.listdir(folder_path):
        full_path = os.path.join(folder_path, entry)
        if os.path.isfile(full_path) and entry.lower().endswith('.txt'):
            count += 1
    return count

def main():
    # 1) Specify the directory path you want to analyze.
    #    Example:  target_dir = r"C:\Users\Example\Documents"
    #    On Linux/Mac: target_dir = "/home/user/some_directory"
    target_dir = r"/Users/jillianness/Desktop/Mouse_SE_birth_analysis"

    if not os.path.isdir(target_dir):
        print(f"Error: The path '{target_dir}' is not a valid directory.")
        return  # Stop execution if the path doesn't exist

    # 2) Get all top-level folders in the specified directory
    top_level_folders = [
        d for d in os.listdir(target_dir)
        if os.path.isdir(os.path.join(target_dir, d))
    ]

    # 3) Print header row
    print(f"{'Folder':<40}  {'Total Subfolders':>15}  {'Subfolders w/ 1 txt':>20}")
    print("-" * 80)

    # 4) For each top-level folder, walk through its subfolders recursively
    for folder in top_level_folders:
        folder_path = os.path.join(target_dir, folder)

        total_subfolders = 0
        single_txt_subfolders = 0

        for root, dirs, files in os.walk(folder_path):
            # Skip the folder_path itself (only count its subfolders)
            if root == folder_path:
                continue

            total_subfolders += 1

            # Count .txt files in this subfolder (non-recursive)
            txt_count = count_txt_files_in_folder(root)
            if txt_count == 1:
                single_txt_subfolders += 1

        # 5) Print the results for the current folder
        print(f"{folder:<40}  {total_subfolders:>15}  {single_txt_subfolders:>20}")

# 6) Execute 'main' if the script is run
if __name__ == "__main__":
    main()
