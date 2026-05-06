import os
import re
import argparse

# Function to process each file
def process_file(file_path, marker_pattern):
    with open(file_path, "r") as f:
        content = f.read()
    # Ignore everything after the marker
    split_content = marker_pattern.split(content, maxsplit=1)
    return split_content[0].rstrip()  # remove trailing spaces

def main():
    parser = argparse.ArgumentParser(description="Combine carb_unique.rtp files from multiple directories into a single file.")
    parser.add_argument("--input", "-i", required=True, help="Root directory containing subdirectories with carb_unique.rtp files.")
    parser.add_argument("--output", "-o", required=True, help="Output file path for the combined carb_unique.rtp.")
    args = parser.parse_args()

    root_dir = args.input
    output_file = args.output

    # Marker pattern to ignore content from
    marker_pattern = re.compile(
        r";\s*={3,}\s*\n;\s*ORIGINAL CHARMM CARBOHYDRATE RESIDUES FOR REFERENCE",
        re.IGNORECASE
    )

    # List all subdirectories matching *_*
    folders = [f for f in os.listdir(root_dir) if os.path.isdir(os.path.join(root_dir, f))]
    folders.sort()  # Simple alphabetical sorting

    all_contents = []

    for folder in folders:
        rtp_file = os.path.join(root_dir, folder, "carb_unique.rtp")
        if os.path.exists(rtp_file):
            processed_content = process_file(rtp_file, marker_pattern)
            all_contents.append(processed_content)

    # Combine everything with two blank lines between each block
    final_content = "\n\n" + "\n\n".join(all_contents) + "\n"

    # Save to output file
    with open(output_file, "w") as f:
        f.write(final_content)

    print(f"All files have been combined into: {output_file}")

if __name__ == "__main__":
    main()

