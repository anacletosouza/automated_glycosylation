import re
import argparse

# --- Command-line argument parsing ---
parser = argparse.ArgumentParser(description="Filter unique RTP blocks by letter.")
parser.add_argument("input_file", help="Path to the input RTP file")
parser.add_argument("output_file", help="Path to the output RTP file")
parser.add_argument("letter", choices=["A","B","C","D", "E", "F", "G", "H", "I", "J", "K", "L", "M", "N", "O", "P", "Q", "R", "S", "T", "U", "V", "X", "Z", "W", "Y", "Z"], help="Letter to filter.")

args = parser.parse_args()
input_file = args.input_file
output_file = args.output_file
letter = args.letter

# --- Regex to find blocks like [NAME_D1], [NAME_E2], etc. ---
block_start_pattern = re.compile(rf"^\[\s*(\w+_{letter}\d+)\s*\]", re.MULTILINE)

# --- Function to process the file ---
def process_rtp_file(input_path, output_path, letter):
    with open(input_path, "r") as f:
        content = f.read()

    # Find all blocks
    blocks = re.split(r"(?=^\[\s*\w+_" + letter + r"\d+\s*\])", content, flags=re.MULTILINE)

    unique_names = set()
    selected_blocks = []

    for block in blocks:
        match = re.match(r"^\[\s*(\w+)_" + letter + r"(\d+)\s*\]", block)
        if match:
            name = match.group(1)  # extract NAME
            if name not in unique_names:
                unique_names.add(name)
                # Rename the block to [NAME] and keep the rest of the content
                new_block = re.sub(r"^\[\s*\w+_" + letter + r"\d+\s*\]", f"[ {name} ]", block, count=1)
                selected_blocks.append(new_block.strip() + "\n")

    # Write final file
    with open(output_path, "w") as f:
        f.write("\n\n".join(selected_blocks))

    print(f"{len(selected_blocks)} unique blocks saved in '{output_path}'.")

# --- Run the function ---
process_rtp_file(input_file, output_file, letter)

