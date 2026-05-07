import re
import argparse

def main():
    parser = argparse.ArgumentParser(description="Remove redundant residue blocks from a total_carb_unique.rtp file.")
    parser.add_argument("--input", "-i", required=True, help="Input total_carb_unique.rtp file")
    parser.add_argument("--output", "-o", required=True, help="Output file with redundant blocks removed")
    args = parser.parse_args()

    input_file = args.input
    output_file = args.output

    # Dictionary to store unique residue blocks
    unique_blocks = {}

    # Regex to detect block start: [ RESNAME ]
    block_start_pattern = re.compile(r"^\[\s*(\w+)\s*\]", re.MULTILINE)

    with open(input_file, "r") as f:
        content = f.read()

    # Find all block starts
    matches = list(block_start_pattern.finditer(content))

    for i, match in enumerate(matches):
        resname = match.group(1)
        start = match.start()
        end = matches[i + 1].start() if i + 1 < len(matches) else len(content)
        block = content[start:end].strip()

        # Add block if it's the first occurrence of the residue
        if resname not in unique_blocks:
            unique_blocks[resname] = block

    # Write unique blocks to output file, separated by two blank lines
    with open(output_file, "w") as f:
        f.write("\n\n".join(unique_blocks.values()) + "\n")

    print(f"Redundant blocks removed. Output saved to: {output_file}")

if __name__ == "__main__":
    main()

