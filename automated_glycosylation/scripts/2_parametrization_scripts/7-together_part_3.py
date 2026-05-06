import argparse

def main():
    parser = argparse.ArgumentParser(description="Append the content of one file to another with two blank lines before.")
    parser.add_argument("--input", "-i", required=True, help="Input file to append (e.g., carb_unique_redundance_removed.rtp)")
    parser.add_argument("--output", "-o", required=True, help="Output file to append to (e.g., carb_unique_total.rtp)")
    args = parser.parse_args()

    input_file = args.input
    output_file = args.output

    # Read content from the input file
    with open(input_file, "r") as f:
        content_to_append = f.read()

    # Append content to the output file, with two blank lines before
    with open(output_file, "a") as f:
        f.write("\n\n" + content_to_append + "\n")

    print(f"Content from '{input_file}' appended to '{output_file}' with two blank lines before.")

if __name__ == "__main__":
    main()

