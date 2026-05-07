import re
import pandas as pd
import argparse

# Configurar argumentos
parser = argparse.ArgumentParser(description="Correct Caselino table sequences")
parser.add_argument("--input", required=True, help="Input TSV file")
parser.add_argument("--output", required=True, help="Output TSV file")
args = parser.parse_args()

input_file = args.input
output_file = args.output

# Read TSV file
df = pd.read_csv(input_file, sep="\t")

def fix_sequence(site, sequence):
    if pd.isna(site) or pd.isna(sequence):
        return sequence

    # Extrair número do site (ex: N26 -> 26, T330 -> 330, S332 -> 332)
    site_match = re.search(r"[A-Z](\d+)", str(site))
    if not site_match:
        return sequence

    site_number = site_match.group(1)

    # Corrigir qualquer sufixo PROx- número final
    fixed_sequence = re.sub(
        r"(PRO[A-Z]*-)\d+$",
        lambda m: m.group(1) + site_number,
        str(sequence)
    )

    return fixed_sequence

# Aplicar a correção
df["sequence"] = df.apply(lambda row: fix_sequence(row["site"], row["sequence"]), axis=1)

# Salvar arquivo corrigido
df.to_csv(output_file, sep="\t", index=False)

print(f"Corrected file saved as: {output_file}")

