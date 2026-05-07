#!/usr/bin/env python3
import re
import argparse
from collections import defaultdict

def parse_rtp(filename):
    residues = {}
    with open(filename) as f:
        lines = f.readlines()

    current = None
    section = None

    for line in lines:
        line = line.strip()
        if not line or line.startswith(";"):
            continue

        m = re.match(r"\[\s*(\w+)\s*\]", line)
        if m:
            key = m.group(1)
            # nome de resíduo
            if key.isupper() and len(key) <= 5:
                current = key
                residues[current] = {"atoms": [], "bonds": []}
                section = None
            else:
                section = key
            continue

        if section == "atoms":
            atom = line.split()[0]
            residues[current]["atoms"].append(atom)

        elif section == "bonds":
            a, b = line.split()[:2]
            residues[current]["bonds"].append((a, b))

    return residues


def build_connectivity(bonds):
    conn = defaultdict(list)
    for a, b in bonds:
        conn[a].append(b)
        conn[b].append(a)
    return conn


def classify_h(parent, neighbors):
    # regras compatíveis com carboidratos CHARMM
    if parent.startswith("O"):
        return 1, 2        # OH
    if len(neighbors) == 3:
        return 3, 4        # CH3
    if len(neighbors) == 2:
        return 2, 6        # CH2
    return 1, 5            # CH


def write_hdb(residues, outfile):
    with open(outfile, "w") as out:
        for res, data in residues.items():
            conn = build_connectivity(data["bonds"])
            hydrogens = [a for a in data["atoms"] if a.startswith("H")]

            out.write(f"{res:<6} {len(hydrogens)}\n")

            for h in hydrogens:
                parent = conn[h][0]
                neighbors = [n for n in conn[parent] if n != h]

                htype, nat = classify_h(parent, neighbors)
                atoms = [h, parent] + neighbors
                atoms = atoms[:nat]

                out.write(
                    f"{htype:1d} {nat:7d} " +
                    " ".join(f"{a:<6}" for a in atoms) + "\n"
                )


def main():
    parser = argparse.ArgumentParser(
        description="Gera arquivo .hdb a partir de um arquivo .rtp (por resíduo)"
    )
    parser.add_argument(
        "-i", "--input",
        required=True,
        help="Arquivo RTP de entrada (ex: carb_redundance_removed.rtp)"
    )
    parser.add_argument(
        "-o", "--output",
        required=True,
        help="Arquivo HDB de saída (ex: carb_redundance_removed.hdb)"
    )

    args = parser.parse_args()

    residues = parse_rtp(args.input)
    write_hdb(residues, args.output)


if __name__ == "__main__":
    main()

