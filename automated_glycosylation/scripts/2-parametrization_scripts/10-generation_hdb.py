#!/usr/bin/env python3
"""
Script to extract matches between residues from .rtp and .hdb files
and generate a new .hdb file with the corresponding data.
"""

import re
import sys
import argparse


def extrair_correspondencias_rtp(arquivo_rtp):
    """
    Extracts matches between PDB codes and CHARMM names from the .rtp file.
    Keeps all entries, even if duplicates exist.
    """
    correspondencias = []

    with open(arquivo_rtp, 'r') as f:
        conteudo = f.read()

    padrao = r'\[ (\w{3}) \]\s*;.*?Original CHARMM name:\s*(\w+)'
    matches = re.findall(padrao, conteudo, re.DOTALL)

    for codigo_pdb, nome_charmm in matches:
        correspondencias.append((codigo_pdb, nome_charmm.strip()))

    return correspondencias


def ler_arquivo_hdb(arquivo_hdb):
    """
    Reads the .hdb file and stores data indexed by CHARMM residue name.
    Compatible with alphanumeric residue names (e.g. ANE5AC).
    """
    dados_hdb = {}

    with open(arquivo_hdb, 'r') as f:
        linhas = f.readlines()

    i = 0
    total_linhas = len(linhas)

    while i < total_linhas:
        linha = linhas[i].strip()

        # Ignore comments and empty lines
        if not linha or linha.startswith(';') or linha.startswith('#'):
            i += 1
            continue

        partes = linha.split()

        # Accept alphanumeric residue names
        if partes and re.match(r'^[A-Za-z0-9]+$', partes[0]):
            nome_charmm = partes[0]

            # Check if the number of entries is on the same line
            if len(partes) > 1 and partes[1].isdigit():
                num_dados = int(partes[1])
                i += 1
            else:
                i += 1
                continue

            dados_residuo = [linha]
            dados_coletados = 0

            while i < total_linhas and dados_coletados < num_dados:
                linha_atual = linhas[i].strip()

                if linha_atual and not linha_atual.startswith(';') and not linha_atual.startswith('#'):
                    partes_linha = linha_atual.split()
                    if partes_linha and partes_linha[0].isdigit():
                        dados_residuo.append(linha_atual)
                        dados_coletados += 1

                i += 1

            if dados_coletados == num_dados:
                dados_hdb[nome_charmm] = dados_residuo
                print(f"  Read: {nome_charmm} with {num_dados} entries")
            else:
                print(
                    f"  WARNING: {nome_charmm} - Expected {num_dados} entries, "
                    f"found {dados_coletados}"
                )
        else:
            i += 1

    return dados_hdb


def gerar_hdb_correspondente(correspondencias, dados_hdb, arquivo_saida):
    """
    Generates the new .hdb file using PDB codes instead of CHARMM names.
    Keeps duplicates and uses all available data.
    """
    with open(arquivo_saida, 'w') as f:
        contador = 0
        contador_nao_encontrados = 0

        for codigo_pdb, nome_charmm in correspondencias:
            if nome_charmm in dados_hdb:
                dados = dados_hdb[nome_charmm]

                primeira_linha = dados[0]
                partes = primeira_linha.split()

                # Replace CHARMM name with PDB code
                resto = ' '.join(partes[1:])
                f.write(f"{codigo_pdb} {resto}\n")

                for linha in dados[1:]:
                    f.write(linha + "\n")

                f.write("\n")

                contador += 1
                print(f"  Processed: {codigo_pdb} -> {nome_charmm}")
            else:
                print(f"  WARNING: Data not found for {codigo_pdb} -> {nome_charmm}")
                contador_nao_encontrados += 1

        return contador, contador_nao_encontrados


def main():
    parser = argparse.ArgumentParser(
        description="Generate a modified .hdb file based on residue mappings from an .rtp file."
    )
    parser.add_argument("rtp", help="Input .rtp file")
    parser.add_argument("hdb", help="Input .hdb file")
    parser.add_argument(
        "-o", "--output", default="carb_modified.hdb",
        help="Output .hdb file (default: carb_modified.hdb)"
    )

    args = parser.parse_args()

    arquivo_rtp = args.rtp
    arquivo_hdb = args.hdb
    arquivo_saida = args.output

    try:
        print("=" * 60)
        print("Extracting matches from .rtp file...")
        print("=" * 60)

        correspondencias = extrair_correspondencias_rtp(arquivo_rtp)

        print(f"\nFound {len(correspondencias)} matches:")
        for i, (pdb, chm) in enumerate(correspondencias, 1):
            print(f"  {i:3d}. {pdb} -> {chm}")

        print("\n" + "=" * 60)
        print("Reading .hdb file data...")
        print("=" * 60)

        dados_hdb = ler_arquivo_hdb(arquivo_hdb)

        print(f"\nFound {len(dados_hdb)} residues in the .hdb file")

        print("\n" + "=" * 60)
        print(f"Generating output file {arquivo_saida}...")
        print("=" * 60)

        ok, nao_ok = gerar_hdb_correspondente(
            correspondencias, dados_hdb, arquivo_saida
        )

        print("\n" + "=" * 60)
        print("SUMMARY")
        print("=" * 60)
        print(f"- Total entries in .rtp: {len(correspondencias)}")
        print(f"- Successfully processed: {ok}")
        print(f"- Not found in .hdb: {nao_ok}")

        print(f"\nFile '{arquivo_saida}' generated successfully.")

    except FileNotFoundError as e:
        print(f"\nERROR: File not found - {e}")
        sys.exit(1)
    except Exception:
        print("\nUnexpected ERROR:")
        import traceback
        traceback.print_exc()
        sys.exit(1)


if __name__ == "__main__":
    main()

