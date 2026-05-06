#!/usr/bin/env python3
"""
Script for residue renumbering and labeling in PDB output files from Glycosylator.

Functionalities:
1. Renumbers protein (ATOM) and carbohydrate (HETATM) residues
2. Maintains comments and everything in English
3. Provides options for chain labeling
4. Preserves connectivity (CONECT lines)
5. Identifies glycan blocks based on sequence in the file

Author: Anacleto SIlva de Souza
"""

import sys
import os
from typing import Dict, Tuple, List, Optional, Set

def parse_pdb_line(line: str) -> Optional[Dict]:
    """Parse a PDB line and return a dictionary with the fields."""
    if len(line) < 6:
        return None
    
    record = line[0:6].strip()
    
    if record in ["ATOM", "HETATM"]:
        try:
            return {
                "record": record,
                "serial": int(line[6:11].strip()) if line[6:11].strip() else 0,
                "name": line[12:16].strip(),
                "altLoc": line[16:17].strip(),
                "resName": line[17:20].strip(),
                "chain": line[21:22].strip(),
                "resSeq": int(line[22:26].strip()) if line[22:26].strip() else 0,
                "iCode": line[26:27].strip(),
                "x": line[30:38].strip(),
                "y": line[38:46].strip(),
                "z": line[46:54].strip(),
                "occupancy": line[54:60].strip(),
                "tempFactor": line[60:66].strip(),
                "element": line[76:78].strip(),
                "charge": line[78:80].strip(),
                "line": line.rstrip('\n')
            }
        except (ValueError, IndexError):
            return None
    elif line.startswith("CONECT"):
        return {"record": "CONECT", "line": line.rstrip('\n')}
    else:
        return {"record": "OTHER", "line": line.rstrip('\n')}

def format_pdb_line(data: Dict) -> str:
    """Format a data dictionary into a PDB line."""
    if data["record"] in ["ATOM", "HETATM"]:
        serial = str(data.get("serial", 1))
        name = data.get("name", "")
        altLoc = data.get("altLoc", "")
        resName = data.get("resName", "")
        chain = data.get("chain", "")
        resSeq = str(data.get("resSeq", 1))
        iCode = data.get("iCode", "")
        x = data.get("x", "0.000")
        y = data.get("y", "0.000")
        z = data.get("z", "0.000")
        occupancy = data.get("occupancy", "1.00")
        tempFactor = data.get("tempFactor", "0.00")
        element = data.get("element", "")
        charge = data.get("charge", "")
        
        # Format fields in standard PDB format
        line = f"{data['record']:6s}{serial:>5s} {name:<4s}{altLoc:1s}{resName:>3s} {chain:1s}{resSeq:>4s}{iCode:1s}   {x:>8s}{y:>8s}{z:>8s}{occupancy:>6s}{tempFactor:>6s}          {element:>2s}{charge:>2s}"
        return line
    else:
        return data["line"]

def analyze_structure(pdb_lines: List[str]) -> Tuple[Dict, Dict, List[Dict], List[Tuple]]:
    """Analyzes the PDB file structure to identify patterns."""
    chains_info = {}
    atom_data = []
    structure_blocks = []  # List of blocks (type, chain, start_line, end_line)
    
    current_block_type = None
    current_block_chain = None
    block_start = 0
    
    for i, line in enumerate(pdb_lines):
        data = parse_pdb_line(line)
        
        if data and data["record"] in ["ATOM", "HETATM"]:
            atom_data.append(data)
            
            chain = data["chain"]
            res_seq = data["resSeq"]
            res_name = data["resName"]
            record_type = data["record"]
            
            if chain not in chains_info:
                chains_info[chain] = {
                    "residues": set(),
                    "min_res": float('inf'),
                    "max_res": float('-inf'),
                    "types": set(),
                    "first_record": record_type
                }
            
            chains_info[chain]["residues"].add(res_seq)
            chains_info[chain]["min_res"] = min(chains_info[chain]["min_res"], res_seq)
            chains_info[chain]["max_res"] = max(chains_info[chain]["max_res"], res_seq)
            chains_info[chain]["types"].add(record_type)
            
            # Detects block changes
            block_type = "protein" if record_type == "ATOM" else "glycan"
            
            if current_block_type is None:
                current_block_type = block_type
                current_block_chain = chain
                block_start = i
            elif block_type != current_block_type or chain != current_block_chain:
                # Finishes previous block
                structure_blocks.append((current_block_type, current_block_chain, block_start, i-1))
                # Starts new block
                current_block_type = block_type
                current_block_chain = chain
                block_start = i
    
    # Adds last block
    if current_block_type is not None:
        structure_blocks.append((current_block_type, current_block_chain, block_start, len(pdb_lines)-1))
    
    # Classifies chains as protein or glycan
    chain_types = {}
    for chain in chains_info:
        if "ATOM" in chains_info[chain]["types"]:
            chain_types[chain] = "protein"
        else:
            chain_types[chain] = "glycan"
    
    return chains_info, chain_types, atom_data, structure_blocks

def get_user_options() -> Dict:
    """Requests options from the user."""
    print("=" * 60)
    print("PDB Residue Renumbering and Relabeling Tool")
    print("=" * 60)
    
    options = {}
    
    # Fixed option for glycan chain labeling - Cadeia A para proteina e carboidrato, cadeia B para proteina e carboidrato
    print("\nChain labeling for glycans: Fixed mapping (A->A, B->B, etc.)")
    options["glycan_labeling"] = 2  # Fixed to option 2: A->A, B->B
    
    # Fixed option for protein renumbering - manter numero original da proteina
    print("\nProtein residue renumbering: Keep original numbering")
    options["protein_renumbering"] = 4  # New option to keep original numbers
    
    # Fixed option for glycan renumbering - cada bloco começa com 1
    print("\nGlycan residue renumbering: Each block starts from 1")
    options["glycan_renumbering"] = 4  # New option for block restart
    
    return options

def generate_glycan_chain_mapping(protein_chains: List[str], options: Dict) -> Dict:
    """Generates glycan chain mapping based on user options."""
    mapping = {}
    
    # Fixed to option 2: A->A, B->B, C->C, etc.
    for prot_chain in sorted(protein_chains):
        mapping[prot_chain] = prot_chain
    
    return mapping

def identify_glycan_blocks(structure_blocks: List[Tuple], chain_types: Dict) -> List[Tuple]:
    """Identifies glycan blocks and associates them with preceding proteins."""
    glycan_blocks = []
    current_protein_chain = None
    
    for i, (block_type, chain, start, end) in enumerate(structure_blocks):
        if block_type == "protein":
            current_protein_chain = chain
        elif block_type == "glycan" and current_protein_chain:
            # Finds all consecutive glycan blocks
            glycan_start = i
            glycan_end = i
            
            # Checks if there are more consecutive glycan blocks
            j = i + 1
            while j < len(structure_blocks) and structure_blocks[j][0] == "glycan":
                glycan_end = j
                j += 1
            
            glycan_blocks.append({
                "protein_chain": current_protein_chain,
                "glycan_chain": chain,
                "blocks": structure_blocks[glycan_start:glycan_end+1]
            })
    
    return glycan_blocks

def process_pdb_file(input_file: str, output_file: str, options: Dict) -> None:
    """Processes the PDB file with the provided options."""
    print(f"\nProcessing file: {input_file}")
    
    # Reads the file
    with open(input_file, 'r') as f:
        lines = f.readlines()
    
    # Analyzes the structure
    chains_info, chain_types, atom_data, structure_blocks = analyze_structure(lines)
    
    print(f"\nFound {len(chains_info)} chains:")
    for chain in sorted(chains_info.keys()):
        chain_type = chain_types.get(chain, "unknown")
        residues = chains_info[chain]["residues"]
        print(f"  Chain {chain}: {chain_type}, {len(residues)} residues "
              f"(range: {min(residues)}-{max(residues)})")
    
    # Identifies protein chains
    protein_chains = [c for c in chain_types if chain_types[c] == "protein"]
    print(f"\nProtein chains: {', '.join(sorted(protein_chains))}")
    
    # Identifies glycan blocks
    glycan_blocks = identify_glycan_blocks(structure_blocks, chain_types)
    
    print(f"\nFound {len(glycan_blocks)} glycan blocks:")
    for i, block_info in enumerate(glycan_blocks):
        total_residues = 0
        for _, _, start, end in block_info["blocks"]:
            # Counts residues in this block
            residues_in_block = set()
            for j in range(start, end + 1):
                data = parse_pdb_line(lines[j])
                if data and data["record"] == "HETATM":
                    residues_in_block.add(data["resSeq"])
            total_residues += len(residues_in_block)
        
        print(f"  Block {i+1}: Protein chain {block_info['protein_chain']} -> "
              f"{len(block_info['blocks'])} glycan sub-blocks, {total_residues} total residues")
    
    # Generates glycan chain mapping
    glycan_chain_map = generate_glycan_chain_mapping(protein_chains, options)
    print(f"\nGlycan chain mapping:")
    for prot, glycan in sorted(glycan_chain_map.items()):
        print(f"  Protein {prot} -> Glycan {glycan}")
    
    # Processes the file
    output_lines = []
    conect_lines = []
    
    # Dictionaries for renumbering mapping
    residue_mapping = {}  # (old_chain, old_residue) -> (new_chain, new_residue)
    
    # Counters for renumbering - each glycan block starts at 1
    glycan_block_counter = 1
    
    # Process protein chains - keep original numbering
    for chain in sorted(chains_info.keys()):
        if chain_types[chain] != "protein":
            continue
            
        residues = sorted(chains_info[chain]["residues"])
        
        # Keep original residue numbers for protein
        for old_res in residues:
            residue_mapping[(chain, old_res)] = (chain, old_res)
    
    # Process glycan blocks - each block starts at 1
    for block_info in glycan_blocks:
        protein_chain = block_info["protein_chain"]
        old_glycan_chain = block_info["glycan_chain"]
        
        if protein_chain in glycan_chain_map:
            new_glycan_chain = glycan_chain_map[protein_chain]
            
            # Reset counter for each glycan block
            current_glycan_residue = 1
            
            # For each glycan block
            for _, _, start, end in block_info["blocks"]:
                # Collects all residues in this block
                residues_in_block = set()
                for j in range(start, end + 1):
                    data = parse_pdb_line(lines[j])
                    if data and data["record"] == "HETATM":
                        residues_in_block.add(data["resSeq"])
                
                # Sorts residues
                sorted_residues = sorted(residues_in_block)
                
                # Creates mapping for these glycan residues - starts at 1 for each block
                new_residue = 1
                for old_res in sorted_residues:
                    residue_mapping[(old_glycan_chain, old_res)] = (new_glycan_chain, new_residue)
                    new_residue += 1
    
    # Processes each line of the file
    atom_counter = 1
    for line in lines:
        data = parse_pdb_line(line)
        
        if not data:
            output_lines.append(line.rstrip('\n'))
            continue
        
        if data["record"] in ["ATOM", "HETATM"]:
            old_chain = data["chain"]
            old_res = data["resSeq"]
            
            # Applies mapping
            if (old_chain, old_res) in residue_mapping:
                new_chain, new_res = residue_mapping[(old_chain, old_res)]
                data["chain"] = new_chain
                data["resSeq"] = new_res
            
            # Renumbers atoms sequentially
            data["serial"] = atom_counter
            atom_counter += 1
            
            output_lines.append(format_pdb_line(data))
            
        elif data["record"] == "CONECT":
            # Saves CONECT lines for later processing
            conect_lines.append(data["line"])
        else:
            output_lines.append(data["line"])
    
    # Processes CONECT lines
    for conect_line in conect_lines:
        output_lines.append(conect_line)
    
    # Writes the output file
    with open(output_file, 'w') as f:
        for line in output_lines:
            f.write(line + '\n')
    
    print(f"\nSuccessfully processed {len(lines)} lines.")
    print(f"Output written to: {output_file}")
    
    # Summary report
    print("\n" + "=" * 60)
    print("SUMMARY")
    print("=" * 60)
    
    # Counts residues per chain in the output file
    output_chains = {}
    for line in output_lines:
        if len(line) >= 22 and line[0:6].strip() in ["ATOM", "HETATM"]:
            chain = line[21:22].strip() or " "
            res_seq = int(line[22:26].strip()) if line[22:26].strip() else 0
            if chain not in output_chains:
                output_chains[chain] = set()
            output_chains[chain].add(res_seq)
    
    for chain in sorted(output_chains.keys()):
        residues = output_chains[chain]
        print(f"Chain {chain or ' '}: {len(residues)} residues "
              f"(range: {min(residues)}-{max(residues)})")

def main():
    """Main function."""
    if len(sys.argv) < 2:
        print("Usage: python pdb_renumber.py input.pdb [output.pdb]")
        print("If output file is not specified, '_renumbered' will be appended to input filename.")
        sys.exit(1)
    
    input_file = sys.argv[1]
    
    if not os.path.exists(input_file):
        print(f"Error: Input file '{input_file}' not found.")
        sys.exit(1)
    
    if len(sys.argv) >= 3:
        output_file = sys.argv[2]
    else:
        base, ext = os.path.splitext(input_file)
        output_file = f"{base}_renumbered{ext}"
    
    # Gets user options (now fixed, no user input needed)
    options = get_user_options()
    
    # Processes the file
    process_pdb_file(input_file, output_file, options)
    
    print("\nDone!")

if __name__ == "__main__":
    main()
