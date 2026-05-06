#!/usr/bin/env python3
"""
Script for residue renumbering and labeling in PDB output files from Glycosylator.

Functionalities:
1. Renumbers protein (ATOM) and carbohydrate (HETATM) residues
2. Maintains comments and everything in English
3. Provides options for chain labeling
4. Preserves connectivity (CONECT lines)
5. Identifies glycan blocks based on sequence in the file

Author: Anacleto Souza
"""

import sys
import os
import argparse
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

def identify_glycan_blocks_by_residue_number(pdb_lines: List[str], residue_initiation: Set[str]) -> List[Dict]:
    """
    Identifies glycan blocks based on residue number changes and initiation residues.
    
    Scans line by line. When it finds a residue initiation (NDG, A2G, or user-defined),
    starts a new block and resets the counter to start_number_carb.
    Continues numbering sequentially until the next initiation residue.
    """
    blocks = []
    current_block = None
    last_resSeq = None
    
    for line_num, line in enumerate(pdb_lines):
        data = parse_pdb_line(line)
        
        if not data or data["record"] not in ["ATOM", "HETATM"]:
            continue
        
        res_name = data["resName"]
        res_seq = data["resSeq"]
        
        # Check if this is an initiation residue
        is_initiation = res_name in residue_initiation
        
        if is_initiation:
            # Start a new block
            if current_block is not None:
                blocks.append(current_block)
            current_block = {
                "start_line": line_num,
                "residues": [],
                "initiation_residue": res_name,
                "initiation_seq": res_seq
            }
            current_block["residues"].append((line_num, data))
            last_resSeq = res_seq
        elif current_block is not None:
            # Check if this is a new residue (different resSeq)
            if res_seq != last_resSeq:
                # Same block, different residue - will be assigned sequential numbers
                current_block["residues"].append((line_num, data))
                last_resSeq = res_seq
            else:
                # Same residue, same block - add to current residue
                current_block["residues"].append((line_num, data))
    
    # Add the last block
    if current_block is not None:
        blocks.append(current_block)
    
    return blocks

def process_pdb_file(input_file: str, output_file: str, 
                     start_number_protein: int = 1,
                     start_number_carb: int = 1,
                     residue_carb_initiation: str = "NDG,A2G") -> None:
    """Processes the PDB file with the standard chain pattern."""
    
    print(f"\nProcessing file: {input_file}")
    print(f"Protein start number: {start_number_protein}")
    print(f"Carbohydrate start number: {start_number_carb}")
    print(f"Carbohydrate initiation residues: {residue_carb_initiation}")
    
    # Parse initiation residues
    init_residues = set(r.strip() for r in residue_carb_initiation.split(','))
    
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
    protein_chains = sorted([c for c in chain_types if chain_types[c] == "protein"])
    print(f"\nProtein chains: {', '.join(protein_chains)}")
    
    # Generate chain mapping: protein chains get A, B, C, etc.
    # Carbohydrate chains get the same letters as their associated protein chains
    chain_mapping = {}
    protein_chain_letters = {}
    
    # Assign letters to protein chains
    for idx, prot_chain in enumerate(protein_chains):
        chain_letter = chr(ord('A') + idx)
        protein_chain_letters[prot_chain] = chain_letter
        chain_mapping[prot_chain] = chain_letter
    
    print(f"\nChain mapping:")
    for prot_chain, letter in protein_chain_letters.items():
        print(f"  Protein chain {prot_chain} -> {letter}")
    
    # Identify glycan blocks
    glycan_blocks = identify_glycan_blocks_by_residue_number(lines, init_residues)
    
    print(f"\nFound {len(glycan_blocks)} glycan blocks:")
    for i, block in enumerate(glycan_blocks):
        # Count unique residues in this block
        unique_residues = set()
        for _, data in block["residues"]:
            unique_residues.add(data["resSeq"])
        print(f"  Block {i+1}: Initiated with {block['initiation_residue']} "
              f"(original seq {block['initiation_seq']}), {len(unique_residues)} residues")
    
    # Create mapping for protein residues (renumber sequentially from start_number_protein)
    protein_residue_mapping = {}  # (old_chain, old_residue) -> (new_chain, new_residue)
    
    # Process each protein chain
    for prot_chain in protein_chains:
        new_chain = protein_chain_letters[prot_chain]
        residues = sorted(chains_info[prot_chain]["residues"])
        current_res = start_number_protein
        
        for old_res in residues:
            protein_residue_mapping[(prot_chain, old_res)] = (new_chain, current_res)
            current_res += 1
    
    # Create mapping for glycan residues
    glycan_residue_mapping = {}  # (old_chain, old_residue) -> (new_chain, new_residue)
    
    # For each glycan block, assign sequential numbers starting from start_number_carb
    block_counter = 0
    for block in glycan_blocks:
        block_counter += 1
        current_res = start_number_carb
        last_res_seq = None
        
        # The chain for this glycan block is determined by the associated protein chain
        # Since we don't track which protein each glycan is attached to in this simplified version,
        # we need to infer it. For now, we'll use chain A for all glycans,
        # but we can look at the original chain to determine.
        # Better approach: find the protein chain that precedes this glycan block in the file
        
        # Find which protein chain this glycan is associated with
        # Look backwards for the nearest protein chain
        associated_protein_chain = None
        block_start_line = block["start_line"]
        
        for i in range(block_start_line - 1, -1, -1):
            data = parse_pdb_line(lines[i])
            if data and data["record"] == "ATOM":
                associated_protein_chain = data["chain"]
                break
        
        # If we found an associated protein chain, use its mapped chain letter
        if associated_protein_chain and associated_protein_chain in protein_chain_letters:
            new_glycan_chain = protein_chain_letters[associated_protein_chain]
        else:
            # Default to the first protein chain's letter
            new_glycan_chain = protein_chain_letters.get(protein_chains[0], "A") if protein_chains else "A"
        
        # Map residues in this block
        for _, data in block["residues"]:
            old_chain = data["chain"]
            old_res = data["resSeq"]
            
            # Check if this is a new residue
            if old_res != last_res_seq:
                if last_res_seq is not None:
                    current_res += 1
                last_res_seq = old_res
            
            glycan_residue_mapping[(old_chain, old_res)] = (new_glycan_chain, current_res)
    
    # Process each line of the file
    output_lines = []
    conect_lines = []
    atom_counter = 1
    
    for line in lines:
        data = parse_pdb_line(line)
        
        if not data:
            output_lines.append(line.rstrip('\n'))
            continue
        
        if data["record"] in ["ATOM", "HETATM"]:
            old_chain = data["chain"]
            old_res = data["resSeq"]
            
            # Apply mapping
            if (old_chain, old_res) in protein_residue_mapping:
                new_chain, new_res = protein_residue_mapping[(old_chain, old_res)]
                data["chain"] = new_chain
                data["resSeq"] = new_res
            elif (old_chain, old_res) in glycan_residue_mapping:
                new_chain, new_res = glycan_residue_mapping[(old_chain, old_res)]
                data["chain"] = new_chain
                data["resSeq"] = new_res
            
            # Renumber atoms sequentially
            data["serial"] = atom_counter
            atom_counter += 1
            
            output_lines.append(format_pdb_line(data))
            
        elif data["record"] == "CONECT":
            # Save CONECT lines for later processing
            conect_lines.append(data["line"])
        else:
            output_lines.append(data["line"])
    
    # Process CONECT lines
    for conect_line in conect_lines:
        output_lines.append(conect_line)
    
    # Write the output file
    with open(output_file, 'w') as f:
        for line in output_lines:
            f.write(line + '\n')
    
    print(f"\nSuccessfully processed {len(lines)} lines.")
    print(f"Output written to: {output_file}")
    
    # Summary report
    print("\n" + "=" * 60)
    print("SUMMARY")
    print("=" * 60)
    
    # Count residues per chain in the output file
    output_chains = {}
    for line in output_lines:
        if len(line) >= 26 and line[0:6].strip() in ["ATOM", "HETATM"]:
            chain = line[21:22].strip() or " "
            try:
                res_seq = int(line[22:26].strip())
            except ValueError:
                continue
            if chain not in output_chains:
                output_chains[chain] = set()
            output_chains[chain].add(res_seq)
    
    protein_residues = 0
    glycan_residues = 0
    
    for chain in sorted(output_chains.keys()):
        residues = output_chains[chain]
        count = len(residues)
        if chain in protein_chain_letters.values():
            protein_residues += count
            print(f"Chain {chain} (protein): {count} residues "
                  f"(range: {min(residues)}-{max(residues)})")
        else:
            glycan_residues += count
            print(f"Chain {chain} (glycan): {count} residues "
                  f"(range: {min(residues)}-{max(residues)})")
    
    print(f"\nTotal protein residues: {protein_residues}")
    print(f"Total glycan residues: {glycan_residues}")

def main():
    """Main function."""
    parser = argparse.ArgumentParser(
        description='PDB Residue Renumbering and Relabeling Tool',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  python pdb_renumber.py input.pdb output.pdb
  python pdb_renumber.py input.pdb --start_number_protein 100 --start_number_carb 50
  python pdb_renumber.py input.pdb --residue_carb_initiation "NDG,A2G,MAN"

The script automatically uses:
- Chain A for first protein, Chain B for second protein, etc.
- Chain A for glycans attached to first protein, Chain B for glycans attached to second protein, etc.
- Residue numbering starts from 1 for both protein and carbohydrate by default
- Carbohydrate blocks are identified by residues starting with NDG or A2G (configurable)
        """
    )
    
    parser.add_argument('input_file', help='Input PDB file')
    parser.add_argument('output_file', nargs='?', help='Output PDB file (optional, will add _renumbered if not specified)')
    parser.add_argument('--start_number_protein', type=int, default=1,
                        help='Starting residue number for protein chains (default: 1)')
    parser.add_argument('--start_number_carb', type=int, default=1,
                        help='Starting residue number for carbohydrate chains (default: 1)')
    parser.add_argument('--residue_carb_initiation', type=str, default="NDG,A2G",
                        help='Comma-separated list of residue names that initiate a new glycan block (default: NDG,A2G)')
    
    args = parser.parse_args()
    
    input_file = args.input_file
    
    if not os.path.exists(input_file):
        print(f"Error: Input file '{input_file}' not found.")
        sys.exit(1)
    
    if args.output_file:
        output_file = args.output_file
    else:
        base, ext = os.path.splitext(input_file)
        output_file = f"{base}_renumbered{ext}"
    
    # Process the file
    process_pdb_file(
        input_file, 
        output_file,
        start_number_protein=args.start_number_protein,
        start_number_carb=args.start_number_carb,
        residue_carb_initiation=args.residue_carb_initiation
    )
    
    print("\nDone!")

if __name__ == "__main__":
    main()
