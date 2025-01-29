This directory contains a script (`script_to_extract_16s.sh`) that:
1. Detects 16S rRNA annotations from each genome’s GFF file.
2. Extracts the corresponding 16S rRNA sequences from the matching FASTA.
3. Outputs each strain’s 16S rRNA sequences in a separate FASTA file.

## Requirements

- **bedtools** (for `bedtools getfasta`)
- **awk** and **grep** (standard POSIX tools)
- **bash** 4.0 or higher
- **mafft** (for multiple sequence alignment)
- **fasttree** (for phylogenetic tree construction)
- **barrnap** (for 16S rRNA annotation)

## Usage

1. Place all your `.fasta` and `.gff` files in the `genomes/` folder (or adjust paths in the script as needed).
2. Make the script executable:
   ```bash
   chmod +x genomes/script_to_extract_16s.sh