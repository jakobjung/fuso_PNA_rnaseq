#!/usr/bin/env bash

# Requirements:
#   - bedtools (for getfasta)
#   - barrnap (for rRNA detection)
#   - awk / grep
touch final_fuso_16s.fasta

for fasta in *.fasta; do
    # Remove extension to get prefix
    prefix=$(basename "$fasta" .fasta)
    gff="${prefix}.gff"
    
    # We'll store final rRNA annotation here:
    rrna_gff="${prefix}_rRNA.gff"

    # 1. If the GFF doesn't exist, or doesn't have 16S lines, use Barrnap:
    if [[ ! -f "$gff" ]]; then
        echo "No GFF found for ${prefix}. Running Barrnap..."
        barrnap "$fasta" > "$rrna_gff"
    else
        # Check if there's any "16S" mention in the existing GFF
        has_16s=$(grep -i "16S" "$gff" | wc -l)
        if (( has_16s == 0 )); then
            echo "No 16S annotation in ${gff}. Running Barrnap..."
            barrnap "$fasta" > "$rrna_gff"
        else
            echo "Found 16S annotation in ${gff}."
            cp "$gff" "$rrna_gff"
        fi
    fi

    # 2. Extract 16S features from the chosen GFF (either original or Barrnap-generated).
    bed="${prefix}_16S.bed"
    out="${prefix}_16S.fasta"

    # Grep for "16S" in the chosen rRNA GFF
    grep -i "16S" "$rrna_gff" | \
    awk 'BEGIN{FS="\t"; OFS="\t"} {print $1,$4,$5,$9,".",$7}' > "$bed"

    # 3. Run bedtools to get the sequence
    bedtools getfasta -fi "$fasta" -bed "$bed" -fo "$out" -s -name

    # for all 16s fastas, only preserve first 2 lines and save them as final_fuso_16s.fasta
    head -n 2 "$out" >> final_fuso_16s.fasta

    # Clean up if you want:
    #rm -f "$bed"
    # remove the single 16s fasta file
    rm -f "$out"
    # remove out .fai
    rm -f "$out.fai"
    # remove the rRNA gff file
    rm -f "$rrna_gff"

    echo "Extracted 16S rRNA(s) for $prefix into $out."
    
    # Optional: keep or remove the Barrnap GFF if you only needed it temporarily
    # rm -f "$rrna_gff"
done

# use mafft to align the 16s sequences
mafft --auto final_fuso_16s.fasta > final_fuso_16s_aligned.fasta

# ok now use fasttree to build a tree
FastTree final_fuso_16s_aligned.fasta > final_fuso_16s_aligned.tree


