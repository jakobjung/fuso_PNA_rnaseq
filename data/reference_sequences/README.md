```bash
# edit cds from ncbi
sed -Ei "s/.*locus_tag=([^]]+).*/>\\1/" CDS_FNN_23726.fasta
sed -Ei "s/.*locus_tag=([^]]+).*/>\\1/" CDS_FNN_25586.fasta 

# run proteinortho to find orthologues
proteinortho6.pl --project=FNN_mapping ./CDS_FNN_falk.fasta ./CDS_FNN_25586.fasta --p=blastn --e=500 --identity=5 --cov=5
```

```bash
grep -P "old_locus_tag" FNV_A2_w_old_lf.gff3 | sed -E "s/.*;locus_tag=([^;]+);.*old_locus_tag=([^\\%; ]+).*/\1\t\2/" > mappings_FNV.tsv
```



Get genes (cds) from 

