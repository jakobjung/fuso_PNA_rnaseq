#!/bin/bash
main(){
    # go to project directory where the reads and reference
    # sequences are stored:
    PROJECT=../data
    echo "Start trimming"
    #rename_trim_rna_libs
    echo "Trimming done. Start mapping"
    #align_rna_reads_genome
    #align_rna_reads_genome_2
    
    echo "Finished mapping. Start connecting all tab files"
    featureCounts -T 5 -t CDS,ncRNA -g id \
		  -a $PROJECT/reference_sequences/ATCC23726_updated2023.gff \
		  -o $PROJECT/rna_align/counttable_fnn23.txt \
		  $PROJECT/rna_align/*Fnn23*.bam
    featureCounts -T 5 -t rRNA -g ID \
		  -a $PROJECT/reference_sequences/ATCC23726_updated2023.gff\
		  -o $PROJECT/rna_align/counttable_fnn23_rRNAs.txt \
		  $PROJECT/rna_align/*Fnn23*.bam
     featureCounts -T 5 -t tRNA -g ID \
		  -a $PROJECT/reference_sequences/ATCC23726_updated2023.gff\
		  -o $PROJECT/rna_align/counttable_fnn23_tRNAs.txt \
		  $PROJECT/rna_align/*Fnn23*.bam

     
    featureCounts -T 5 -t CDS,ncRNA -g id \
		  -a $PROJECT/reference_sequences/FNV3_1_36A.gff \
		  -o $PROJECT/rna_align/counttable_fnv.txt \
		  $PROJECT/rna_align/*Fnv*.bam
    featureCounts -T 5 -t rRNA -g locus_tag \
		  -a $PROJECT/reference_sequences/FNV3_1_36A.gff\
		  -o $PROJECT/rna_align/counttable_fnv_rRNAs.txt \
		  $PROJECT/rna_align/*Fnv*.bam
     featureCounts -T 5 -t tRNA -g locus_tag \
		  -a $PROJECT/reference_sequences/FNV3_1_36A.gff\
		  -o $PROJECT/rna_align/counttable_fnv_tRNAs.txt \
		  $PROJECT/rna_align/*Fnv*.bam     
    
}

rename_trim_rna_libs(){
    mkdir -pv $PROJECT/libs
    # shellcheck disable=SC2045
    for NAME in $(ls $PROJECT/fastq/*.fq.gz)
    do
        echo "$NAME starts trimming nowwwwwww"
        NEWNAME=${NAME##*/}
        NEWNAME=${NEWNAME%.fastq.gz}_trimmed.fastq.gz
        echo $NEWNAME
       # bbduk trims low quality bases and removes adapters:
        bbduk.sh  in=$NAME \
			     ref=../data/reference_sequences/adapters.fa -Xmx4g t=20\
			     out=$PROJECT/libs/${NEWNAME} ktrim=r k=23 mink=11\
			     hdist=1 qtrim=r trimq=10 ftl=12 
    done
}


align_rna_reads_genome(){
    mkdir -p $PROJECT/rna_align
    DIR=$PROJECT/rna_align
    for i in $(ls $PROJECT/libs/*Fnn23*.fastq.gz)
    do
        NAME=${i##*/}
        NAME=${NAME%_trimmed.fastq.gz}
        echo "Starting mapping for sample: $NAME"
	echo "$i"
        bbmap.sh in=$i trimreaddescription=t  t=20 ref=$PROJECT/reference_sequences/ATCC23726_ncbi.fasta k=12 outm=$DIR/$NAME.sam 
        # sort sam file, create BAM file:
        samtools sort -O BAM -@ 40 $DIR/$NAME.sam > $DIR/$NAME.bam
        # remove sam file: (not actually needed)
        # rm $DIR/$NAME.sam
        # index BAM file:
        samtools index $DIR/$NAME.bam
        # generate coverage statistics with bamcoverage:
        bamCoverage -b $DIR/$NAME.bam -o $DIR/$NAME_forward.bw -bs 1 --filterRNAstrand forward
	bamCoverage -b $DIR/$NAME.bam -o $DIR/$NAME_reverse.bw -bs 1 --filterRNAstrand reverse
	fastqc $i
	#rm $i
    done
}

align_rna_reads_genome_2(){
    mkdir -p $PROJECT/rna_align
    DIR=$PROJECT/rna_align
    for i in $(ls $PROJECT/libs/*Fnv*.fastq.gz)
    do
        NAME=${i##*/}
        NAME=${NAME%_trimmed.fastq.gz}
        echo "Starting mapping for sample: $NAME"
        bbmap.sh in=$i trimreaddescription=t  t=20 \
			     ref=$PROJECT/reference_sequences/FNV3_1_36A.fasta \
			     k=12 outm=$DIR/$NAME.sam 
        # sort sam file, create BAM file:
        samtools sort -O BAM -@ 40 $DIR/$NAME.sam > $DIR/$NAME.bam
        # remove sam file: (not actually needed)
        rm $DIR/$NAME.sam
        # index BAM file:
        samtools index $DIR/$NAME.bam
        # generate coverage statistics with bamcoverage:
        bamCoverage -b $DIR/$NAME.bam -o $DIR/$NAME_forward.bw -bs 5 --filterRNAstrand forward
	bamCoverage -b $DIR/$NAME.bam -o $DIR/$NAME_reverse.bw -bs 5 --filterRNAstrand reverse
	
	fastqc $i
	#rm $i
    done
}

main  
