#!/usr/bin/env nextflow

params.query = "$baseDir/data/*.fastq"
 

 Channel
    .fromPath(params.query)
    .splitFastq(by: 1, file:true)
    .set { fastq_ch }


process filter_length {
    input:
    path "single_read.fastq" from fastq_ch
 
    output:
    tuple path("single_read.fastq"), file("output.fastq") into sequences_ch
 
    """
    cat single_read.fastq | paste - - - - | awk 'length(\$2)  >= 10000 ' | sed 's/\\t/\\n/g' > output.fastq
    """
}


process split_sequence {
    errorStrategy 'ignore'
    input:
    tuple path("single_read.fastq"), file("input.fastq") from sequences_ch
    output:
    tuple path("single_read.fastq"), file("output.fastq") into split_ch
 
    """
    python $baseDir/apps/split_sequence.py input.fastq output.fastq 1500
    """
}


process run_assembly {
    publishDir "output/assemblies/", pattern: "*.fasta", mode: 'copy', overwrite: true
    input:
    tuple path("single_read.fastq"), file("input.fastq") from split_ch
 
    output:
    tuple path("single_read.fastq"), file("plasmid_candidat.fasta") into assembly_ch
 
    """
    $baseDir/apps/minimap2 -x ava-ont input.fastq input.fastq > reads.paf
    $baseDir/apps/miniasm -s 800 -f input.fastq reads.paf > reads.gfa
    awk '/^S/{print ">"\$2"\\n"\$3}' reads.gfa  > plasmid_candidat.fasta
    """
}


process run_map_back {
    publishDir  "output/remap" , mode: 'copy', overwrite: true
    input:
    tuple path("single_read.fastq"), file("plasmid_candidat.fasta") from assembly_ch
 
    output:
    file("remap.tsv") into remap_ch
 
    """
    $baseDir/apps/minimap2 -ax map-ont single_read.fastq plasmid_candidat.fasta -p 0.0 > remap.sam
    bedtools bamtobed -i remap.sam | sort -nk 2 > remap.bed
    Rscript $baseDir/apps/filter_it.R remap.bed remap.tsv 50
    """
}


