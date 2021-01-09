#!/usr/bin/env nextflow

params.query = "$baseDir/data/*.fastq"
 

 Channel
    .fromPath(params.query)
    .map { file -> tuple(file.baseName, file) }
    .splitFastq(by: 1, file:true)
    .set { fastq_files }


Channel
    .fromPath(params.query)
    .splitFastq(by:1, record: true)
    .map{ record -> record.readHeader }
    .merge(fastq_files)
    .into{fastq_ch1; fastq_ch2}

process filter_length {
    input:
    tuple read_id, id , path("single_read.fastq") from fastq_ch1
    output:
    tuple read_id, id, path("${id}.fastq") into sequences_ch
    
    """
    awk 'BEGIN {OFS = "\\n"} {header = \$0 ; getline seq ; getline qheader ; getline qseq ; if (length(seq) >= 10000) {print header, seq, qheader, qseq}}' < single_read.fastq > ${id}.fastq
    """
}


process split_sequence {
    errorStrategy 'ignore'
    input:
    tuple read_id, id, path("input.fastq") from sequences_ch
    output:
    tuple read_id, id, path("${id}.fastq") into split_ch
 
    """
    split_sequence.py input.fastq ${id}.fastq 1500
    """
}

process run_assembly {
    input:
    tuple read_id, id, path("input.fastq") from split_ch
 
    output:
    tuple read_id, path("${id}.fasta") into assembly_ch
    tuple read_id, id, path("${id}.fasta") into assembly_results_ch
 
    """
    minimap2 -x ava-ont input.fastq input.fastq > ${id}.paf
    miniasm -s 800 -f input.fastq ${id}.paf > ${id}.gfa
    awk '/^S/{print ">"\$2"\\n"\$3}' ${id}.gfa  > ${id}.fasta
    sed -i "s/>.*/&_${read_id}/" ${id}.fasta
    """
}

assembly_results_ch
   .collectFile(storeDir: "output/assemblies/") { item ->
       [ "${item[1]}.fasta", item[2]]}
   .println { "Merged assemblies results into " + it }


process run_map_back {
    input:
    tuple read_id, id, path("split.fastq"), path("plasmid_candidat.fasta") from fastq_ch2.join(assembly_ch)
 
    output:
    tuple read_id, id, path("${id}.tsv") into remap_ch
 
    """
    minimap2 -ax map-ont split.fastq plasmid_candidat.fasta -p 0.0 > ${id}.sam
    bedtools bamtobed -i ${id}.sam | sort -nk 2 > ${id}.bed
    Rscript ${baseDir}/scripts/filter_it.R ${id}.bed ${id}.tsv 50
    """
}

remap_ch
   .collectFile(keepHeader: true, skip: 1, storeDir: "output/remapping/") { item ->
       [ "${item[1]}.tsv", item[2]]}
   .println { "Merged remapping results into " + it }

