# plasmidPublication
A Peak into the Plasmidome of Global Sewage

# plasmidPublication

This repository contains the nextflow scripts for our publication 'A Peak into the Plasmidome of Global Sewage'

Abstract of publication:

Advances in genomics have the potential to revolutionize clinical diagnostics. Here, we examine the microbiome of vitreous (intraocular body fluid) from patients who developed endophthalmitis following cataract surgery or intravitreal injection. Endophthalmitis is an inflammation of the intraocular cavity and can lead to a permanent loss of vision. As controls, we included vitreous from endophthalmitis-negative patients, balanced salt solution used during vitrectomy and DNA extraction blanks. We compared two DNA isolation procedures and found that an ultraclean production of reagents appeared to reduce background DNA in these low microbial biomass samples. We created a curated microbial genome database (>5700 genomes) and designed a metagenomics workflow with filtering steps to reduce DNA sequences originating from: (i) human hosts, (ii) ambiguousness/contaminants in public microbial reference genomes and (iii) the environment. Our metagenomic read classification revealed in nearly all cases the same microorganism than was determined in cultivation- and mass spectrometry-based analyses. For some patients, we identified the sequence type of the microorganism and antibiotic resistance genes through analyses of whole genome sequence (WGS) assemblies of isolates and metagenomic assemblies. Together, we conclude that genomics-based analyses of human ocular body fluid specimens can provide actionable information relevant to infectious disease management.

Link to publication: 
TBD

Cite as: 
TBD


## Getting Started

These instructions will let you run the assembly workflow.

###Dependencies

Java 8 <br/>
Nextflow <br/>
Docker <br/>

### Running the workflow

1. Create docker image

  ```
  docker build -t plasmid_publication <dir>/docker
  ```

2. Run workflow
 
  Move fastq-file into the data folder.
 
  ```
  nextflow run main.nf -with-docker plasmid_publication
  ```
  
3. Output

  Folder with candidate assemblies per sample (\<dir\>/output/assemblies) <br/>
  Folder with stats for each candidate per sample (\<dir\>/output/remapping) <br/>

## Author of scripts

* **Philipp Kirstahler** - [philDTU](https://github.com/philDTU)
