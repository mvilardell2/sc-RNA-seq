# Cell Ranger

Cell Ranger is a software developed by 10x Genomics for processing and analyzing data from single-cell RNA sequencing (scRNA-seq) experiments.

Once the FASTQ files are generated by the sequencer, Cell Ranger aligns the sequencing reads to a reference genome. 


## PIPELINE TO PERFORM CELL RANGER

There are two scripts that allow to run cell ranger, both perform the execution of cell ranger for multiple samples (For count and vdj modules. For multi modules run the script for each sample). 

- Script run_cellranger.slm: Samples are processed in a loop (not array task).
  
- Script run_cellRanger_Array.slm: Samples are processed in a array task (so it has to be changed manually). 
  In the job description of the script, modify the SBATCH --array parameter depending on the number of samples to be processed. 

There is a sample.config file that has to be filled with the requiered variables. 

Create a logs folder to save the logs, and optionally a results folder. 



MODULES: 

- count: for gene expression. Allows multiple samples 
- vdj: VDJ recombination in T/B cells. Identifies clonotypes and sequences of VDJ regions. Allows multiple samples
- multi: gene expression + VDJ + antibody capture. Just for one sample each sbatch job 


For count and vdj modules, it is important that all the fastqs are located in the same directory (even dealing with different samples, then cellranger will select the fastqs for the same sample)


For the multi mode, fill the multi_config.csv file, and add the path of this file to the sample.config file. Do not need to add the FASTQS path to the sample.config. In this case, separate the FASTQS in different folders depending on the feature type (ex: gene expression fastq files in a folder, vdj fastq in another)




## Functionality: 

In single-cell experiments, each cell is tagged with a unique barcode. Cell Ranger identifies these barcodes to assign each read to its corresponding cell and it also counts the UMIs associated with each gene in each cell, providing a quantitative 
measure of gene expression --> Cell ranger outputs a filtered gene expression matrix (h5), used for downstream analysis (Seurat). 


Cell Ranger also includes a specialized module called Cell Ranger V(D)J that is designed for processing and analyzing single-cell immune profiling data, specifically focused on the adaptive immune repertoire,
including T-cell and B-cell receptors. In this particular case, instead of aligning reads to the whole transcriptome, Cell Ranger V(D)J aligns them to a specialized reference that includes sequences for V (variable), D (diversity), J (joining),
and C (constant) gene segments of immune receptors. ---> Targeted algiment to VDJ sequences.

The outputs of cell ranger vdj and multi include: annotated V(D)J contigs, clonotype information, and summaries of V(D)J usage. The important ouput files are: Clonotypes.csv, filtered_contig_annotations.csv, and airr_rearangment.tsv.

If combined with gene expression profiling, Cell Ranger V(D)J can also provide information on the gene expression profile of the same cells.

