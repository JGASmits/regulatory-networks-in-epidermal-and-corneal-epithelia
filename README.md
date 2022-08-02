# regulatory-networks-in-epidermal-and-corneal-epithelia
code files related to the manuscript: Cell identity control coordinated by common and distinct transcription factors and associated gene regulatory networks in epidermal and corneal epithelia

For regenerating all the figures from the manuscript, follow the subsequent steps:

1. Download and pre-process all the ATAC, RNA, scRNAseq, Chipseq (narrow) and Chispeq broad data using seq2science (0.7.1). All used sample and config files are present in the respective 'data' folders. Change the path to the fastq directory and other settings where needed for your specific setup. For more information regarding running seq2science, see: https://github.com/vanheeringen-lab/seq2science

2. install conda and generate conda environments from the yaml files present in 'conda yamls/'

3. RNAseq analysis (conda env bulk_rna_seq)
Optional: regenerate the pseudobulk files from the scRNAseq-pseudobulk files and or in invivo samples using the code in: "Rmarkdown/1. Generate_scRNAseq_pseudobulk" and "Rmarkdown/2. Generate_invivo_pseudobulk.Rmd".
Alternatively to skip pre-processsing use the present counttable in data/RNA_seqfiles/merged_counttables.tsv
Follow the steps written in "Rmarkdown/3. Bulk_RNA_seq.Rmd" to perform bulkRNA seq analysis (leading to the figure 1 figures).

4. Identify Cis-Regulatory Elements  (conda env ID_CRE_env)
Optional:Follow the steps written in "notebooks/1_Identify_CREs.ipynb" to identify cis regulatory elements. Change the file locations in data/ID_CRE_files/ where needed.
Alternatively skip this step and use data/Motif_enrichment_files/all_CREs_data.csv as a result for the motif enrichment instead.

5. Perform motif enrichement analysis (conda env ID_CRE_env)
Follow the steps written in "notebooks/2_Motif_enrichment.ipynb" to perform motif enrichment analysis. Change the file locations in data/ID_CRE_files/Motif_enrichment_files where needed.
After running motif enrichment for enrichment of genes close to these peaks activate the conda environment bulk_rna_seq and follow the steps in: 
"Rmarkdowns/4. CREs_MotifEnrichment"

