# Cancer Microbiome Scripts
The Cancer Microbiome Scripts repository is a collection of miscellaneous scripts used during analyzing the microbiome in tumour material based on metatranscriptomic data. 

Created during the BSc graduation internship of Birgit Rijvers at the Clinical Bioinformatics group of Erasmus University Medical Center.

The main goal of this internship was to create a pipeline that: 
- Takes sequencing data (FASTQ) as input
- Performs quality control
- Maps reads to human reference genome (host depletion)
- Removes reads that mapped
- Taxonomically classifies the reads left 
- Reports quality metrics and reads counts at the main steps

The final version of this pipeline is available in another repository: https://github.com/BirgitRijvers/emc-cancermicro, but this repository contains the first versions of the pipeline which were written in Python.

## Scripts included
### Pipeline scripts
#### switchpipe.py
Runs the full pipeline, switches between commands based on single or paired end input data.
#### pipeline.py
Runs pipeline without QC on paired end data.
#### pipeline_single.py
Runs pipeline without QC on single end data.
### Other scripts
#### ncbi_id_fetcher.py
Fetches taxonomy IDs from NCBI based on taxonomic names.
#### kraken_classification_counter.py
Counts certain classifications (Human, root, unclassified & user specified) in Kraken2 output (not report style).
#### bracken_runner.py
Runs Bracken on all Kraken2 reports in a directory.
#### samplesheeter.py
Creates a CSV samplesheet compatible with nf-core pipelines based on (compressed) FASTQ files in a directory.

## Citations
### Tools
[Fastp](https://github.com/OpenGene/fastp)
> Chen S, Zhou Y, Chen Y, Gu J. fastp: an ultra-fast all-in-one FASTQ preprocessor. Bioinformatics. 2018 Sep 1;34(17):i884-i890. doi: 10.1093/bioinformatics/bty560. PMID: 30423086; PMCID: PMC6129281.

[Kraken2](https://github.com/DerrickWood/kraken2) 
> Wood, D.E., Lu, J. & Langmead, B. Improved metagenomic analysis with Kraken 2. Genome Biol 20, 257 (2019). https://doi.org/10.1186/s13059-019-1891-0

[BWA-MEM2](https://github.com/bwa-mem2/bwa-mem2)
>Vasimuddin Md, Sanchit Misra, Heng Li, Srinivas Aluru. Efficient Architecture-Aware Acceleration of BWA-MEM for Multicore Systems. IEEE Parallel and Distributed Processing Symposium (IPDPS), 2019. 10.1109/IPDPS.2019.00041

[Samtools](https://www.htslib.org/doc/samtools.html)
>Twelve years of SAMtools and BCFtools
Petr Danecek, James K Bonfield, Jennifer Liddle, John Marshall, Valeriu Ohan, Martin O Pollard, Andrew Whitwham, Thomas Keane, Shane A McCarthy, Robert M Davies, Heng Li
GigaScience, Volume 10, Issue 2, February 2021, giab008, https://doi.org/10.1093/gigascience/giab008

[kraken-biom](https://github.com/smdabdoub/kraken-biom)
> Dabdoub, SM (2016). kraken-biom: Enabling interoperative format conversion for Kraken results (Version 1.2) [Software]. Available at https://github.com/smdabdoub/kraken-biom.

[eUtils](https://www.ncbi.nlm.nih.gov/books/NBK25501/)
>Entrez Programming Utilities Help [Internet]. Bethesda (MD): National Center for Biotechnology Information (US); 2010-. Available from: https://www.ncbi.nlm.nih.gov/books/NBK25501/
