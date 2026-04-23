# CBTi-seq_tools

The instruction of CBTi-seq_tools

The bioinformatics processing scripts for CBTi-seq, from raw sequencing data to gene expression matrixs.

For more information on CBTi-seq, please refer to the paper "He LY, et al. High-Resolution Multiplexed Sequencing of Single-Cell Full-length Transcriptome Via Combinational Barcoded Tn5 Transposon Insertion, Advanced Sciences, 2026, e16013, doi: 10.1002/advs.202516013."

ⅠDocuments and tools needed
  Picard – 2.13.3, 
  dropseq_tools – 3.0.2, 
  STAR – 2.7.11b, 
  genome.fasta -- GRCh38.fa, 
  genome.refFlat -- GRCh38.99.refFlat, 
  python – 3.12, 
  setuptools – 75.1.0,
  pyyaml – 6.0.2, 
  pysam – 0.22.1,
  pandas – 2.2.3,
  java – java23

ⅡSetting and setup
  i. Open CBTi-seq_tools/config/paths.yaml and replace the paths to the necessary files and tools
  ii. Run the following command in the directory: pip install -e .
  iii. If you see the following message, the installation was successful:

    Successfully built CBTi-seq_tools
    Installing collected packages: CBTi-seq_tools
    Successfully installed CBTi-seq_tools-1.0.0

ⅢHow to Run
  CBTi-seq_tools \
      --workdir <Work Directory> \
      --sample <Sample Name> \
      --r1 <The path of Read1.fastq> \
      --r2 <The path of Read2.fastq> \
      --bcA <The position of BarcodeA in Read1> \
      --bcB <The position of BarcodeB in Read2> \
      --umiA <The position of UMIA in Read1> \
      --umiB <The position of UMIB in Read2> \
      --whitelist <The path of whitelist.txt>

Note: workdir, sample, r1 and r2 are required parameters, the rest are optional.
