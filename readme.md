# Project: Accelerating GWAS with GPU Computing
Project Option: Option #1 - Implementation of GWAS and Application to Real Data
Team Members: Weixiao Zhan (A59023453)

### Description of Methods: 
In PS3, I observed significant computational delays when executing plink –linear for GWAS analysis, primarily due to the absence of multi-threading support in plink 1.9. Our project intends to investigate the capacity of GPU computing to markedly expedite GWAS computations. This will be achieved by utilizing Python in conjunction with libraries that are optimized for GPU operations. We will make use of cyvcf2 for efficient VCF file parsing, PyTorch for in-memory data management, and PyTorch.GPU modules to leverage GPU acceleration capabilities. 
I will also compare the runtime performance across three different platforms: plink 1.9 (lacking multi-threading support), plink 2.0 (which includes multi-threading capabilities), and our own GPU-accelerated implementation. 

### Description of the Dataset: 
For this project, we will utilize the dataset provided in PS3, specifically ps3_gwas.phen and ps3_gwas.vcf.gz.

# Example Usage
[code/main.py](code/main.py])

# results
[report.ipynb](report.ipynb)