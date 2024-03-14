# Project: Accelerating GWAS with Parallel Computing
Project Option: Option #1 - Implementation of GWAS and Application to Real Data
Team Members: Weixiao Zhan (A59023453)

### Description of Methods: 
In PS3, I observed significant computational delays when executing plink –linear for GWAS analysis on DataHub. My project intends to investigate the capacity of GPU computing to expedite GWAS computations. This will be achieved by utilizing Python in conjunction with libraries that are optimized for GPU operations. I used cyvcf2 for efficient VCF file parsing, PyTorch for in-memory data management, and PyTorch.GPU modules to leverage GPU acceleration capabilities. 
I will also compare the runtime performance across three different tools: 
* `PLINK v1.90b6.9 64-bit (4 Mar 2019)` (the default version on DataHub), 
* `PLINK v1.90b7.2 64-bit (11 Dec 2023)`(the plink 1.9 latest stable version),
* and my implementation:
    * CPU only
    * with GPU-acceleration

### Description of the Dataset: 
The performance is evaluated on the dataset provided in PS3, specifically `ps3_gwas.phen` and `ps3_gwas.vcf.gz`.
Meanwhile, [`data/example.phen`](data/example.phen) and [`data/example.vcf`](data/example.vcf) are extracted from the PS3 data set for demonstration.

# Get started
### clone this repo
```
git clone https://github.com/weixiao-zhan/CSE284_GWAS_on_GPU.git --depth=1
cd CSE284_GWAS_on_GPU
```

### install PyTorch dependency
CPU only:
```
pip3 install torch --index-url https://download.pytorch.org/whl/cpu
```
CUDA supported *(SKIP if just using CPU version)*:
```
conda install pytorch pytorch-cuda=12.1 -c pytorch -c nvidia
```
### run example code
```
cd code
python3 main.py
```
### run plink
```
plink --linear \
--vcf data/example.vcf --pheno data/example.phen --allow-no-sex --maf 0.05 --out data/plink

```

# Report
[report.pdf](report/report.pdf)