# Project: Accelerating GWAS with Parallel Computing
Project Option: Option #1 - Implementation of GWAS and Application to Real Data
Team Members: Weixiao Zhan (A59023453)

### Report
[report.pdf](report/report.pdf)

# Get started
### clone this repo
```
git clone https://github.com/weixiao-zhan/CSE284_ParallelGWAS.git --depth=1
cd CSE284_ParallelGWAS
```

### install dependency

1. PyTorch
    * CPU only: 
    `pip3 install torch --index-url https://download.pytorch.org/whl/cpu`
    
    * CUDA supported *(SKIP if just using CPU version)*:
    `conda install pytorch pytorch-cuda=12.1 -c pytorch -c nvidia`

2. others
`pip install numpy tqdm scipy cyvcf2`
Note: cyvcf2 would run into issue on Mac (with Arm cpu), recommend to to use datahub.

### run example code
```
cd code
python3 main.py
```

### run plink
```
plink --linear 
    --vcf data/example.vcf 
    --pheno data/example.phen 
    --allow-no-sex --maf 0.05 
    --out data/plink
```