import gwas
import torch
import time

if __name__ == "__main__":
    pheno_pth = '../data/ps3_gwas.phen'
    vcf_path = "../data/ps3_gwas.vcf.gz"
    outfile = "../data/testout"
    
    start = time.time()
    gwas.ACTIVE_DEVICE = torch.device("cpu")
    gwas.gwas(pheno_pth, vcf_path, outfile)
    
    end = time.time()
    print(end - start)