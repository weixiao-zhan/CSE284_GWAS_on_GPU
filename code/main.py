import gwas
import torch
import time

def test_all():
    pheno_path = '../data/ps3_gwas.phen'
    vcf_path = "../data/ps3_gwas.vcf.gz"
    # out_path = "../data/testout"
    
    batch_size_powers = [2, 6, 10, 14]
    if_use_cpu = [True, False]
    if_compute_p = [True, False]

    for batch_size_power in batch_size_powers:
        for use_cpu in if_use_cpu:
            for compute_p in if_compute_p:
                start = time.time()
                if use_cpu:
                    gwas.ACTIVE_DEVICE = torch.device("cpu")
                gwas.gwas(pheno_path, vcf_path, "", 
                          batch_size=2**batch_size_power, compute_pval=compute_p)
                end = time.time()
                print(f"batch size 2**{batch_size_power}, using {'CPU' if use_cpu else 'CUDA'}, {'computing' if compute_p else 'skipping'} p-value:\n\t{round(end - start, 4)}")

def test_one():
    pheno_path = '../data/ps3_gwas.phen'
    vcf_path = "../data/ps3_gwas.vcf.gz"
    out_path = "../data/testout"
    gwas.ACTIVE_DEVICE = torch.device("cpu")
    start = time.time()
    gwas.gwas(pheno_path, vcf_path, out_path, batch_size=2**16)
    end = time.time()
    print(round(end - start, 4))

if __name__ == "__main__":
    test_one()