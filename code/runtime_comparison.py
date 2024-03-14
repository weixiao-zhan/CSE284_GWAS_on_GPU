import gwas
import torch
import time

if __name__ == "__main__":
    pheno_path = "../data/ps3_gwas.phen"
    vcf_path = "../data/ps3_gwas.vcf.gz"
    out_path = ""
    
    batch_size_powers = [2, 6, 10, 14]
    if_use_cpu = [True, False]
    if_compute_p = [True, False]

    for batch_size_power in batch_size_powers:
        for use_cpu in if_use_cpu:
            for compute_p in if_compute_p:
                print(f"batch size 2**{batch_size_power}, using {'CPU' if use_cpu else 'CUDA'}, {'computing' if compute_p else 'skipping'} p-value:\n")
                if use_cpu:
                    gwas.ACTIVE_DEVICE = torch.device("cpu")
                gwas.gwas(pheno_path, vcf_path, "", 
                          batch_size=2**batch_size_power, compute_pval=compute_p)
