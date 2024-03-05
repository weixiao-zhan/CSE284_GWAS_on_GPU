from cyvcf2 import VCF
import numpy as np
import torch
from scipy.stats import t
from tqdm import tqdm

ACTIVE_DEVICE = torch.device("cuda" if torch.cuda.is_available() else "cpu")

def get_id_phen(pheno_path):
    """
    Reads a .phen file and converts it into a dictionary with individual ID as keys and phenotype as values.
    
    Parameters:
    - pheno_path: path to the .phen file
    
    Returns:
    - A dictionary where keys are individual IDs and values are phenotypes.
    """
    id_to_phen_dict = {}
    with open(pheno_path, 'r') as file:
        for line in file:
            parts = line.strip().split('\t')
            if len(parts) >= 3:
                individual_id, phenotype = parts[0], parts[2]
                id_to_phen_dict[individual_id] = float(phenotype)
    return id_to_phen_dict

def build_Y(vcf, id_to_phen_dict, centered = True):
    """
    Parameters:
    - vcf: cyvcf2.VCF object
    - id_to_phen_dict: A dictionary maps individual IDs to values are phenotypes
    
    Returns:
    - (d,) shape PyTorch array representing the Y
    - (d+,) mask indicate sample has phenotype
    """
    pheno = []
    samples_missing_pheno = []
    for sample in vcf.samples:
        if sample in id_to_phen_dict:
            pheno.append(id_to_phen_dict[sample])
        else:
            samples_missing_pheno.append(sample)
    if samples_missing_pheno:
        print("Existing samples without phenotype. Abort", samples_missing_pheno)
    pheno_torch = torch.tensor(pheno, dtype=torch.float32, device=ACTIVE_DEVICE)
    if centered:
        mean = torch.mean(pheno_torch)
        pheno_torch = (pheno_torch - mean)
    return pheno_torch

def nonzero_genotype_mask(genotype):
    pass

def vcf_batch_iter(vcf, batch_size=1024):
    """
    Iterates over a VCF file in batches.
    
    Parameters:
    - vcf: : cyvcf2.VCF object
    - batch_size: Number of variants per batch.
    
    Yields:
    - (batch_size, # of sample, 3) torch tensor.
    """

    batch_CHR = []
    batch_SNP = []
    batch_BP = []
    batch_A1 = []
    batch_genotypes = []

    for snp in vcf:
        
        batch_CHR.append(snp.CHROM)
        batch_SNP.append(snp.ID)
        batch_BP.append(snp.POS)
        batch_A1.append(snp.ALT[0])
        
        batch_genotypes.append(snp.genotypes)
        if len(batch_genotypes) == batch_size:
            np_CHR = np.array(batch_CHR)
            np_SNP = np.array(batch_SNP)
            np_BP = np.array(batch_BP)
            np_A1 = np.array(batch_A1)
            np_genotypes = np.sum(
                np.array(batch_genotypes, dtype=np.uint8)[:, :, :2], 
                axis=2)
            np_genotypes = np_genotypes - np.mean(np_genotypes, axis=1, keepdims=True)
            valid_mask = np.any(np_genotypes != 0, axis=1)
            yield np_CHR[valid_mask], np_SNP[valid_mask], np_BP[valid_mask], np_A1[valid_mask], np_genotypes[valid_mask,:]
            # Reset batch list after yielding
            batch_CHR = []
            batch_SNP = []
            batch_BP = []
            batch_A1 = []
            batch_genotypes = []
    
    # Yield any remaining variants in the last batch (if not empty)
    if batch_genotypes:
            np_CHR = np.array(batch_CHR)
            np_SNP = np.array(batch_SNP)
            np_BP = np.array(batch_BP)
            np_A1 = np.array(batch_A1)
            np_genotypes = np.sum(
                np.array(batch_genotypes, dtype=np.uint8)[:, :, :2], 
                axis=2)
            np_genotypes = np_genotypes - np.mean(np_genotypes, axis=1, keepdims=True)
            valid_mask = np.all(np_genotypes == 0, axis=1)
            yield np_CHR[valid_mask], np_SNP[valid_mask], np_BP[valid_mask], np_A1[valid_mask], np_genotypes[valid_mask]

def init_output(out_path):
    outfile = open(out_path, 'w')
    header = "CHR\tSNP\tBP\tA1\tBETA\tSTAT\tP\n"
    outfile.write(header)
    return outfile

def batch_output(outfile, batch_CHR, batch_SNP, batch_BP, batch_A1, batch_Beta, batch_STAT, batch_P):
    for i in range(len(batch_CHR)):
        outfile.write(f"{batch_CHR[i]}\t{batch_SNP[i]}\t{batch_BP[i]}\t{batch_A1[i]}\t{batch_Beta[i]:.4f}\t{batch_STAT[i]}\t{batch_P[i]}\n")

def gwas(pheno_pth, vcf_path, outfile):
    id_phen_dict = get_id_phen(pheno_pth)
    vcf = VCF(vcf_path)
    num_samples = len(vcf.samples)
    df = num_samples - 2

    y = build_Y(vcf, id_phen_dict, centered=True)

    outfile = init_output(outfile)

    bar = tqdm(vcf_batch_iter(vcf, batch_size=1024))
    for batch_CHR, batch_SNP, batch_BP, batch_A1, batch_genotypes in bar:
        X = torch.tensor(batch_genotypes, dtype=torch.float32, device=ACTIVE_DEVICE)
        Xy = torch.sum(X * y, dim=1, keepdim=True)
        XX = torch.sum(X * X, dim=1, keepdim=True)
        beta = (Xy / XX)
        stat = beta / torch.sqrt(torch.var(y -beta*X) / XX)

        beta_cpu = beta.cpu().detach().numpy().squeeze()
        stat_cpu = stat.cpu().detach().numpy().squeeze()
        p = 2 * t.sf(np.abs(stat_cpu), df)
        batch_output(outfile, batch_CHR, batch_SNP, batch_BP, batch_A1, beta_cpu, stat_cpu, p)
    
    outfile.close()

if __name__ == "__main__":
    pass