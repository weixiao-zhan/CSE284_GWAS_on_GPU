import gwas

if __name__ == "__main__":
    # set the file paths
    pheno_path = '../data/example.phen'
    vcf_path = "../data/example.vcf"
    out_path = "../data/out.assoc.linear" 

    # Run GWAS
    gwas.gwas(pheno_path, vcf_path, out_path)