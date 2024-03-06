def parse_gt_line(line):
    fields = line.split()
    return {
        'CHR': fields[0],
        'SNP': fields[1],
        'BP': int(fields[2]),
        'A1': fields[3],
        'BETA': float(fields[6]),
        'STAT': float(fields[7]),
        'P': float(fields[8]),
    }

def parse_output_line(line):
    fields = line.split()
    return {
        'CHR': fields[0],
        'SNP': fields[1],
        'BP': int(fields[2]),
        'A1': fields[3],
        'BETA': float(fields[4]),
        'STAT': float(fields[5]),
        'P': float(fields[6]),
    }

def check_line(gt_line, out_line, delta = 0.1):
    """Check the fields of a line from the ground truth and output files."""
    gt_fields = parse_gt_line(gt_line)
    out_fields = parse_output_line(out_line)
    
    for key in ['CHR', 'SNP']:  # Check string fields for exact match
        if gt_fields[key] != out_fields[key]:
            return False
    
    # for key in ['BP', 'BETA', 'STAT', 'P']:  # Check numerical fields with a threshold
    #     if key == 'BP':  # 'BP' is an integer, so check for exact equality
    #         if gt_fields[key] != out_fields[key]:
    #             return False
    #     else:  # For 'BETA', 'STAT', and 'P', allow a small error
    #         if abs(gt_fields[key] - out_fields[key]) > delta:
    #             return False
    
    return True

def check_output(out_path, gt_path):
    with open(out_path, 'r') as out_file, open(gt_path, 'r') as gt_file:
        # skip the title
        next(out_file, None)
        next(gt_file, None)

        for out_line, gt_line in zip(out_file, gt_file):
            if not check_line(out_line.strip(), gt_line.strip()):
                print(f"Mismatch found:\nOutput: \n{out_line}\nGround Truth: \n{gt_line}")
                return False
            
        for out_line in out_file:
            print(f"Extra line in output file: {out_line}")
            return False
        for gt_line in gt_file:
            print(f"Missing line in output, present in ground truth: {gt_line}")
            return False
    return True
# Example usage
# outpath = "path/to/your/output/file"
# ground_truth = "path/to/your/ground/truth/file"
# check_output(outpath, ground_truth)
if __name__ == '__main__':
    out_path = "../data/testout"
    gt_path = "../data/ps3_gwas.assoc.linear"
    if check_output(gt_path, out_path):
        print('passed')