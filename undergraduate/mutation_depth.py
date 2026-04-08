############################################################
# Mutation Depth Calculator
# Description:
# Computes read depth over genomic regions using samtools.
############################################################

import subprocess

def get_depth(chrom, start, stop, bam):
    region = f"{chrom}:{start}-{stop}"
    
    cmd = ["samtools", "depth", "-r", region, bam]
    process = subprocess.Popen(cmd, stdout=subprocess.PIPE)
    
    total_depth = 0
    
    for line in process.stdout:
        depth = int(line.split()[2])
        total_depth += depth
    
    length = stop - start + 1
    
    return total_depth / length

# =========================
# Input
# =========================
regions_file = "data/input/regions.out"

for line in open(regions_file):
    
    A = line.split()
    
    chrom = A[4]
    start = int(A[5])
    stop = int(A[3])
    
    strains = A[7].split(',')
    
    for strain in strains:
        bam = f"data/input/bams/{strain}.bam"
        
        avg_depth = get_depth(chrom, start, stop, bam)
        
        print(strain, chrom, start, stop, avg_depth)
