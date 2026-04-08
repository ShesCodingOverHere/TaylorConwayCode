############################################################
# FST and Delta P Calculation
# Description:
# Computes FST and allele frequency differences between populations.
############################################################

import math

# =========================
# Input files
# =========================
file_pop1 = "data/input/pop1_freqs.out"
file_pop2 = "data/input/pop2_freqs.out"

output_file = "results/output/fst_deltaP.out"

# =========================
# Read data
# =========================
def read_dict(file):
    d = {}
    
    for line in open(file):
        A = line.split()
        
        span = A[0].split('.')
        chrom, start, stop = span[0], int(span[1]), int(span[2])
        
        num, denom = map(int, A[1].split('/'))
        freq = num / denom if denom > 0 else 0
        
        d[A[0]] = [chrom, start, stop, num, denom, freq]
    
    return d

p1 = read_dict(file_pop1)
p2 = read_dict(file_pop2)

# =========================
# Compute statistics
# =========================
out = open(output_file, 'w')

for key in p1:
    
    chrom, start, stop, num1, denom1, freq1 = p1[key]
    num2, denom2, freq2 = p2[key][3:]
    
    total = denom1 + denom2
    if total == 0:
        continue
    
    p_bar = (num1 + num2) / total
    q_bar = 1 - p_bar
    
    if p_bar * q_bar == 0:
        fst = 0
    else:
        c1 = denom1 / total
        c2 = denom2 / total
        
        fst = ((p_bar * q_bar) - ((freq1*(1-freq1)*c1) + (freq2*(1-freq2)*c2))) / (p_bar * q_bar)
    
    delta_p = freq1 - freq2
    
    midpoint = math.trunc((start + stop) / 2)
    
    out.write(f"{chrom}\t{start}\t{stop}\t{midpoint}\t{delta_p}\t{fst}\n")
