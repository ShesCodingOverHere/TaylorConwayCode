############################################################
# Site Frequency Spectrum (SFS)
############################################################

from scipy import special
import numpy as np

def choose(n, k):
    return special.binom(n, k)

max_sample = 19
SFS = [0] * (max_sample + 1)

for line in open("data/input/allele_freqs.out"):
    
    if line.startswith("Chr"):
        continue
    
    num, denom = map(int, line.split()[1].split('/'))
    
    if denom >= max_sample:
        
        for i in range(max_sample + 1):
            SFS[i] += (
                choose(num, i) *
                choose(denom - num, max_sample - i) /
                choose(denom, max_sample)
            )

SFS = SFS[1:]
total = sum(SFS)

for x in SFS:
    print(x / total)
