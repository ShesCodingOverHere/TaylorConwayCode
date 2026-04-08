############################################################
# Allele Frequency Calculator
# Description:
# Computes allele frequencies from haplotype output files.
############################################################

# =========================
# Input file (EDIT PATH)
# =========================
input_file = "data/input/haplotypes.out"

# =========================
# Read file
# =========================
InFile = open(input_file, 'r')

# Dictionary: key = chrom.start.stop
MyDict = {}

# =========================
# Parse input
# =========================
for line in InFile:
    A = line.split()
    
    key = A[1] + '.' + A[2] + '.' + A[3]
    
    value = (A[-2], A[-1], A[0])
    
    if key in MyDict:
        MyDict[key].append(value)
    else:
        MyDict[key] = [value]

# =========================
# Compute allele frequencies
# =========================
for key in MyDict:
    
    num = []
    denom = []
    
    for entry in MyDict[key]:
        
        n = int(entry[0])
        d = int(entry[1])
        
        # Fix impossible case (2/1)
        if n == 2 and d == 1:
            num.append(1)
            denom.append(1)
        else:
            num.append(n)
            denom.append(d)
    
    numerator = sum(num)
    denominator = sum(denom)
    
    # Skip zero-frequency sites (optional)
    if numerator == 0:
        continue
    
    print(key, f"{numerator}/{denominator}")
