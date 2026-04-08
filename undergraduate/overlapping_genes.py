############################################################
# Gene Overlap Annotation
############################################################

genes = [line.split() for line in open("data/input/gene_coords.txt")]
mutations = open("data/input/mutations.out")

for line in mutations:
    
    A = line.split()
    chrom = A[1]
    mut_start = int(A[5])
    mut_stop = int(A[3])
    
    for g in genes:
        
        gene_start = int(g[4])
        gene_stop = int(g[5])
        gene_chr = g[2].split("=")[1]
        
        if gene_chr != chrom:
            continue
        
        if gene_start < mut_start and gene_stop > mut_stop:
            label = "within"
        elif mut_start > gene_start and mut_start < gene_stop:
            label = "partial_overlap"
        else:
            continue
        
        print(g[0], chrom, mut_start, mut_stop, label)
