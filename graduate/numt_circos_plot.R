############################################################
# NUMT Circos Plot
# Author: [Your Name]
# Description:
# Visualize NUMTs by linking mitochondrial genome regions
# to nuclear genome insertion sites using circos.
############################################################

# =========================
# Load libraries
# =========================
library(circlize)
library(dplyr)
library(stringr)

# =========================
# Input files (EDIT PATHS)
# =========================
mt_hits_file <- "data/input/mt_alignments.tsv"
nuclear_fasta_headers <- "data/input/nuclear_numts.fa"
mt_annotation_file <- "data/input/mt_annotation.bed"

# =========================
# Load data
# =========================

# mtDNA alignment data
mt_hits <- read.table(mt_hits_file, header = FALSE)
colnames(mt_hits) <- c("numt_id", "flag", "mt_chr", "mt_start", "cigar")

# nuclear genome coordinates (from fasta headers)
headers <- readLines(nuclear_fasta_headers)
headers <- headers[grepl("^>", headers)]

# mtDNA gene annotation
mt_genes <- read.table(mt_annotation_file, header = FALSE)
colnames(mt_genes) <- c("chr", "start", "end", "gene", "score", "strand")

# =========================
# Function: parse CIGAR string
# =========================
get_ref_length <- function(cigar) {
  ops <- str_extract_all(cigar, "[0-9]+[MDN=X]")[[1]]
  lengths <- as.numeric(str_extract(ops, "[0-9]+"))
  return(sum(lengths))
}

# =========================
# Process mtDNA alignments
# =========================

# Clean IDs
mt_hits$numt_id <- str_extract(mt_hits$numt_id, "FRAG_[0-9]+")

# Compute alignment lengths
mt_hits$ref_length <- sapply(mt_hits$cigar, get_ref_length)
mt_hits$mt_end <- mt_hits$mt_start + mt_hits$ref_length - 1

# Normalize coordinates
mt_clean <- mt_hits %>%
  mutate(
    mt_start = pmin(mt_start, mt_end),
    mt_end   = pmax(mt_start, mt_end)
  ) %>%
  select(numt_id, mt_start, mt_end)

# =========================
# Parse nuclear coordinates
# =========================
nuclear <- tibble(header = headers) %>%
  mutate(
    numt_id = str_extract(header, "FRAG_[0-9]+"),
    chr = str_extract(header, "(?<=::)[^:]+"),
    coords = str_extract(header, "[0-9]+-[0-9]+"),
    start = as.numeric(str_extract(coords, "^[0-9]+")),
    end   = as.numeric(str_extract(coords, "[0-9]+$"))
  ) %>%
  select(numt_id, chr, start, end)

# =========================
# Merge datasets
# =========================
final_numts <- nuclear %>%
  inner_join(mt_clean, by = "numt_id") %>%
  filter(!is.na(mt_start) & !is.na(mt_end))

cat("Total NUMTs:", nrow(final_numts), "\n")

# =========================
# Genome setup
# =========================
chrom_sizes <- final_numts %>%
  group_by(chr) %>%
  summarize(length = max(end), .groups = "drop")

mt_size <- max(mt_genes$end)

sectors <- c("mtDNA", chrom_sizes$chr)

xlims <- rbind(
  c(0, mt_size),
  cbind(rep(0, nrow(chrom_sizes)), chrom_sizes$length)
)

# =========================
# Initialize circos
# =========================
circos.clear()

circos.par(
  start.degree = 90,
  gap.degree = 2
)

circos.initialize(
  factors = sectors,
  xlim = xlims
)

# =========================
# Color scheme
# =========================
chrom_colors <- rep("grey70", length(chrom_sizes$chr))
names(chrom_colors) <- chrom_sizes$chr

# Highlight Y chromosome if present
if ("Y" %in% chrom_sizes$chr) {
  chrom_colors["Y"] <- "#D55E00"
}

mt_color <- "#1f3b7a"

# =========================
# Outer ring (chromosomes)
# =========================
circos.trackPlotRegion(
  ylim = c(0, 1),
  track.height = 0.08,
  panel.fun = function(x, y) {
    
    sector <- get.cell.meta.data("sector.index")
    
    circos.rect(
      get.cell.meta.data("xlim")[1], 0,
      get.cell.meta.data("xlim")[2], 1,
      col = ifelse(sector == "mtDNA", mt_color, chrom_colors[sector]),
      border = NA
    )
    
    circos.text(
      get.cell.meta.data("xcenter"),
      1.2,
      labels = sector,
      facing = "inside",
      cex = 0.7
    )
  }
)

# =========================
# mtDNA gene annotation
# =========================
mt_genes$type <- ifelse(grepl("^trn", mt_genes$gene), "tRNA",
                       ifelse(grepl("^rrn", mt_genes$gene), "rRNA", "protein"))

gene_colors <- c(
  "protein" = "#D55E00",
  "rRNA"    = "#E69F00",
  "tRNA"    = "#F0E442"
)

circos.trackPlotRegion(
  factors = "mtDNA",
  ylim = c(0, 1),
  track.height = 0.05,
  bg.border = NA
)

for (i in 1:nrow(mt_genes)) {
  
  circos.rect(
    mt_genes$start[i],
    0,
    mt_genes$end[i],
    1,
    sector.index = "mtDNA",
    col = gene_colors[mt_genes$type[i]],
    border = NA
  )
}

# =========================
# NUMT fragments (stacked)
# =========================
if (nrow(final_numts) > 0) {
  
  final_numts <- final_numts %>% arrange(mt_start)
  final_numts$stack <- seq_len(nrow(final_numts))
  
  circos.trackPlotRegion(
    factors = "mtDNA",
    ylim = c(0, max(final_numts$stack)),
    track.height = 0.1,
    bg.border = NA
  )
  
  for (i in 1:nrow(final_numts)) {
    
    circos.rect(
      final_numts$mt_start[i],
      final_numts$stack[i] - 0.5,
      final_numts$mt_end[i],
      final_numts$stack[i],
      sector.index = "mtDNA",
      col = "#2c7fb8",
      border = NA
    )
  }
}

# =========================
# Links: mtDNA → nuclear genome
# =========================
for (i in 1:nrow(final_numts)) {
  
  circos.link(
    "mtDNA",
    c(final_numts$mt_start[i], final_numts$mt_end[i]),
    final_numts$chr[i],
    c(final_numts$start[i], final_numts$end[i]),
    col = adjustcolor(chrom_colors[final_numts$chr[i]], alpha.f = 0.4),
    border = NA
  )
}

############################################################
# END
############################################################
