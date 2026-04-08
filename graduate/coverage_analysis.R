############################################################
# Coverage Analysis (Male vs Female)
############################################################

library(tidyverse)

process_coverage <- function(male_file, female_file) {
  
  male <- read.table(male_file)
  female <- read.table(female_file)
  
  names(male) <- paste("m", c("chrom","start","end","reads","cov"), sep=".")
  names(female) <- paste("f", c("chrom","start","end","reads","cov"), sep=".")
  
  df <- cbind(male, female) %>%
    filter(m.cov + f.cov > 0) %>%
    mutate(log2_ratio = log2((m.cov + 1) / (f.cov + 1)))
  
  return(df)
}

df <- process_coverage(
  "data/input/male.cov",
  "data/input/female.cov"
)

ggplot(df, aes(x = m.chrom, y = log2_ratio)) +
  geom_boxplot() +
  theme_minimal()
