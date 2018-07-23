library(tidyverse)

# set working directory
working.dir <- '~/AndersenLab/Github_Repos/Computational-Bootcamp/'
setwd(working.dir)

# load data
genomics_data <- readr::read_tsv(glue::glue("{working.dir}/Data/subset_Soft_Filter_Training.tsv"),
                                 col_names = F)

# make column names
colnames(genomics_data) <- c("chrom", "pos", "ref", "alt", "qual", "variant_filter", "allele_ct", "allele_freq", "strain", "genotype", "depth", "strain_filter", "allele_depth")

# 
head(genomics_data)

smaller_data <- dplyr::filter(genomics_data, chrom %in% c("IV","V"))

dir.create("Processed_Data")

save(smaller_data, file = glue::glue("{working.dir}/Processed_Data/smaller_dataset.Rda"))

rm(genomics_data)

# how many strains do we have in the data set 

length(unique(smaller_data$strain))

# remove variants that have high levels of heterozygosity
unique(smaller_data$variant_filter)

no_het_variants <- smaller_data %>%
  dplyr::ungroup() %>%
  dplyr::filter(variant_filter != "high_heterozygosity")

# remove multi-allelic variants

biallelic_variants <- no_het_variants %>%
  dplyr::ungroup() %>%
  dplyr::filter(!grepl(",", ref), !grepl(",", alt))

# remove multi-allelic variants

passing_variants <- biallelic_variants %>%
  dplyr::mutate(GT = ifelse(strain_filter == "PASS", genotype, "./.")) %>%
  dplyr::select(-genotype)


# remove variants not present in 12 strains

variable_variants <- passing_variants %>%
  dplyr::group_by(chrom, pos, ref, alt) %>%
  dplyr::filter(sum(grepl("1/1", GT)) > 1) %>%
  dplyr::ungroup()

# re-calculate allele frequencies for remaining strains

n_strains <- length(unique(variable_variants$strain))

passing_variants_af <- variable_variants %>%
  dplyr::select(-allele_freq) %>%
  dplyr::group_by(chrom, pos, ref, alt) %>%
  dplyr::mutate(allele_freq = sum(grepl("1/1", GT))/n_strains)

# plot allele frequency

hist(passing_variants_af$allele_freq)

ggplot(passing_variants_af) + 
  aes(x = allele_freq) +
  geom_histogram() +
  theme_bw(18) + 
  labs(x = "Allele Frequency", y = "Count") + 
  NULL

# split allele depth column into two columns - ref_depth and alt_depth

allele_dp_variants <- passing_variants_af %>%
  dplyr::ungroup() %>%
  tidyr::separate(allele_depth, into = c("ref_depth","alt_depth"), sep = ",", remove = F) %>%
  dplyr::mutate(ref_depth = ifelse(is.na(alt_depth), depth, ref_depth),
                alt_depth = ifelse(is.na(alt_depth), 0, alt_depth))

# plot depth per strain

ggplot(allele_dp_variants) +
  aes(x = pos/1e6, y = as.numeric(depth)) +
  geom_line() +
  facet_grid(strain~chrom) +
  theme_bw(15) +
  labs(x = "Genomic Position (Mb)", y = "Depth") + 
  NULL

# make new column that differentiates SNVs from indels

indel_snv_variants <- allele_dp_variants %>%
  dplyr::mutate(variant_type = ifelse(nchar(ref) == nchar(alt), "SNV", "indel"))

# print out counts for indels and SNVs
indel_snv_variants %>%
  dplyr::distinct(chrom, pos, ref, alt, .keep_all = T) %>%
  group_by(variant_type) %>%
  dplyr::summarise(variant_ct = n())
  
# plot depth per strain by variant type

ggplot(indel_snv_variants) +
  aes(x = pos/1e6, y = as.numeric(depth), color = variant_type) +
  geom_line() +
  scale_color_manual(values = c("cadetblue3", "hotpink3"), name = "Variant\nType")+
  facet_grid(strain~chrom) +
  theme_bw(15) +
  labs(x = "Genomic Position (Mb)", y = "Depth") + 
  NULL

# generate a genotype matrix

geno_matrix <- indel_snv_variants %>%
  tidyr::unite(col = "marker", chrom, pos) %>%
  dplyr::select(marker, strain, GT) %>%
  dplyr::distinct(marker, strain, .keep_all = T) %>% 
  tidyr::spread(strain, GT) %>%
  








# what strains have the most variants?
# what chromosome has the most variation?
# what strains are most similar to each other?
# how many transition vs transversion variants are there?
# what do the distribution of ALT calls for each strain?







