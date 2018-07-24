library(tidyverse)
library(rlang)

# generate data using:
# bcftools query --samples CB4856,CX11314,DL238,ED3017,EG4725,JT11398,JU258,JU775,LKC34,MY16,MY23,N2 WI.20180510.soft-filter.vcf.gz -f '[%CHROM\t%POS\t%REF\t%ALT\t%QUAL\t%FILTER\t%AC\t%AF\t%SAMPLE\t%GT\t%DP\t%FT\t%AD\n]' > Soft_Filter_Training.tsv
# awk '$2 > 1e6 && $2 < 2e6 { print }' Soft_Filter_Training.tsv > subset_Soft_Filter_Training.tsv


# set working directory
working.dir <- '~/github_repos/Computational-Bootcamp/'
setwd(working.dir)

# load data
genomics_data <- readr::read_tsv(glue::glue("{working.dir}/Data/subset_Soft_Filter_Training.tsv.gz"),
                                 col_names = F)

# add column names
colnames(genomics_data) <- c("chrom", "pos", "ref", "alt", "qual", "variant_filter", "allele_ct", "allele_freq", "strain", "genotype", "depth", "strain_filter", "allele_depth")

# inspect data
head(genomics_data)
str(genomics_data)

# generate a smaller data set for faster processing
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
  dplyr::filter_all(all_vars(. != "./."))

# what strains are most similar to each other?

geno_matrix[geno_matrix == "0/0"] <- 0
geno_matrix[geno_matrix == "1/1"] <- 1

geno_pc_input <- geno_matrix %>%
  modify_at(c(2:ncol(geno_matrix)), as.numeric) %>%
  dplyr::select(-marker)

pc_results <- prcomp(t(geno_pc_input))

pc_df <- data.frame(strain = row.names(pc_results$x), pc_results$x)

ggplot(pc_df) +
  aes(x = PC1, y = PC2) +
  geom_point() +
  ggrepel::geom_label_repel(aes(label=strain)) +
  theme_bw(15)

# what strains have the most variants?

# facetted
indel_snv_variants %>%
  dplyr::filter(GT != "./.", GT != "0/0") %>%
  dplyr::group_by(strain, variant_type) %>%
  dplyr::summarise(n_v = n()) %>%
  dplyr::ungroup() %>%
  dplyr::arrange(variant_type, n_v) %>%
  dplyr::mutate(strain_f = factor(strain, 
                                  levels = unique(strain), 
                                  labels = unique(strain))) %>%
  ggplot() +
  aes(x = strain_f, y = n_v) +
  geom_bar(stat="identity")+
  facet_grid(variant_type~., scale = "free") +
  labs(x = "Strain", y = "Count") +
  theme_bw(15)

# stacked bar
indel_snv_variants %>%
  dplyr::filter(GT != "./.", GT != "0/0") %>%
  dplyr::group_by(strain, variant_type) %>%
  dplyr::summarise(n_v = n()) %>%
  dplyr::group_by(strain) %>%
  dplyr::mutate(tot_var = sum(n_v)) %>%
  dplyr::arrange(tot_var) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(strain_f = factor(strain, 
                                  levels = unique(strain), 
                                  labels = unique(strain))) %>%
  ggplot() +
  aes(x = strain_f, y = n_v, fill = variant_type) +
  geom_bar(stat="identity", color = "black")+
  scale_fill_manual(values = c("hotpink3", "gray50"), name = "Variant\nType") +
  labs(x = "Strain", y = "Count") +
  theme_bw(15)


# what chromosome has the most variation?

indel_snv_variants %>%
  dplyr::filter(GT != "./.", GT != "0/0") %>%
  dplyr::distinct(chrom, pos, ref, alt, .keep_all = T) %>%
  dplyr::group_by(chrom) %>%
  dplyr::summarise(ct = n())

# what is the correlation between SNV and indel ct?

indel_snv_variants %>%
  dplyr::filter(GT != "./.", GT != "0/0") %>%
  dplyr::group_by(strain, variant_type) %>%
  dplyr::summarise(n_v = n()) %>%
  dplyr::group_by(strain) %>%
  tidyr::spread(variant_type, n_v) %>%
  ggplot() +
  aes(x = SNV, y = indel)+
  geom_point() +
  theme_bw(15) +
  ggrepel::geom_label_repel(aes(label=strain)) +
  labs(x = "SNV Count", y = "Indel Count")

# add linear model equation on plot

snv_indel_cts <- indel_snv_variants %>%
  dplyr::filter(GT != "./.", GT != "0/0") %>%
  dplyr::group_by(strain, variant_type) %>%
  dplyr::summarise(n_v = n()) %>%
  dplyr::group_by(strain) %>%
  tidyr::spread(variant_type, n_v) %>%
  na.omit()

mod <- lm(indel~SNV,data=snv_indel_cts)
s_i_eq <- transform(snv_indel_cts, Fitted = fitted(mod))
eq <- substitute(italic(r)~"="~rvalue*","~italic(p)~"="~pvalue, 
                 list(rvalue = sprintf("%.2f",sign(coef(mod)[2])*sqrt(summary(mod)$r.squared)), 
                      pvalue = format(summary(mod)$coefficients[2,4], digits = 3)))

dftext <- data.frame(SNV = 1000, indel = 750, eq = as.character(as.expression(eq)))


ggplot(snv_indel_cts) +
  aes(x = SNV, y = indel)+
  geom_point() +
  theme_bw(15) +
  ggrepel::geom_label_repel(aes(label=strain)) +
  labs(x = "SNV Count", y = "Indel Count") +
  geom_text(aes(label = eq), data = dftext, parse = TRUE, size  = 6)

# what do the distribution of ALT calls for each strain?

indel_snv_variants %>%
  dplyr::filter(GT != "./.", GT != "0/0") %>%
  ggplot() +
  aes(x = pos/1e6, y = strain, color = variant_type, size = variant_type) +
  geom_point() +
  scale_color_manual(values = c("red", "black"), name = "Variant\nType") +
  scale_size_manual(values = c(0.3, 0.1), name = "Variant\nType") +
  facet_grid(chrom~.) +
  theme_bw(15)

# what variants are unique to a given strain?

# how many transition vs transversion variants are there?

  
  
  
  
  create_mapper_gt <- function(.p){
    glue::glue("~ ({f_text(.p)})") %>% 
      as.formula() %>%
      as_mapper()
  }

create_mapper_gt(~ .x == "0/0")

gt_fix <- function(vec, .p) {
  modify_if(vec, create_mapper_gt(.p) , ~ 0) %>% 
    reduce(c)
}

geno_matrix$CB4856[1:10]

gt_fix(geno_matrix$CB4856[1:10], ~ .x == "0/0")

gt_fix_df <- function(tbl, .p) {
  map_df(tbl, ~ gt_fix(.x, .p) )
} 

geno_matrix[1:10,]

gt_fix_df(geno_matrix[1:10,], ~ .x == "0/0")
