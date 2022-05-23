###############################################################################
##
## Differential Analysis for Apsis Project
##
## Author: Oliver Hunewald
##
## This workflow has been adapted from its original version developed 
## by Nowicka et al.
##
## https://f1000research.com/articles/6-748/v2
##
###############################################################################

library(readxl)
library(ggplot2)
library(reshape2)
library(dplyr)
library(matrixStats)
library(RColorBrewer)
library(pheatmap)
library(mvtnorm)
library(multcomp)
library(emmeans)
library(lme4)
library(ggrepel)


# add function helper script
source("myWrapperFunctions.R")
# set your data directory
setwd("your_data_path")
## ----------------------------------------------------------------------------
## Loading the metadata
## ----------------------------------------------------------------------------
md <- read_excel("Metadata.xlsx")
## Define colors for conditions
color_conditions <- c("#6A3D9A", "#FF7F00", "#008000", "#6A0D9A", "#FF7A00", "#0080A0")
names(color_conditions) <- levels(md$condition)
## Make sure condition variables are factors with the right levels
md$condition <- factor(md$condition, levels = c("yes", "none"))
md$timepoint <- factor(md$timepoint)

head(data.frame(md))

# load the tables clustersizes and proportions
df_counts <- read.csv("cluster_counts_1024.csv")
df_counts$X <- NULL

# load the cluster merging 
cm <- read.csv("annotated_df_1024.csv")
cm$X <- NULL

df_counts['cell'] <- cm$cell
# merge and aggregate by cell
df_counts <- aggregate(. ~ cell, df_counts, sum)

rownames(df_counts) <- df_counts$cell
df_counts$cell <- NULL

props_table <- t(t(df_counts) / colSums(df_counts)) * 100

# use a copy of the original propotion table for later use
my_props <- as.data.frame(props_table)

counts <- as.data.frame.matrix(df_counts)
props <- as.data.frame.matrix(props_table)

props_table <- melt(data.frame(cluster = rownames(props), props),
             id.vars = "cluster", value.name = "proportion", variable.name = "sample_id")

mm <- match(props_table$sample_id, md$sample_id)
props_table$condition <- factor(md$condition[mm])
props_table$timepoint <- factor(md$timepoint[mm])
props_table$patient_id <- factor(md$patient_id[mm])

df_prop <- props_table

mm <- unique(df_prop$sample_id)
mm <- na.omit(mm)
mm<-as.character(mm)

## ----------------------------------------------------------------------------
## Model
## 
## This is an example of sub-setting for BAT,
## comparison high vs. low at timepoint 0
## 
## any other comparison should be adopted below
## ----------------------------------------------------------------------------
library(limma)

head(df_counts)

md$BAT <- factor(md$BAT)
levels(md$BAT)
# subset the metadata
md_subset <- md[md$BAT %in% c("high", "low"), ]
md_subset_t0 <- md_subset[md_subset$timepoint == "t0" & md_subset$condition == "yes", ]

my_props_sub <- my_props[, md_subset_t0$sample_id]

df_counts_sub <- df_counts[md_subset_t0$sample_id]
md_subset_t0$BAT <- factor(md_subset_t0$BAT)


model.matrix( ~ BAT, data = md_subset_t0)
# example for contrast matrix
K <- rbind("HighvsLow" = c(-1, 1))

FDR_cutoff <- 0.05

# formula_glmer_binomial2 <- y/total ~ condition + (1|patient_id) + (1|sample_id)
formula_glmer_binomial2 <- y/total ~ BAT  + (1|sample_id)

counts <- df_counts_sub
ntot <- colSums(counts)

fit_binomial <- lapply(1:nrow(counts), function(i){

  data_tmp <- data.frame(y = as.numeric(counts[i, md_subset_t0$sample_id]),
                         total = ntot[md_subset_t0$sample_id], md_subset_t0)
  
  fit_tmp <- glmer(formula_glmer_binomial2, weights = total, family = binomial, 
                   data = data_tmp)

  ## Fit contrasts one by one
  out <- apply(K, 1, function(k){ #what is K?
    
    contr_tmp <- glht(fit_tmp, linfct = mcp(BAT = k))
    summ_tmp <- summary(contr_tmp, test = adjusted("none"))
    out <- c(summ_tmp$test$coefficients, summ_tmp$test$pvalues)
    names(out) <- c("coeff", "pval")
    return(out)
  })

  return(t(out))
})

### Extract fitted contrast coefficients
coeffs <- lapply(fit_binomial, function(x){
  x[, "coeff"]
})
coeffs <- do.call(rbind, coeffs)
colnames(coeffs) <- paste0("coeff_", rownames(K))
rownames(coeffs) <- rownames(counts)

### Extract p-values
pvals <- lapply(fit_binomial, function(x){
  x[, "pval"]
})
pvals <- do.call(rbind, pvals)
colnames(pvals) <- paste0("pval_", rownames(K))
rownames(pvals) <- rownames(counts)

### Adjust the p-values
adjp <- apply(pvals, 2, p.adjust, method = "BH")
colnames(adjp) <- paste0("adjp_", rownames(K))

da_out1 <- list(coeffs = coeffs, pvals = pvals, adjp = adjp)

## ----------------------------------------------------------------------------
## ----------------------------------------------------------------------------
apply(da_out1$adjp < FDR_cutoff, 2, table)

# props <- df_prop_subset
da_output <- data.frame(cluster = rownames(props),  
                        props, 
                        da_out1$coeffs,
                        da_out1$pvals, 
                        da_out1$adjp, 
                        row.names = NULL)

print(head(da_output), digits = 2)

counts_table <- counts
## Apply the arcsine-square-root transformation
asin_table <- asin(sqrt((t(t(counts_table) / colSums(counts_table)))))
asin <- as.data.frame.matrix(asin_table)

head(da_out1$adjp)
s <- "adjp_HighvsLow"

sign_clusters <- names(which(sort(da_out1$adjp[, s]) < FDR_cutoff))
## Get the adjusted p-values
sign_adjp <- da_out1$adjp[sign_clusters , s, drop=FALSE]
## Keep samples for condition A and normalize to mean = 0 and sd = 1
asin_norm <- as.data.frame(normalization_wrapper(asin[sign_clusters, ]))

plot_differential_heatmap_wrapper(expr_norm = asin_norm, 
                                  sign_adjp = sign_adjp,
                                  condition = md_subset_t0$BAT, 
                                  color_conditions = color_conditions)


props <- my_props_sub


ggdf <- melt(data.frame(cluster = rownames(props), props),
             id.vars = "cluster", value.name = "proportion", variable.name = "sample_id")
ggdf$cluster <- factor(ggdf$cluster)
## Add condition info
mm <- match(ggdf$sample_id, md$sample_id)
ggdf$Param <- factor(md$BAT[mm])

## ----------------------------------------------------------------------------
## Relative abundance of the populations in each sample
## ----------------------------------------------------------------------------
ggdf$patient_id <- factor(md$patient_id[mm])
ggplot(ggdf) +
  geom_boxplot(aes(x = Param, y = proportion, color = Param,
                   fill = Param), position = position_dodge(), alpha = 0.5,
               outlier.color = NA) +
  geom_point(aes(x = Param, y = proportion, color = Param), alpha = 0.8, position = position_jitterdodge()) +
  facet_wrap(~ cluster, scales = "free", nrow = 6) +
  theme_bw() +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(),
        axis.title.x = element_blank(), strip.text = element_text(size = 6)) +
  scale_color_manual(values = color_conditions) +
  scale_fill_manual(values = color_conditions) +
  scale_shape_manual(values = c(16, 17, 8, 3, 12, 0, 1, 2))

