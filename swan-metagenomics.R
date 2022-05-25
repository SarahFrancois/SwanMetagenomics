# Swans feces metagenomic sequencing data analysis - R code

# Author: Sarah François
# Date of last change: 25/05/2022


###~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~###
# Set the working directory and data loading ----
###~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~###

#
# Set the working directory
#
setwd("/path/to/the/working/directory")

#
# Data loading
#
input_file <- "/path/to/the/contingency_table.tab" # contingency table: rows=taxa, columns=samples
input_data <-read.table(file=input_file, header=T, row.names = 1,dec=".",sep="\t")
input_data[input_data==""] <- 0 # replace empty cells by 0
metadata<-read.table(file="/path/to/the/metadata.tab", header=T, row.names = 1,dec=".",sep="\t") # metadata table: rows=samples, columns=variables
output_file <- "/path/to/the/working/directory/results.tab"



###~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~###
# Data formatting ----
###~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~###

#
# Data filtering: abundance threshold (taxonomic attribution artifacts removal) 
#

# Set abundance threshold
abundance_threshold <- metadata$number_of_cleaned_reads /1000000 # set abundance threshold per sample = 1/1,000,000 reads for bacteria, 1/10,000,000 reads for viruses
abundance_threshold_rounded <- round(abundance_threshold, digits = 0) # round to the specified number of decimals
# Data filtering
input_data_t <- t(input_data)
input_data_t_ab <-ifelse(input_data_t<abundance_threshold_rounded,0,input_data_t) # replace counts < threshold by 0
filtered_data_t <- input_data_t_ab[, colSums(input_data_t_ab != 0) > 0] # remove empty columns
filtered_data <- t(filtered_data_t)
# Export table
write.table(filtered_data, output_file, sep="\t") # export table

#
# 1- Removal of samples which do not contain any counts
#

# Libraries
library(dplyr)

# Contingency table
filtered_data <- filtered_data[, colSums(filtered_data != 0) > 0]
# Metadata table
metadata_df <- as.data.frame(t(metadata))
filtered_data_df <- as.data.frame(filtered_data)
metadata_t_cleaned <- metadata_df[,(colnames(metadata_df) %in% colnames(filtered_data_df))] # delete uncommon columns between the two files
metadata <- as.data.frame(t(metadata_t_cleaned))

#
# 2 - Removal of samples with no associated metadata
#
VARIABLE <- metadata$VAR # VAR: variable
NA_ROWS <- which(is.na(VARIABLE))
filtered_data_t <- as.data.frame(t(filtered_data))
filtered_data_t_woNA <- filtered_data_t[-NA_ROWS, ] 
metadata_woNA <- metadata[-NA_ROWS, ] 

#
# 3 - Differential abundance and prevalence testing: removal of taxa present in >5 samples
# 
presence_absence_data <- ifelse(filtered_data_t_woNA>0,1,0) # replace counts by 0/1 values
data_sup_5_samples_t <- filtered_data_t_woNA[, colSums(presence_absence_data != 0) > 5]



###~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~###
# Accumulation and Rarefaction curves  ----
###~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~###

# Libraries
library(vegan)

# Rarefaction curves
pdf(file = "rarefaction_curve.pdf")
rarecurve(t(filtered_data), step = 100000, ylab = "number of taxa", xlab = "number of reads", col = "blue", label = FALSE, cex = 1) # step = Step size for sample sizes in rarefaction curves.
dev.off() # save to /path/to/the/working/directory

# Accumulation curve
pdf(file = "accumulation_curve.pdf")
plot(specaccum(t(filtered_data), xlab = "Number of samples", ylab = "Number of taxa"))
dev.off()



###~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~###
# Data summary - basic statistics, abundance visualisation ----
###~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~###

# Load library
library(DescTools)

#
# Taxa summary of abundance and prevalence in samples
#
mean<-apply(filtered_data,1,mean) # mean
  mean_round <- round(mean, digits = 0) # rounded mean
median<-apply(filtered_data,1,median) # median
sd<- apply(filtered_data,1,sd) # sd
  sd_round <- round(sd, digits = 0) # rounded median
min<- apply(filtered_data,1,min) # min
max<- apply(filtered_data,1,max) # max
  # Prevalence and 95% CI 
data_01_counts<-ifelse(filtered_data>0,1,0)
row_sums <- rowSums(data_01_counts)
prevalence_CI_95_lists  <- lapply(1:nrow(data_01_counts), function(x) BinomCI(row_sums[[x]], ncol(data_01_counts), conf.level = 0.95, method = "wald"))
prevalence_CI_95_df <- data.frame(matrix(unlist(prevalence_CI_95_lists), nrow=length(prevalence_CI_95_lists), byrow=T))
prevalence_CI_95_df_rounded <- round(prevalence_CI_95_df, digits = 2) # rounded prevalence, 2 digits
colnames(prevalence_CI_95_df_rounded) <- c("prevalence_estimate", "prevalence_lower_95%CI", "prevalence_upper_95%CI")
  # Table creation 
basic_stats <- cbind(mean_round, sd_round, median, min, max, prevalence_CI_95_df_rounded, deparse.level=1) # tables binding
write.table(basic_stats, output_file, sep="\t") # export table

#
# Abundance visualisation using barplots
#
  # Libraries
library(phyloseq)
library(ggplot2)
library(ggsignif)
  # Construction of a phyloseq object
filtered_data_t <- t(filtered_data)
OTU <- otu_table(filtered_data_t, taxa_are_rows=FALSE)
MET <- sample_data(metadata)
physeq <- phyloseq(OTU, MET)
  # Plot most abundant taxa
par(mar = c(10, 4, 4, 2) + 0.1) # make more room on bottom margin
N <- 30
barplot(sort(taxa_sums(physeq), TRUE)[1:N], las=2)
barplot(sort(taxa_sums(physeq), TRUE)[1:N]/nsamples(physeq), las=2) # divide by number of samples
  # Estimate density of first sample
densityOfSample1 = density(filtered_data_t[,1])
  # plot densitygraph / histogram of abundances
plot(densityOfSample1)
hist(filtered_data_t[,1], nclass = 50)
  # Estimate and plot density distributions of all samples
S_densities = apply(filtered_data_t,2,density)
lapply(S_densities,lines)
S_hist = apply(filtered_data_t,2,hist)
lapply(S_hist)
X11()
  # Create barplot of original data
pdf(file = "Abundance_barplot.pdf")
barplot(filtered_data_t, las = 2)
dev.off()
  # Taxa distribution visualisation
pdf(file = "Abundance_boxplot.pdf")
boxplot(filtered_data_t, las = 2)
dev.off()



###~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~###
# Normalisation  ----
###~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~###

# Launch libraries
library(edgeR)
library("metagenomeSeq")

#
# min-max normalisation 
#
  # normalisation
min <- apply(filtered_data,1,min) # 1 = by row = taxa
max <- apply(filtered_data,1,max) #  = by row = taxa
minmax_normalised <- (filtered_data - min) / (max - min)
  # export table
write.table(minmax_normalised, output_file, sep="\t")

#
# cpm normalisation 
#
  # normalisation
cpm_normalised <- cpm(filtered_data, normalized.lib.sizes = FALSE)
  # export table
write.table(cpm_normalised, output_file, sep="\t") 

#
# Trimmed Mean of M-values (TMM) Normalisation 
#
  # data formatting - DGEList object creation
d <- DGEList(counts=filtered_data)
d <- estimateCommonDisp(d)
  # normalisation
d <- calcNormFactors(d, method = "TMM")
TMM_normalised <- filtered_data * d$samples$norm.factors
TMM_cpm_normalised <- cpm(TMM_normalised, normalized.lib.sizes = FALSE)
  # export table
write.table(TMM_normalised, output_file, sep="\t") # export table TMM normalised
write.table(TMM_cpm_normalised, output_file, sep="\t") # export table TMM + cpm normalised

#
# Relative Log Expression (RLE) Normalisation 
#
  # data formatting - DGEList object creation
filtered_data_2 <- filtered_data
filtered_data_2[filtered_data_2 < 1] <- 1
d <- DGEList(counts=filtered_data_2)
d <- estimateCommonDisp(d)
  # normalisation
d <- calcNormFactors(d, method = "RLE")
RLE_normalised <- filtered_data * d$samples$norm.factors
RLE_cpm_normalised <- cpm(RLE_normalised, normalized.lib.sizes = FALSE)
  # export table
write.table(RLE_normalised, output_file, sep="\t") # export table RLE normalised
write.table(RLE_cpm_normalised, output_file, sep="\t") # export table RLE + cpm normalised

#
# Upper quantile (q75) Normalisation 
#
  # data formatting - DGEList object creation
filtered_data_2 <- filtered_data
filtered_data_2[filtered_data_2 < 1] <- 1
d <- DGEList(counts=filtered_data_2)
d <- estimateCommonDisp(d)
  # normalisation
d <- calcNormFactors(d, method = "upperquartile", p = 0.75)
Q75_normalised <- filtered_data * d$samples$norm.factors
Q75_cpm_normalised <- cpm(Q75_normalised, normalized.lib.sizes = FALSE)
  # export table
write.table(Q75_normalised, output_file, sep="\t") # export table Q75 normalised
write.table(Q75_cpm_normalised, output_file, sep="\t") # export table Q75 + cpm normalised

#
# CSS Normalisation
#
  # data formatting
MRE <- newMRexperiment(filtered_data)
MRE <- cumNorm(MRE, p=cumNormStatFast(MRE))
  # normalisation (remove samples for which colsum = 0 or 1 before normalisation step)
CSS_normalised <- MRcounts(MRE, norm = T)
  # export table
write.table(CSS_normalised, output_file, sep="\t")



###~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~###
# Diversity analysis ----
###~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~###

# Launch libraries
library(phyloseq)
library(ggplot2)
library(ggsignif)
library(vegan)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~# 
# 1- Data formating
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~# 

# Construction of a phyloseq object
OTU <- otu_table(filtered_data_t_woNA, taxa_are_rows=FALSE)
MET <- sample_data(metadata_woNA)
physeq <- phyloseq(OTU, MET)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~# 
# 2- Diversity analyses and visualisation
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~# 

# Alpha diversity indexes and boxplot visualisation
alpha_meas = c("Observed", "Shannon", "Simpson") 
richness <- plot_richness(physeq, x="VARIABLE", measures=alpha_meas) # x= variable 
boxplot <- richness + geom_boxplot(data=richness$data, aes(x=VARIABLE, y=value, color=NULL), alpha=0.1) +
  theme_classic() + 
  scale_x_discrete(limits=c("VARIABLE_VALUE1", "VARIABLE_VALUE2", "VARIABLE_VALUE3", "VARIABLE_VALUE4")) +
  theme(axis.text.x = element_text(angle=45, hjust=1)) +
  geom_signif(test="wilcox.test", comparisons = combn(levels(metadata_woNA$VARIABLE),2, simplify = F),
              step_increase = 0.2)
# Save plot
pdf(file = "Diversity_indexes_boxplot.pdf")
boxplot
dev.off()

# Bray-Curtis beta diversity analysis and visualisation
comm.mds <- metaMDS(comm = filtered_data_t_woNA, distance = "bray", trace = FALSE, autotransform = FALSE) # Bray-Curtis distance measure without transformation of the data
plot(comm.mds$points)
MDS_xy <- data.frame(comm.mds$points)
MDS_xy$VARIABLE <- metadata_woNA$VARIABLE # set variable of interest
ggplot(MDS_xy, aes(MDS1, MDS2, color = VARIABLE)) + geom_point() + stat_ellipse() + theme_bw()


# PERMANOVA significance test for group-level differences based on Bray-Curtis index
library(vegan)
permanova <- adonis(filtered_data_t_woNA ~ VARIABLE, # set variable of interest
                    data = metadata_woNA, permutations=10000, method = "bray")
print(permanova)

#Checking the homogeneity condition (Check that variance homogeneity assumptions hold (to ensure the reliability of the results))
  ## Bray-Curtis distances between samples
dist <- vegdist(filtered_data_t_woNA) 
## Calculate multivariate dispersions
anova(betadisper(dist, metadata_woNA$VARIABLE))
#Investigate the top 10 driving taxa, and show coefficients for the top taxa separating the groups
coef <- coefficients(permanova)["VARIABLE",] 
top.coef <- coef[rev(order(abs(coef)))[1:10]]
par(mar = c(3, 14, 2, 1))
barplot(sort(top.coef), horiz = T, las = 1, main = "Top taxa")

###~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~###
# Differential abundance analysis ----
### Wilcoxon test [1]
### Kruskal-Wallis test [2]
  ### Differential abundance testing 
  ### Data visualisation using boxplots 
###~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~###

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~# 
# 1- Wilcoxon tests
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~# 

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~# 
# Wilcoxon tests - table
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~# 

# Variable selection
VARIABLE <- metadata_woNA$VARIABLE # set variable

# Wilcoxon tests
wilcoxon_test_output <- lapply(data_sup_5_samples_t, function(x) wilcox.test(x ~ VARIABLE, p.adjust.method = "none", mu=0, alt="two.sided", paired=F, conf.int=T, conf.level=0.95, exact=F, correct=F))
names(wilcoxon_test_output) <- colnames(data_sup_5_samples_t)
  # No adjustment for multiple comparisons
wilcoxon_pval_none<- sapply(wilcoxon_test_output, function(x) {x$p.value})
  # Multiple comparisons test adjustments
wilcoxon_pval_adjusted_BY <- p.adjust(wilcoxon_pval_none, method="BY")
  # Confidence interval
wilcoxon_pval_confint <- sapply(wilcoxon_test_output, function(x) {x$conf.int})
wilcoxon_pval_confint_t <- t(wilcoxon_pval_confint)
colnames(wilcoxon_pval_confint_t) <- c("lower_estimate", "upper_estimate")
  # Estimate
wilcoxon_pval_estimate <- sapply(wilcoxon_test_output, function(x) {x$estimate})

# Table creation
wilcoxon_test_results_final <- cbind(wilcoxon_pval_none, wilcoxon_pval_adjusted_BY, wilcoxon_pval_estimate, wilcoxon_pval_confint_t, deparse.level=1) # bind tables
colnames(wilcoxon_test_results_final) <- c("taxa", "wilcoxon_pvalue", "wilcoxon_pvalue_adjusted_BY", "estimate", "lower_interval", "upper_interval")
# Export table
write.table(wilcoxon_test_results_final, output_file, sep="\t")

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~# 
# Data visualisation - boxplots
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~# 

# Launch libraries
library(ggplot2)
library(gplots)
library(ggsignif)
library(dplyr)

# Variable selection
VARIABLE <- metadata_woNA$VARIABLE # set variable

# Data formatting
pseudocount_log10 <- log10(data_sup_5_samples_t +1) # convertion of counts to Log10

# Boxplots creation 
plots = list()
for(nm in names(pseudocount_log10)) {
  plots[[nm]] <- 
    ggplot(pseudocount_log10, aes_string(x=VARIABLE, y=nm)) +
    geom_boxplot(fill="gray") +
    theme_classic() + 
    scale_y_continuous(name="Number of reads") +
    scale_x_discrete(name=nm, limits=c("VARIABLE_VALUE1", "VARIABLE_VALUE2")) +
    geom_jitter(shape=16, position=position_jitter(0.1)) +
    geom_signif(comparison = list(c("VARIABLE_VALUE1","VARIABLE_VALUE2", test = "wilcox.test", map_signif_level=TRUE)))
  ggsave(plots[[nm]], device = "png", filename=paste("/path/to/the/working/directory/boxplot_", nm, sep=""))
}


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~# 
# 2- Kruskal-Wallis tests
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~# 

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~# 
# Kruskal-Wallis tests - table
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~# 

# Variable selection
VARIABLE <- metadata_woNA$VARIABLE # set variable

# Kruskal-Wallis tests
kruskal_wallis_test_output <- lapply(1:ncol(data_sup_5_samples_t), function(x) kruskal.test(data_sup_5_samples_t[[x]], p.adjust.method = "none", VARIABLE))
names(kruskal_wallis_test_output) <- colnames(data_sup_5_samples_t)
  # No adjustment for multiple comparisons
kruskal_wallis_pval_none<- sapply(kruskal_wallis_test_output, function(x) {x$p.value})
  # Multiple comparisons test adjustments
kruskal_wallis_pval_adjusted_BY <- p.adjust(kruskal_wallis_pval_none, method="BY")

# Table creation
kruskal_wallis_test_results_final <- cbind(kruskal_wallis_pval_none, kruskal_wallis_pval_adjusted_BY, deparse.level=1) # bind tables
colnames(kruskal_wallis_test_results_final) <- c("taxa", "kruskal_wallis_pvalue", "kruskal_wallis_pvalue_adjusted_BY")
# Export table
write.table(kruskal_wallis_test_results_final, output_file, sep="\t")


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#    
# Data visualisation - boxplots  
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~# 

# Launch libraries
library(ggplot2)
library(gplots)
library(ggsignif)
library(dplyr)

# Variable selection
VARIABLE <- metadata_woNA$VARIABLE # set variable

# Data formatting
pseudocount_log10 <- log10(data_sup_5_samples_t +1) # convertion of counts to Log10

# Boxplots creation 
plots = list()
for(nm in names(pseudocount_log10)) {
  plots[[nm]] <- 
    ggplot(pseudocount_log10, aes_string(x=VARIABLE, y=nm)) +
    geom_boxplot(fill="gray") +
    theme_classic() + 
    scale_y_continuous(name="Log10 number of reads") +
    scale_x_discrete(name=nm, limits=c("VARIABLE_VALUE1", "VARIABLE_VALUE2", "VARIABLE_VALUE3", "VARIABLE_VALUE4")) +
    geom_jitter(shape=16, position=position_jitter(0.1)) +
    theme(axis.text.x = element_text(angle=45, hjust=1)) +
    geom_signif(test="wilcox.test", comparisons = combn(levels(metadata_woNA$VARIABLE),2, simplify = F),
                step_increase = 0.1)
  ggsave(plots[[nm]], device = "png", filename=paste("/path/to/the/working/directory/boxplot_", nm, sep=""))
}



###~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~###
# Differential prevalence analysis ----
###~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~###

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~# 
# 1- Data formatting
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~# 

# Variable selection
VARIABLE <- metadata_woNA$VARIABLE # set variable

# Creation of a table of presence/absence table per taxa
rownames(data_sup_5_samples_t) <- VARIABLE # variable values as row names
occ_per_var <- as.data.frame(table(VARIABLE)) # number of occurences per variable value - metadata table
occ_per_var_per_taxa_t <- as.data.frame(t(rowsum(data_sup_5_samples_t, VARIABLE))) # number of occurences per variable value per taxa - contingency table
df_VARIABLE_VALUE1 <- occ_per_var_per_taxa_t$VARIABLE_VALUE1
df_VARIABLE_VALUE2 <- occ_per_var_per_taxa_t$VARIABLE_VALUE2
occ_per_var # displays number of samples in VARIABLE_VALUE1 and VARIABLE_VALUE2
absence_VARIABLE_VALUE1 <- lapply(df_VARIABLE_VALUE1, FUN=function(x) "number_of_samples_VARIABLE_VALUE1"-x) # number_of_samples_VARIABLE_VALUE1: in occ_per_var
absence_VARIABLE_VALUE2 <- lapply(df_VARIABLE_VALUE2, FUN=function(x) "number_of_samples_VARIABLE_VALUE2"-x) # number_of_samples_VARIABLE_VALUE2: in occ_per_var
  # Table creation
fisher_table <- cbind(occ_per_var_per_taxa_t$VARIABLE_VALUE2, absence_VARIABLE_VALUE2, occ_per_var_per_taxa_t$VARIABLE_VALUE1, absence_VARIABLE_VALUE1, deparse.level=1) # bind tables
colnames(fisher_table) <- c("presence_VARIABLE_VALUE2", "absence_VARIABLE_VALUE2", "presence_VARIABLE_VALUE1", "absence_VARIABLE_VALUE1")
rownames(fisher_table) <- colnames(data_sup_5_samples_t)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~# 
# 2- Fisher tests
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~# 

# Data formatting
fisher_table_t <- t(fisher_table)
fisher_results_p <- data.frame(matrix(nrow= ncol(fisher_table_t), ncol =1)) # p-values
  colnames(fisher_results_p) <- c("p.value")
  rownames(fisher_results_p) <- colnames(fisher_table_t)
fisher_results_c <- data.frame(matrix(nrow= ncol(fisher_table_t), ncol =2)) # confidence interval
  colnames(fisher_results_c) <- c("confidence_interval_min", "confidence_interval_max")
  rownames(fisher_results_c) <- colnames(fisher_table_t)
fisher_results_e <- data.frame(matrix(nrow= ncol(fisher_table_t), ncol =1)) # estimate
  colnames(fisher_results_e) <- c("estimate")
  rownames(fisher_results_e) <- colnames(fisher_table_t)

# Fisher tests for all taxa
for(nm in 1:ncol(fisher_table_t)) {
  challenge.df = matrix(unlist(fisher_table_t[,nm]), ncol=2)
  fisher_test <- fisher.test(challenge.df)
  fisher_results_p[nm, 1] <- fisher_test$p.value
  fisher_results_c[nm,] <- fisher_test$conf.int 
  fisher_results_e[nm, 1] <- fisher_test$estimate
} 
p.value_BY <- p.adjust(fisher_results_p$p.value, method="BY")

# Table creation
table <- cbind(fisher_table, fisher_results_p, p.value_BY, fisher_results_c, fisher_results_e)
table <- as.matrix(table)

# Export table
write.table(table, output_file, sep="\t")
