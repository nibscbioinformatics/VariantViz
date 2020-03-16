args = commandArgs(trailingOnly=TRUE)

vcf_file_path = args[1]
output_name = args[2]

options(warn=-1)

library(vcfR)
library(tidyr)
library(stringr)
library(tidyverse)

vcf <- read.vcfR(vcf_file_path, verbose = FALSE)

# Formating dataframe
ad = as.data.frame(extract.gt(vcf, element = 'AD'))
dp = as.data.frame(extract.gt(vcf, element = 'DP'))
ad_separated = separate(ad, colnames(ad), into = c('REF', 'ALT'), sep = ',')
ad_ref = ad_separated$REF
ref_alt_dataframe = dp
ref_alt_dataframe$ad_ref = ad_ref
colnames(ref_alt_dataframe) = c('DP', 'AD_REF')
ref_alt_dataframe$AD_REF = as.numeric(ref_alt_dataframe$AD_REF)
ref_alt_dataframe$DP = as.numeric(levels(ref_alt_dataframe$DP)[as.numeric(ref_alt_dataframe$DP)])
ref_alt_dataframe$FQ_REF = ref_alt_dataframe$AD_REF / ref_alt_dataframe$DP
ref_alt_dataframe$FQ_ALT = 1 - ref_alt_dataframe$FQ_REF

# Get ClinVar information
vcf_info = getINFO(vcf)
clinvar_annotations = c()
rs_annotations = c()
gene_annotations = c()
for (SNP_info in vcf_info) {
  raw_clinvar_annotation = str_extract(SNP_info, regex('CLNSIG=.*?;{1}'))
  raw_rs_annotation =str_extract(SNP_info, regex('RS=.*'))
  raw_gene_annotation = str_extract(SNP_info, regex('GENEINFO=.*?:{1}'))
  
  clinvar_annotation = substr(raw_clinvar_annotation, 8, nchar(raw_clinvar_annotation) - 1)
  rs_annotation = substr(raw_rs_annotation, 4, nchar(raw_rs_annotation))
  rs_annotation = paste('rs', rs_annotation, sep = '')
  gene_annotation = substr(raw_gene_annotation, 10, nchar(raw_gene_annotation) - 1)
  
  clinvar_annotations = c(clinvar_annotations, clinvar_annotation)
  rs_annotations = c(rs_annotations, rs_annotation)
  gene_annotations = c(gene_annotations, gene_annotation)
}

ref_alt_dataframe$CLNSIG = clinvar_annotations
ref_alt_dataframe$GENE = gene_annotations
ref_alt_dataframe$RS = rs_annotations
ref_alt_dataframe$CLNSIG <- gsub('_', ' ', ref_alt_dataframe$CLNSIG)
ref_alt_dataframe$RS <- gsub('rsNA', 'No rs ID', ref_alt_dataframe$RS)



# Plot Circos plots
data_circos = subset(ref_alt_dataframe, select = -c(DP, AD_REF, FQ_REF))
data_circos = data_circos[complete.cases(data_circos), ]
data_circos$individual = rownames(data_circos)
colnames(data_circos) = c('value', 'group', 'gene', 'rs', 'individual')
rownames(data_circos) = NULL
data = data_circos


# Set a number of 'empty bar' to add at the end of each group
empty_bar <- 4
to_add <- data.frame( matrix(NA, empty_bar*length(unique(data$group)), ncol(data)) )
colnames(to_add) <- colnames(data)
to_add$group <- rep(unique(data$group), each=empty_bar)
data <- rbind(data, to_add)
data = data %>% arrange(match(group, c('Benign', 'Benign/Likely benign', 'Likely benign', 'Conflicting interpretations of pathogenicity', 'Pathogenic', 'Uncertain significance')))
data$id <- seq(1, nrow(data))

# Get the name and the y position of each label
label_data <- data
number_of_bar <- nrow(label_data)
angle <- 90 - 360 * (label_data$id-0.5) /number_of_bar     # I substract 0.5 because the letter must have the angle of the center of the bars. Not extreme right(1) or extreme left (0)
label_data$hjust <- ifelse(angle < -90, 1, 0)
label_data$angle <- ifelse(angle < -90, angle+180, angle)

# Make the plot
p <- ggplot(data, aes(x=as.factor(id), y=value*100, fill=group)) +       # Note that id is a factor. If x is numeric, there is some space between the first bar
  geom_bar(stat="identity", alpha=0.5) +
  ylim(-100,120) +
  theme_minimal() +
  theme(
    legend.title = element_text(size = 20),
    legend.text = element_text(size = 13),
    axis.text = element_blank(),
    axis.title = element_blank(),
    panel.grid = element_blank()
  ) + geom_text(aes(x=0, y=-100, label = dim(data_circos)[1]), color='black', size = 40) +
  
  coord_polar() +
  geom_text(data=label_data, aes(x=id, y=120, label=gene, hjust=hjust), color="black", fontface="bold",alpha=0.6, size=3.5, angle= label_data$angle, inherit.aes = FALSE ) +
  ggtitle(output_name) +  theme(plot.title = element_text(hjust = 0.5, size = 50))


# Export data

pdf(paste(output_name, '_SNP_inhouse_circos_plot.pdf', sep = ''), width = 16, height = 8)
p + labs(fill = 'ClinVar tiers:') + scale_fill_manual(breaks=c('Benign', 'Benign/Likely benign', 'Likely benign', 'Conflicting interpretations of pathogenicity', 'Pathogenic', 'Uncertain significance'), 
                                                      values=c('Benign' = 'green4', 'Benign/Likely benign'= 'turquoise', 'Likely benign'= 'steelblue', 'Conflicting interpretations of pathogenicity'= 'yellow2', 'Pathogenic'= 'red3', 'Uncertain significance'= 'grey'))
dev.off()

data = data[complete.cases(data),]
data$id = output_name

write_tsv(data, path = paste(output_name, '_SNP_inhouse_table.tsv', sep = ''))
