#load packages

#if (!require("BiocManager", quietly = TRUE))
 # install.packages("BiocManager")

#BiocManager::install("phyloseq")

library(tidyverse)
library(phyloseq)
library(csv)

#set seed for reproducible data
set.seed(81)

#metadata
sdata <- as.csv("/Users/lgschaer/Desktop/Permafrost_180/SIP_16S/metadata_v2.csv", row.names = 1, header = TRUE, sep = ",", check.names = TRUE, stringsAsFactors = TRUE) %>%
  mutate(Fraction = factor(Fraction, levels = c("H", "MH", "M", "ML", "L")),
         Day = as.character(Day),
         Temperature = as.character(Temperature))
head(sdata)

#sequence table
sequence_table <- read_table("/Users/lgschaer/Desktop/Permafrost_180/SIP_16S/qiime_out/sequence-table.txt", skip = 1, col_names = TRUE) %>% 
  column_to_rownames(var = "SampleID") %>%
  t() 
sequence_table[1:5,1:5] 

#taxa table
taxa <- read_table("/Users/lgschaer/Desktop/Permafrost_180/SIP_16S/qiime_out/taxonomy.tsv") %>%
  select(Feature, ID) %>%
  separate(ID, into = c(NA, "Kingdom", NA, "Phylum", NA, "Class", NA, "Order", NA, "Family", NA, "Genus", NA, "Species"), sep = ";|__") %>%
  column_to_rownames(var = "Feature") %>%
  as.matrix()
head(taxa)
#View(taxa)

#format for phyloseq
samdata = sample_data(sdata)
seqtab = otu_table(sequence_table, taxa_are_rows = FALSE)
taxtab = tax_table(taxa)

ps = phyloseq(otu_table(seqtab), tax_table(taxtab), sample_data(samdata))
ps

write_rds(ps, "/Users/lgschaer/Desktop/Permafrost_180/SIP_16S/phyloseq_05132024.rds")

###### Phyloseq object is made! On to alpha diversity
library(FSA)

head(sample_data(ps))

#sample colors
sample_colors <- c("lightblue", "maroon", "olivedrab", "purple", "orange")

#plot
ps %>%                                                              #phyloseq object
  plot_richness(
    x = "Fraction",                                                    #compare diversity of datatype
    measures = c("Observed", "Shannon")) +                              #choose diversity measures
  geom_boxplot(aes(fill = Fraction), show.legend = FALSE)+             #make violin plot, set fill aes to sampletype
  theme_linedraw()+                                                     #change theme to classic
  xlab(NULL)+                                                           #no label on x-axis
  theme(axis.text.y.left = element_text(size = 20),                     #adjust y-axis text
        axis.text.x = element_text(size = 20, hjust = 0.5, vjust = 1, angle = 0),           #adjust x-axis label position
        axis.title.y = element_text(size = 20))+                        #adjust y-axis title
  theme(strip.text = element_text(face = "bold", size = 25))+           #adjust headings
  scale_fill_manual(values = sample_colors)+                            #set fill colors
  theme(plot.title=element_text(size = 25, face = "bold", hjust = 0.5)) #change title size, face and position



ps %>%                                                              #phyloseq object
  plot_richness(
    x = "Isotope",                                                    #compare diversity of datatype
    measures = c("Observed", "Shannon")) +                              #choose diversity measures
  geom_boxplot(aes(fill = Isotope), show.legend = FALSE)+             #make violin plot, set fill aes to sampletype
  theme_linedraw()+                                                     #change theme to classic
  xlab(NULL)+                                                           #no label on x-axis
  theme(axis.text.y.left = element_text(size = 20),                     #adjust y-axis text
        axis.text.x = element_text(size = 20, hjust = 0.5, vjust = 1, angle = 0),           #adjust x-axis label position
        axis.title.y = element_text(size = 20))+                        #adjust y-axis title
  theme(strip.text = element_text(face = "bold", size = 25))+           #adjust headings
  scale_fill_manual(values = sample_colors)+                            #set fill colors
  theme(plot.title=element_text(size = 25, face = "bold", hjust = 0.5)) #change title size, face and position

ps %>%                                                              #phyloseq object
  plot_richness(
    x = "Day",                                                    #compare diversity of datatype
    measures = c("Observed", "Shannon")) +                              #choose diversity measures
  geom_boxplot(aes(fill = Day), show.legend = FALSE)+             #make violin plot, set fill aes to sampletype
  theme_linedraw()+                                                     #change theme to classic
  xlab(NULL)+                                                           #no label on x-axis
  theme(axis.text.y.left = element_text(size = 20),                     #adjust y-axis text
        axis.text.x = element_text(size = 20, hjust = 0.5, vjust = 1, angle = 0),           #adjust x-axis label position
        axis.title.y = element_text(size = 20))+                        #adjust y-axis title
  theme(strip.text = element_text(face = "bold", size = 25))+           #adjust headings
  scale_fill_manual(values = sample_colors)+                            #set fill colors
  theme(plot.title=element_text(size = 25, face = "bold", hjust = 0.5)) #change title size, face and position

ps %>%                                                              #phyloseq object
  plot_richness(
    x = "Temperature",                                                    #compare diversity of datatype
    measures = c("Observed", "Shannon")) +                              #choose diversity measures
  geom_boxplot(aes(fill = Temperature), show.legend = FALSE)+             #make violin plot, set fill aes to sampletype
  theme_linedraw()+                                                     #change theme to classic
  xlab(NULL)+                                                           #no label on x-axis
  theme(axis.text.y.left = element_text(size = 20),                     #adjust y-axis text
        axis.text.x = element_text(size = 20, hjust = 0.5, vjust = 1, angle = 0),           #adjust x-axis label position
        axis.title.y = element_text(size = 20))+                        #adjust y-axis title
  theme(strip.text = element_text(face = "bold", size = 25))+           #adjust headings
  scale_fill_manual(values = sample_colors)+                            #set fill colors
  theme(plot.title=element_text(size = 25, face = "bold", hjust = 0.5)) #change title size, face and position

#### Beta Diversity

library(vegan)
library(csv)
library(ape)
library(pairwiseAdonis)


#Bray curtis PCoA ordination

bray <- ordinate(
  physeq = ps, #change this to your phyloseq
  method = "PCoA", 
  distance = "bray" 
)


#plot
colors <- c("lightblue", "maroon", "olivedrab", "purple", "orange")

plot_ordination(
  physeq = ps,                                                          #phyloseq object
  ordination = bray)+                                                        #ordination
  geom_point(aes(fill = Fraction, shape = Isotope), size = 6) +                         
  scale_fill_manual(values = colors) +
  scale_shape_manual(values = c(21, 22, 23, 24, 25))+
  theme_linedraw() +                                                         #changes theme, removes grey background
  theme(                             
    legend.title = element_blank(),                                          #removes legend title
    legend.position = "bottom",
    legend.text = element_text(size = 20, face = "bold"),                                 
    axis.text.y.left = element_text(size = 10),
    axis.title.y = element_text(size = 20),
    axis.text.x = element_text(size = 10),
    axis.title.x = element_text(size = 20),
    strip.text = element_text(face = "bold", size = 20))+
  guides(fill = guide_legend(override.aes = list(shape = 21)))          #fills legend points based on the fill command


# Unifrac PCoA ordination

random_tree = rtree(ntaxa(ps), rooted=TRUE, tip.label=taxa_names(ps))
plot(random_tree)

ps2 = merge_phyloseq(ps, random_tree)
ps2


#ordination
uni <- ordinate(
  physeq = ps2, 
  method = "PCoA", 
  distance = "unifrac"
)
#summary(distance)
#distance

#plot
plot_ordination(
  physeq = ps2,                                                          #phyloseq object
  ordination = uni)+                                                #ordination
  geom_point(aes(fill = Fraction, shape = Isotope), size = 6) +                         #sets fill color to sampletype
  scale_fill_manual(values = colors) +
  scale_shape_manual(values = c(21, 22, 23, 24, 25))+
  theme_linedraw() +                                                      #changes theme, removes grey background
  theme(                             
    legend.title = element_blank(),                                      #removes legend title
    #legend.background = element_rect(fill = "white", color = "black"),  #adds black boarder around legend
    legend.position = "bottom",
    legend.text = element_text(size = 20, face = "bold"),                                 
    axis.text.y.left = element_text(size = 10),
    axis.title.y = element_text(size = 20),
    axis.text.x = element_text(size = 10),
    axis.title.x = element_text(size = 20),
    strip.text = element_text(face = "bold", size = 20))+
  guides(fill = guide_legend(override.aes = list(shape = 21)))          #fills legend points based on the fill command


# MAKING A TAXA PLOT
genusabundance <- ps %>%
  tax_glom(taxrank = "Genus") %>%                        # Set to smallest taxonomic level you are interested in
  transform_sample_counts(function(x) {x/sum(x)} ) %>%   # Transform to rel. abundance
  psmelt()                                               # Melt to long format
head(genusabundance)

family <- genusabundance %>%
  mutate(Order = as.character(Order)) %>%
  select(Fraction, Temperature, Isotope, Day, Order, Abundance) %>%
  group_by(Fraction, Temperature, Isotope, Day) %>%
  mutate(totalSum = sum(Abundance)) %>%
  ungroup() %>%
  group_by(Fraction, Temperature, Isotope, Day, Order, totalSum) %>%
  summarise(
    Abundance = sum(Abundance),
    Order = ifelse(Abundance < 0.05, "< 5 %", Order)) %>%               #change Genus label to group low abundance taxa together
  group_by(Fraction, Temperature, Isotope, Day, Order, totalSum) %>%  #now group and summarize again to group newly labeled low abundance taxa together
  summarise(
    Abundance = sum(Abundance),
    RelAb = Abundance/totalSum) %>%
  unique()
head(family)

colFunc <- colorRampPalette(c("purple4", "lightpink", "firebrick", "orange", "lemonchiffon", "olivedrab4", "darkcyan", "lightblue", "darkblue"))
color_list <- colFunc(length(unique(family$Order)))
genus_colors <- c("black", color_list)
length(genus_colors)

ggplot(family)+
  geom_col(mapping = aes(x = Fraction, y = RelAb, fill = Order), color = "black", position = "stack", show.legend = TRUE)+
  facet_grid(cols = vars(Isotope, Day, Temperature))+
  ylab("Proportion of Community") +
  xlab(NULL)+
  scale_fill_manual(values = genus_colors) +
  theme_linedraw()+
  theme(axis.text.y = element_text(size = 20, color = "black"),
        axis.title.y = element_text(size = 18, color = "black"),
        axis.text.x = element_text(size = 20, angle = 0, vjust = 1, hjust = 0.5, color = "black"),
        legend.text = element_text(size = 13),
        legend.position = "bottom",
        legend.spacing.x = unit(0.1, 'mm'),
        legend.spacing.y = unit(0.05, 'mm'),
        plot.margin=grid::unit(c(0.1,0.1,0.1,0.1), "mm"),
        strip.text = element_text(size = 18, face = "bold", angle = 0),
        legend.title = element_text(face="bold", size = 22))+
  guides(fill=guide_legend(ncol=6,byrow=TRUE))


test <- genusabundance %>%
  filter(Fraction == "H" | Fraction == "L") %>%
  filter(Day == "184") %>%
  filter(Isotope == "18O") %>%
  filter(Temperature == "-1.5") %>%
  filter(Genus == "Methanobacterium_B") %>%
  group_by(Sample) %>%
  summarise(RelAb = sum(Abundance)/n())
head(test)

## taxonomy heatmap

list <- filter(family, RelAb > 0.01)
head(list$Order)

family <- genusabundance %>%
  mutate(Order = as.character(Order)) %>%
  select(Fraction, Temperature, Isotope, Day, Phylum, Order, Abundance) %>%
  group_by(Fraction, Temperature, Isotope, Day) %>%
  mutate(totalSum = sum(Abundance)) %>%
  ungroup() %>%
  group_by(Fraction, Temperature, Isotope, Day, Phylum, Order, totalSum) %>%
 # summarise(
    #Abundance = sum(Abundance),
    #Order = ifelse(Abundance < 0.05, "< 5 %", Order)) %>%               #change Genus label to group low abundance taxa together
  #group_by(Fraction, Temperature, Isotope, Day, Order, totalSum) %>%  #now group and summarize again to group newly labeled low abundance taxa together
  summarise(
    Abundance = sum(Abundance),
    RelAb = Abundance/totalSum) %>%
  unique() %>%
  filter(Day == "184" & Temperature == "-1.5") %>%
  filter(Order %in% list$Order)
head(family)

hist(family$RelAb)

ggplot(family)+
  geom_tile(aes(x = Isotope, y = Order, fill = RelAb))+
  #geom_col(mapping = aes(x = Fraction, y = RelAb, fill = Order), color = "black", position = "stack", show.legend = TRUE)+
  facet_grid(cols = vars(Fraction), rows = vars(Phylum), shrink = TRUE, scales = "free", space = "free")+
  ylab("Proportion of Community") +
  xlab(NULL)+
  scale_fill_gradient2(high = "orange", mid = "firebrick", low = "darkblue", midpoint = 0.15) +
  theme_linedraw()+
  theme(axis.text.y = element_text(size = 18, color = "black"),
        axis.title.y = element_blank(), #text(size = 18, color = "black"),
        axis.text.x = element_text(size = 18, angle = 0, vjust = 1, hjust = 0.5, color = "black"),
        legend.text = element_text(size = 18),
        legend.position = "bottom",
        #legend.spacing.x = unit(0.1, 'mm'),
        #legend.spacing.y = unit(0.05, 'mm'),
        #plot.margin=grid::unit(c(0.1,0.1,0.1,0.1), "mm"),
        strip.text.x = element_text(size = 18, face = "bold", angle = 0),
        strip.text.y = element_text(size = 18, face = "bold", angle = 0),
        legend.title = element_text(face="bold", size = 22))#+
  #guides(fill=guide_legend(ncol=3,byrow=TRUE))


### Differential Abundance

#install DESeq2: https://bioconductor.org/packages/release/bioc/html/DESeq2.html
#if (!requireNamespace("BiocManager", quietly = TRUE))
 # install.packages("BiocManager")

#BiocManager::install("DESeq2")


library(tidyverse)
#install.packages("matrixStats")
library(matrixStats)
library(DESeq2)
library(phyloseq)
library(csv)

## Make a new phyloseq object with the unrarefied table
#metadata
sdata <- as.csv("/Users/lgschaer/Desktop/Permafrost_180/SIP_16S/metadata_v2.csv", row.names = 1, header = TRUE, sep = ",", check.names = TRUE, stringsAsFactors = TRUE) %>%
  mutate(Fraction = factor(Fraction, levels = c("H", "MH", "M", "ML", "L")),
         Day = as.character(Day),
         Temperature = as.character(Temperature))
head(sdata)

#sequence table
sequence_table <- read_table("/Users/lgschaer/Desktop/Permafrost_180/SIP_16S/qiime_out/unrare_table/table-dn-99.txt", skip = 1, col_names = TRUE) %>% 
  column_to_rownames(var = "SampleID") %>%
  t() 
sequence_table[1:5,1:5] 

#taxa table
taxa <- read_table("/Users/lgschaer/Desktop/Permafrost_180/SIP_16S/qiime_out/taxonomy.tsv") %>%
  select(Feature, ID) %>%
  separate(ID, into = c(NA, "Kingdom", NA, "Phylum", NA, "Class", NA, "Order", NA, "Family", NA, "Genus", NA, "Species"), sep = ";|__") %>%
  column_to_rownames(var = "Feature") %>%
  as.matrix()
head(taxa)
#View(taxa)

#format for phyloseq
samdata = sample_data(sdata)
seqtab = otu_table(sequence_table, taxa_are_rows = FALSE)
taxtab = tax_table(taxa)

ps = phyloseq(otu_table(seqtab), tax_table(taxtab), sample_data(samdata))
ps

write_rds(ps, "/Users/lgschaer/Desktop/Permafrost_180/SIP_16S/phyloseq_norare_05162024.rds")

#Heaviest 18O vs Lightest 18O

subPs <- ps %>%
  subset_samples(Fraction == "H" | Fraction == "L") %>%
  subset_samples(Day == "184") %>%
  subset_samples(Isotope == "18O") %>%
  subset_samples(Temperature == "-1.5")
subPs

subPs_test <- subPs %>%
  subset_taxa(Genus == "Methanobacterium_B")
sample_sums(subPs_test)

subPs_dsq <- phyloseq_to_deseq2(subPs, ~ Fraction)

subPs_dsq$Fraction<-factor(subPs_dsq$Fraction, levels = c("L", "H"))

dds_subPs_dsq <- DESeq(subPs_dsq)

resultsNames(dds_subPs_dsq) 

res <- as.data.frame(results(dds_subPs_dsq, contrast=c("Fraction","H","L"))) %>% mutate(Comparison="Lightest vs Heaviest") %>% rownames_to_column(var = "taxon")
res_H_vs_L <- as.data.frame(tax_table(subPs)) %>% rownames_to_column(var = "taxon") %>% full_join(res) %>% 
  filter(!is.na(padj)) %>%
  mutate(
    threshold = ifelse(padj <= 0.001 & abs(log2FoldChange) >= 2, "Enriched", "Not_Enriched"),
    Enriched_Genus = ifelse(threshold == "Enriched", as.character(Genus), "Not Enriched"),
    Enriched_Genus = ifelse(is.na(Enriched_Genus), "Unclassified Genus", Enriched_Genus),
    Enriched_Genus = factor(Enriched_Genus, levels = c("Not Enriched", "Unclassified Genus", "LD21", "Methanobacterium_B", "Fen-1298", "DUES01",            
                                                       "Smithella", "Terracidiphilus", "SURF-40", "JAJYIG01"))
  ) 
head(res_H_vs_L)
#View(res_H_vs_L)

# Are any ASVs enriched?
any(res_H_vs_L$threshold == "Enriched")
sum(res_H_vs_L$threshold == "Enriched")

colFunc <- colorRampPalette(c("purple4", "lightpink", "firebrick", "orange", "lemonchiffon", "olivedrab4", "darkcyan", "lightblue", "darkblue"))
color_list <- colFunc(length(unique(res_H_vs_L$Enriched_Genus)))
genus_colors <- c("white", "black", color_list)
length(genus_colors)

ggplot(data=res_H_vs_L, aes(x=log2FoldChange, y=-log10(padj))) +
  geom_point(aes(fill=Enriched_Genus), shape = 21, color = "black", size=6) +
  facet_grid(cols = vars(Comparison))+
  scale_fill_manual(values=genus_colors) +
  #scale_shape_manual(values=shapes) +
  labs(x = "log2 fold change", 
       y = "-log10 p-value") +
  theme_classic(base_size = 14)+
  geom_hline(yintercept = -log10(0.001), colour="#990000", linetype="dashed") + 
  geom_vline(xintercept = 2, colour="#990000", linetype="dashed") + 
  geom_vline(xintercept = -2, colour="#990000", linetype="dashed")+ 
  theme(axis.text.y.left = element_text(size = 20),
        axis.text.x = element_text(size = 20, angle = 0, vjust = 0.5, hjust = 0.5),
        axis.title.y = element_text(size = 25),
        axis.title.x = element_text(size = 25),
        legend.text = element_text(size = 18),
        legend.title = element_blank(),
        strip.text.x = element_text(size = 25, face = "bold", angle = 0),
        legend.position = "bottom",
        title = element_text(size = 18))



## Heaviest 18O vs Heaviest 16O

subPs <- ps %>%
  subset_samples(Fraction == "H") %>%
  subset_samples(Day == "184") %>%
  subset_samples(Isotope == "18O" | Isotope == "16O") %>%
  subset_samples(Temperature == "-1.5")
subPs

subPs_test <- subPs %>%
  subset_taxa(Genus == "Paludibacter")
sample_sums(subPs_test)

#add pseudo count of 1
#subPs@otu_table <- as.matrix(subPs@otu_table)+1

subPs_dsq <- phyloseq_to_deseq2(subPs, ~ Isotope)

#subPs_dsq$Fraction<-factor(subPs_dsq$Fraction, levels = c("L", "H"))

dds_subPs_dsq <- DESeq(subPs_dsq)

resultsNames(dds_subPs_dsq) 

res <- as.data.frame(results(dds_subPs_dsq, contrast=c("Isotope","18O","16O"))) %>% mutate(Comparison="16O vs 18O") %>% rownames_to_column(var = "taxon")
res_18O_vs_16O <- as.data.frame(tax_table(subPs)) %>% rownames_to_column(var = "taxon") %>% full_join(res) %>% 
  filter(!is.na(padj)) %>%
  mutate(
    threshold = ifelse(padj <= 0.001 & abs(log2FoldChange) >= 2, "Enriched", "Not_Enriched"),
    Enriched_Genus = ifelse(threshold == "Enriched", as.character(Genus), "Not Enriched"),
    Enriched_Genus = ifelse(is.na(Enriched_Genus), "Unclassified Genus", Enriched_Genus),
    Enriched_Genus = factor(Enriched_Genus, levels = c("Not Enriched", "Unclassified Genus", "Terracidiphilus", "Paludibacter", "Sediminibacterium"))
  ) 
head(res_18O_vs_16O)
#View(res_18O_vs_16O)

unique(res_18O_vs_16O$Enriched_Genus)

# Are any ASVs enriched?
any(res_18O_vs_16O$threshold == "Enriched")
sum(res_18O_vs_16O$threshold == "Enriched")

colFunc <- colorRampPalette(c("purple4", "lightpink", "firebrick", "orange", "lemonchiffon", "olivedrab4", "darkcyan", "lightblue", "darkblue"))
color_list <- colFunc(length(unique(res_H_vs_L$Enriched_Genus)))
genus_colors <- c("white", "black", color_list)
length(genus_colors)

ggplot(data=res_18O_vs_16O, aes(x=log2FoldChange, y=-log10(padj))) +
  geom_point(aes(fill=Enriched_Genus), shape = 21, color = "black", size=6) +
  facet_grid(cols = vars(Comparison))+
  scale_fill_manual(values=genus_colors) +
  #scale_shape_manual(values=shapes) +
  labs(x = "log2 fold change", 
       y = "-log10 p-value") +
  theme_classic(base_size = 14)+
  geom_hline(yintercept = -log10(0.001), colour="#990000", linetype="dashed") + 
  geom_vline(xintercept = 2, colour="#990000", linetype="dashed") + 
  geom_vline(xintercept = -2, colour="#990000", linetype="dashed")+ 
  theme(axis.text.y.left = element_text(size = 20),
        axis.text.x = element_text(size = 20, angle = 0, vjust = 0.5, hjust = 0.5),
        axis.title.y = element_text(size = 25),
        axis.title.x = element_text(size = 25),
        legend.text = element_text(size = 18),
        legend.title = element_blank(),
        strip.text.x = element_text(size = 25, face = "bold", angle = 0),
        legend.position = "bottom",
        title = element_text(size = 18))


