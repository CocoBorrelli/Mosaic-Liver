## ANALYSIS PIPELINE

#1. demultiplexing of FASTQS with Bcl2Fastq (Illumina), each sample has a unique barcode in the P7 primer. 

#---------------------------––--#

#2. trimming with cutadapt > cutadapt loop
cutadapt -g CACCG -o B11_trimup.fastq ../demux_fastqs/B11_merged.fastq.gz
cutadapt -a GTTTT -o B11_trimmed.fastq B11_trimup.fastq

#---------------------------––--#

#3. bowtie alignment
awk -F ',' '{print ">"$1"\n"$2}' library2.csv > library2.fa #build bowtie index 
bowtie2-build library2.fa bowtie2_ind_library
bowtie2 -x bowtie2_ind_library -U ../cutadapt_output/B11_trimmed.fastq --norc | samtools view -bS - > B11.bam

#---------------------------––--#

#4. check library retention with correlation between plasmid prep and post injection library (fig 1)
mageck count -l library2.csv -n premets_SPH1 --sample-label "premets,SPH1" --fastq 3071.bam SPH1.bam --norm-method total

premets_SPH1.count_normalized <- read.delim("/media/Coco/MOSAIC LIVER/Experiments/Screen1/Fastqs/MAGECK/Bam_files/premets_SPH1.count_normalized.txt")
View(premets_SPH1.count_normalized)

library("ggpubr")
ggscatter(premets_SPH1.count_normalized, x = "preinj", y = "plasmid", 
          add = "reg.line", conf.int = TRUE, conf.int.level = 0.95, color = "black",
          fill = "#83C441",
          cor.coef = TRUE, cor.method = "pearson",
          xlab = "sgRNA count hepatocytes", ylab = "sgRNA count plasmid")+
          ggsave("plasmid_preinj.pdf", width = 5, height = 5)

#---------------------------––--#

#5. check correlation of library batches (3 independent plasmid preps)
cd /NAS/Coco/MOSAIC\ LIVER/Experiments/Screen1/Fastqs/MAGECK/Bam_files
mageck count -l library.csv -n SPH1 --sample-label "SPH1" --fastq A3.bam --norm-method total
mageck count -l library.csv -n SPH2 --sample-label "SPH2" --fastq SPH2.bam --norm-method total
mageck count -l library.csv -n SPH3 --sample-label "SPH3" --fastq SPH3.bam --norm-method total

#merge all counts
`SPH1.count` <- read.delim("/media/Coco/MOSAIC LIVER/Experiments/Screen1/Fastqs/MAGECK/Bam_files/SPH1.count.txt")
View(SPH1.count)
`SPH2.count` <- read.delim("/media/Coco/MOSAIC LIVER/Experiments/Screen1/Fastqs/MAGECK/Bam_files/SPH2.count.txt")
View(SPH2.count)
`SPH3.count` <- read.delim("/media/Coco/MOSAIC LIVER/Experiments/Screen1/Fastqs/MAGECK/Bam_files/SPH3.count.txt")
View(SPH3.count)

SPH1_2_3 <- list(`SPH1.count` ,`SPH2.count`, `SPH3.count`) %>% 
  lapply(tibble::rownames_to_column) %>% purrr::reduce(full_join, by="sgRNA") %>% 
  mutate_all(~replace(., is.na(.), 0))
colnames(SPH1_2_3)
head(SPH1_2_3)
SPH1_2_3_counts <- SPH1_2_3[, c(2,3,4,7,10)]
colnames(SPH1_2_3_counts) <- c( "sgRNA" ,"Gene", "SPH1", "SPH2", "SPH3")

View(SPH1_2_3_counts)

#check correlation
library("ggpubr")
a <- ggscatter(SPH1_2_3_counts, x = "SPH1", y = "SPH2", 
               add = "reg.line", conf.int = TRUE, conf.int.level = 0.95, color = "black",
               fill = "#83C441",
               cor.coef = TRUE, cor.method = "pearson",
               xlab = "HI1_counts", ylab = "HI2_counts")#+ggsave("plasmid_preinj.pdf", width = 5, height = 5)

b <- ggscatter(SPH1_2_3_counts, x = "SPH1", y = "SPH3", 
               add = "reg.line", conf.int = TRUE, conf.int.level = 0.95, color = "black",
               fill = "#83C441",
               cor.coef = TRUE, cor.method = "pearson",
               xlab = "HI1_counts", ylab = "HI3_counts")#+ggsave("plasmid_preinj.pdf", width = 5, height = 5)

c <- ggscatter(SPH1_2_3_counts, x = "SPH2", y = "SPH3", 
               add = "reg.line", conf.int = TRUE, conf.int.level = 0.95, color = "black",
               fill = "#83C441",
               cor.coef = TRUE, cor.method = "pearson",
               xlab = "HI2_counts", ylab = "HI3_counts")#+ggsave("plasmid_preinj.pdf", width = 5, height = 5)
ggarrange(a,b,c,ncol = 3, nrow = 1) +ggsave("/media/Coco/MOSAIC LIVER/Manuscript/Graphs/plasmid_SPH123.pdf", width = 10, height = 3)

#---------------------------––--#

#6. coverage plot sup fig 3
Coverageplot <- read.csv("/media/Coco/MOSAIC LIVER/Experiments/Screen1/Fastqs/Coverageplot.csv")
View(Coverageplot)

ggplot(data=Coverageplot, aes(x=Sample, y=Total, fill=library)) +
  geom_bar(stat="identity", color="white", position=position_dodge())+
  theme_classic() + scale_fill_manual(values=c('#81C341','#818641', "#2F8641"))+
  scale_x_discrete(guide = guide_axis(angle = 45))  +ggsave("/media/Coco/MOSAIC LIVER/Manuscript/Graphs/coverage.pdf", width = 4, height = 3)
 
#---------------------------––--#

#7. check of all mice individually 

####SPH1####
mageck count -l library2.csv -n 3068 --sample-label "distal,proximal" --fastq 3068distal.bam 3068proximal.bam --norm-method total
mageck count -l library2.csv -n 3174 --sample-label "distal,proximal" --fastq 3174distal.bam 3174proximal.bam --norm-method total

#nocre
mageck count -l library2.csv -n 3070 --sample-label "distal,proximal" --fastq 3070distal.bam 3070proximal.bam --norm-method total

mageck test -k 3068.count_normalized.txt -t proximal -c distal -n 3068
`3068.gene_summary` <- read.delim("/media/Coco/MOSAIC LIVER/Experiments/Screen1/Fastqs/MAGECK/Bam_files/3068.gene_summary.txt")
View(`3068.gene_summary`)
#negatives: Yap, Vegfa
#positives: Plxnb2, Gprc5c/b, Psen1

mageck test -k 3174.count_normalized.txt -t proximal -c distal -n 3174
`3174.gene_summary` <- read.delim("/media/Coco/MOSAIC LIVER/Experiments/Screen1/Fastqs/MAGECK/Bam_files/3174.gene_summary.txt")
View(`3174.gene_summary`)
#negatives:Gpr5c5/b 

mageck test -k 3070.count_normalized.txt -t proximal -c distal -n 3070
`3070.gene_summary` <- read.delim("/media/Coco/MOSAIC LIVER/Experiments/Screen1/Fastqs/MAGECK/Bam_files/3070.gene_summary.txt")
View(`3070.gene_summary`)

####paired analysis of cre mice####
`3068.count_normalized` <- read.delim("/media/Coco/MOSAIC LIVER/Experiments/Screen1/Fastqs/MAGECK/Bam_files/3068.count_normalized.txt")
`3174.count_normalized` <- read.delim("/media/Coco/MOSAIC LIVER/Experiments/Screen1/Fastqs/MAGECK/Bam_files/3174.count_normalized.txt")
`3070.count_normalized` <- read.delim("/media/Coco/MOSAIC LIVER/Experiments/Screen1/Fastqs/MAGECK/Bam_files/3070.count_normalized.txt") 

SPH1 <- list(`3174.count_normalized` ,`3068.count_normalized`) %>% 
  lapply(tibble::rownames_to_column) %>% purrr::reduce(full_join, by="sgRNA") %>% 
  # NB: below we are IMPUTING the genes where they miss in some samples, with value 0:
  # many ways possible, I choose one simple enough but also ok for big data sets:
  # https://stackoverflow.com/questions/8161836/how-do-i-replace-na-values-with-zeros-in-an-r-dataframe
  mutate_all(~replace(., is.na(.), 0))
head(SPH1)
SPH1_counts <- SPH1[, c(2,3,4,5,8,9)]
colnames(SPH1_counts) <- c( "sgRNA" ,"Gene", "3174d","3174p","3068d","3068p")
head(SPH1_counts)
write.table(SPH1_counts, "/media/Coco/MOSAIC LIVER/Experiments/Screen1/Fastqs/MAGECK/Bam_files/SPH1_counts.txt",append = FALSE, sep= "\t",  row.names = FALSE, quote = FALSE, col.names = TRUE)

mageck test -k SPH1_counts.txt -t 3174p,3068p -c 3174d,3068d -n SPH1_paired --paired
SPH1_paired.gene_summary <- read.delim("/media/Coco/MOSAIC LIVER/Experiments/Screen1/Fastqs/MAGECK/Bam_files/SPH1_paired.gene_summary.txt")
View(SPH1_paired.gene_summary)

ggplot(data=SPH1_paired.gene_summary, aes(x=neg.lfc, y=-log10(neg.p.value), label=id)) + 
  geom_vline(xintercept=c(-0.5, 0.5), col="red") + geom_text()+
  geom_hline(yintercept=-log10(0.05), col="red")+
  geom_point() + 
  theme_classic()
ggplot(data=SPH1_paired.gene_summary, aes(x=pos.lfc, y=-log10(pos.p.value), label=id)) + 
  geom_vline(xintercept=c(-0.5, 0.5), col="red") + geom_text()+
  geom_hline(yintercept=-log10(0.05), col="red")+
  geom_point() + 
  theme_classic()

#barplot of pvalue
#has both negative and positive pvalue > take positive p val for positive logfc and neg for neg
neg<-which(SPH1_paired.gene_summary$neg.lfc<0)
pos<-which(SPH1_paired.gene_summary$pos.lfc>0)
colnames(SPH1_paired.gene_summary)
neg_data <- SPH1_paired.gene_summary[neg, c(1,4, 8)]
head(neg_data)
pos_data <- SPH1_paired.gene_summary[pos, c(1,10, 14)]
head(pos_data)
colnames(pos_data) <-  c("id" ,"p.value" ,"lfc"   )
colnames(neg_data) <-  c("id" ,"p.value" ,"lfc"   )
pos_data$log10pval <- -log10(pos_data$p.value)
neg_data$log10pval <- -log10(neg_data$p.value)
neg_data$log10pval <- neg_data$log10pval*-1
plot_data <-rbind(neg_data, pos_data)
#View(plot_data)
plot_data$group <- ifelse(plot_data$lfc < 0, "neg", "pos")
reduced_plot_data <- plot_data %>%
  group_by(group) %>%
  top_n(n = 15, wt = abs(log10pval))
reduced_plot_data

ggplot(data=reduced_plot_data, aes(x=reorder(id, log10pval),y= log10pval,  fill=group)) +
  geom_bar(stat="identity")+
  theme_classic() + scale_fill_manual(values=c('#81C341','#D12026'))+ coord_flip()#+
 # ggsave("/media/Coco/MOSAIC LIVER/Manuscript/Graphs/SPH1.pdf", width = 5, height = 4)

gdata1 = ReadRRA(SPH1_paired.gene_summary)
View(gdata1)
sdata1 = ReadsgRRA(read.delim("/media/Coco/MOSAIC LIVER/Experiments/Screen1/Fastqs/MAGECK/Bam_files/SPH1_paired.sgrna_summary.txt"))
sgRankView(sdata1, top = 2, gene = c("Psen1", "Plxnb2", "Gpc4", "App" )) +ggsave("/media/Coco/MOSAIC LIVER/Manuscript/Graphs/sgRNAplotSPH1.pdf", width = 5, height = 4)
gdata1$LogFDR = -log10(gdata1$FDR)
ScatterView(gdata1, x = "Score", y = "LogFDR", label = "id", 
            model = "volcano", top = 20, y_cut=0.2, x_cut = 0.05, max)
VolcanoView(gdata1, x = "Score", y = "FDR", Label = "id", top = 10, x_cut = 0.05, y_cut = 0.8, alpha=1, max.overlaps=Inf)+ 
  theme_classic() + scale_fill_manual(values=c('#81C341', "grey", '#D12026')) +#+ ylim(0,3)+ xlim(-2,2) #+
  ggsave("/media/Coco/MOSAIC LIVER/Manuscript/Graphs/VolcanoSPH1_paired.pdf", width = 5, height = 4)

####SPH2####
####check mice individually####
mageck count -l library2.csv -n 3039 --sample-label "distal,proximal" --fastq 3039distal.bam 3039proximal.bam --norm-method total
mageck count -l library2.csv -n 3379 --sample-label "distal,proximal" --fastq 3379distal.bam 3379proximal.bam --norm-method total
mageck count -l library2.csv -n 3425 --sample-label "distal,proximal" --fastq 3425distal.bam 3425proximal.bam --norm-method total

#nocre
mageck count -l library2.csv -n 3378 --sample-label "distal,proximal" --fastq 3378distal.bam 3378proximal.bam --norm-method total
mageck count -l library2.csv -n 3381 --sample-label "distal,proximal" --fastq 3381distal.bam 3381proximal.bam --norm-method total

mageck test -k 3039.count_normalized.txt -t proximal -c distal -n 3039
`3039.gene_summary` <- read.delim("/media/Coco/MOSAIC LIVER/Experiments/Screen1/Fastqs/MAGECK/Bam_files/3039.gene_summary.txt")
View(`3039.gene_summary`)
#negatives: Yap, Vegfa
#positives: Plxnb2, Gprc5c/b, Psen1

mageck test -k 3379.count_normalized.txt -t proximal -c distal -n 3379
`3379.gene_summary` <- read.delim("/media/Coco/MOSAIC LIVER/Experiments/Screen1/Fastqs/MAGECK/Bam_files/3379.gene_summary.txt")
View(`3379.gene_summary`)
#ok

mageck test -k 3378.count_normalized.txt -t proximal -c distal -n 3378
`3378.gene_summary` <- read.delim("/media/Coco/MOSAIC LIVER/Experiments/Screen1/Fastqs/MAGECK/Bam_files/3378.gene_summary.txt")
View(`3378.gene_summary`)

#nocres are weird
mageck test -k 3381.count_normalized.txt -t proximal -c distal -n 3381
`3381.gene_summary` <- read.delim("/media/Coco/MOSAIC LIVER/Experiments/Screen1/Fastqs/MAGECK/Bam_files/3381.gene_summary.txt")
View(`3381.gene_summary`)

mageck test -k 3381.count_normalized.txt -t proximal -c distal -n 3381
`3381.gene_summary` <- read.delim("/media/Coco/MOSAIC LIVER/Experiments/Screen1/Fastqs/MAGECK/Bam_files/3381.gene_summary.txt")
View(`3381.gene_summary`)

#####paired analysis of cre mice#####
`3039.count_normalized` <- read.delim("/media/Coco/MOSAIC LIVER/Experiments/Screen1/Fastqs/MAGECK/Bam_files/3039.count_normalized.txt")
`3379.count_normalized` <- read.delim("/media/Coco/MOSAIC LIVER/Experiments/Screen1/Fastqs/MAGECK/Bam_files/3379.count_normalized.txt")
`3425.count_normalized` <- read.delim("/media/Coco/MOSAIC LIVER/Experiments/Screen1/Fastqs/MAGECK/Bam_files/3425.count_normalized.txt") 

SPH2 <- list(`3039.count_normalized` ,`3379.count_normalized`, `3425.count_normalized`) %>% 
  lapply(tibble::rownames_to_column) %>% purrr::reduce(full_join, by="sgRNA") %>% 
  # NB: below we are IMPUTING the genes where they miss in some samples, with value 0:
  # many ways possible, I choose one simple enough but also ok for big data sets:
  # https://stackoverflow.com/questions/8161836/how-do-i-replace-na-values-with-zeros-in-an-r-dataframe
  mutate_all(~replace(., is.na(.), 0))
head(SPH2)
SPH2_counts <- SPH2[, c(2,3,4,5,8,9,12,13)]
colnames(SPH2_counts) <- c( "sgRNA" ,"Gene", "3039d","3039p","3379d","3379p", "3425d", "3425p" )
head(SPH2_counts)
write.table(SPH2_counts, "/media/Coco/MOSAIC LIVER/Experiments/Screen1/Fastqs/MAGECK/Bam_files/SPH2_counts.txt",append = FALSE, sep= "\t",  row.names = FALSE, quote = FALSE, col.names = TRUE)

mageck test -k SPH2_counts.txt -t 3039p,3379p,3425p -c 3039d,3379d,3425d -n SPH2_paired --paired
SPH2_paired.gene_summary <- read.delim("/media/Coco/MOSAIC LIVER/Experiments/Screen1/Fastqs/MAGECK/Bam_files/SPH2_paired.gene_summary.txt")
View(SPH2_paired.gene_summary)

ggplot(data=SPH2_paired.gene_summary, aes(x=neg.lfc, y=-log10(neg.p.value), label=id)) + 
  geom_vline(xintercept=c(-0.5, 0.5), col="red") + geom_text()+
  geom_hline(yintercept=-log10(0.05), col="red")+
  geom_point() + 
  theme_classic()
ggplot(data=SPH2_paired.gene_summary, aes(x=pos.lfc, y=-log10(pos.p.value), label=id)) + 
  geom_vline(xintercept=c(-0.5, 0.5), col="red") + geom_text()+
  geom_hline(yintercept=-log10(0.05), col="red")+
  geom_point() + 
  theme_classic()

#barplot of pvalue
#has both negative and positive pvalue > take positive p val for positive logfc and neg for neg
neg<-which(SPH2_paired.gene_summary$neg.lfc<0)
pos<-which(SPH2_paired.gene_summary$pos.lfc>0)
colnames(SPH2_paired.gene_summary)
neg_data <- SPH2_paired.gene_summary[neg, c(1,4, 8)]
head(neg_data)
pos_data <- SPH2_paired.gene_summary[pos, c(1,10, 14)]
head(pos_data)
colnames(pos_data) <-  c("id" ,"p.value" ,"lfc"   )
colnames(neg_data) <-  c("id" ,"p.value" ,"lfc"   )
pos_data$log10pval <- -log10(pos_data$p.value)
neg_data$log10pval <- -log10(neg_data$p.value)
neg_data$log10pval <- neg_data$log10pval*-1
plot_data <-rbind(neg_data, pos_data)
#View(plot_data)
plot_data$group <- ifelse(plot_data$lfc < 0, "neg", "pos")
reduced_plot_data <- plot_data %>%
  group_by(group) %>%
  top_n(n = 15, wt = abs(log10pval))
reduced_plot_data

ggplot(data=reduced_plot_data, aes(x=reorder(id, log10pval),y= log10pval,  fill=group)) +
  geom_bar(stat="identity")+
  theme_classic() + scale_fill_manual(values=c('#81C341','#D12026'))+ coord_flip()
 # ggsave("/media/Coco/MOSAIC LIVER/Manuscript/Graphs/SPH2.pdf", width = 5, height = 4)

gdata2 = ReadRRA(SPH2_paired.gene_summary)
View(gdata2)
sdata2 = ReadsgRRA(read.delim("/media/Coco/MOSAIC LIVER/Experiments/Screen1/Fastqs/MAGECK/Bam_files/SPH2_paired.sgrna_summary.txt"))
sgRankView(sdata2, top = 2, gene = c("Psen1", "Plxnb2", "Gpc4", "App" )) +ggsave("/media/Coco/MOSAIC LIVER/Manuscript/Graphs/sgRNAplot_SPH2.pdf", width = 5, height = 4)
gdata2$LogFDR = -log10(gdata2$FDR)
ScatterView(gdata2, x = "Score", y = "LogFDR", label = "id", 
            model = "volcano", top = 20, y_cut=0.2, x_cut = 0.05, max)
VolcanoView(gdata2, x = "Score", y = "FDR", Label = "id", top = 10, x_cut = 0.05, y_cut = 0.8, alpha=1, max.overlaps=Inf)+ 
  theme_classic() + scale_fill_manual(values=c('#81C341', "grey", '#D12026')) +#+ ylim(0,3)+ xlim(-2,2) #+
  ggsave("/media/Coco/MOSAIC LIVER/Manuscript/Graphs/VolcanoSPH2_paired.pdf", width = 5, height = 4)

####SPH3####
#####check mice individually####
mageck count -l library2.csv -n 3434 --sample-label "distal,proximal" --fastq 3434distal.bam 3434proximal.bam --norm-method total
mageck count -l library2.csv -n 3436 --sample-label "distal,proximal" --fastq 3436distal.bam 3436proximal.bam --norm-method total

#nocre
mageck count -l library2.csv -n 3428 --sample-label "distal,proximal" --fastq 3428distal.bam 3428proximal.bam --norm-method total
mageck count -l library2.csv -n 3431 --sample-label "distal,proximal" --fastq 3431distal.bam 3431proximal.bam --norm-method total

mageck test -k 3434.count_normalized.txt -t proximal -c distal -n 3434
`3434.gene_summary` <- read.delim("/media/Coco/MOSAIC LIVER/Experiments/Screen1/Fastqs/MAGECK/Bam_files/3434.gene_summary.txt")
View(`3434.gene_summary`)
#also weird

mageck test -k 3436.count_normalized.txt -t proximal -c distal -n 3436
`3436.gene_summary` <- read.delim("/media/Coco/MOSAIC LIVER/Experiments/Screen1/Fastqs/MAGECK/Bam_files/3436.gene_summary.txt")
View(`3436.gene_summary`)
#switched!! recount
mageck count -l library2.csv -n 3436 --sample-label "distal,proximal" --fastq 3436proximal.bam 3436distal.bam --norm-method total
mageck test -k 3436.count_normalized.txt -t proximal -c distal -n 3436
`3436.gene_summary` <- read.delim("/media/Coco/MOSAIC LIVER/Experiments/Screen1/Fastqs/MAGECK/Bam_files/3436.gene_summary.txt")
View(`3436.gene_summary`)

#no cres are weird here too
mageck test -k 3428.count_normalized.txt -t proximal -c distal -n 3428
`3428.gene_summary` <- read.delim("/media/Coco/MOSAIC LIVER/Experiments/Screen1/Fastqs/MAGECK/Bam_files/3428.gene_summary.txt")
View(`3428.gene_summary`)

mageck test -k 3431.count_normalized.txt -t proximal -c distal -n 3431
`3431.gene_summary` <- read.delim("/media/Coco/MOSAIC LIVER/Experiments/Screen1/Fastqs/MAGECK/Bam_files/3431.gene_summary.txt")
View(`3431.gene_summary`)


#####paired analysis of cre mice#####
`3434.count_normalized` <- read.delim("/media/Coco/MOSAIC LIVER/Experiments/Screen1/Fastqs/MAGECK/Bam_files/3434.count_normalized.txt")
`3436.count_normalized` <- read.delim("/media/Coco/MOSAIC LIVER/Experiments/Screen1/Fastqs/MAGECK/Bam_files/3436.count_normalized.txt")

SPH3 <- list(`3434.count_normalized`,`3436.count_normalized`) %>% 
  lapply(tibble::rownames_to_column) %>% purrr::reduce(full_join, by="sgRNA") %>% 
  # NB: below we are IMPUTING the genes where they miss in some samples, with value 0:
  # many ways possible, I choose one simple enough but also ok for big data sets:
  # https://stackoverflow.com/questions/8161836/how-do-i-replace-na-values-with-zeros-in-an-r-dataframe
  mutate_all(~replace(., is.na(.), 0))
head(SPH3)
SPH3_counts <- SPH3[, c(2,3,4,5,8,9)]
colnames(SPH3_counts) <- c( "sgRNA" ,"Gene", "3434d","3434p","3436d","3436p")
head(SPH3_counts)
write.table(SPH3_counts, "/media/Coco/MOSAIC LIVER/Experiments/Screen1/Fastqs/MAGECK/Bam_files/SPH3_counts.txt",append = FALSE, sep= "\t",  row.names = FALSE, quote = FALSE, col.names = TRUE)

mageck test -k SPH3_counts.txt -t 3434p,3436p -c 3434d,3436d -n SPH3_paired --paired
SPH3_paired.gene_summary <- read.delim("/media/Coco/MOSAIC LIVER/Experiments/Screen1/Fastqs/MAGECK/Bam_files/SPH3_paired.gene_summary.txt")
head(SPH3_paired.gene_summary)
View(SPH3_paired.gene_summary)

ggplot(data=SPH3_paired.gene_summary, aes(x=neg.lfc, y=-log10(neg.p.value), label=id)) + 
  geom_vline(xintercept=c(-0.5, 0.5), col="red") + geom_text()+
  geom_hline(yintercept=-log10(0.05), col="red")+
  geom_point() + 
  theme_classic()
ggplot(data=SPH3_paired.gene_summary, aes(x=pos.lfc, y=-log10(pos.p.value), label=id)) + 
  geom_vline(xintercept=c(-0.5, 0.5), col="red") + geom_text()+
  geom_hline(yintercept=-log10(0.05), col="red")+
  geom_point() + 
  theme_classic()

#barplot of pvalue
#has both negative and positive pvalue > take positive p val for positive logfc and neg for neg
neg<-which(SPH3_paired.gene_summary$neg.lfc<0)
pos<-which(SPH3_paired.gene_summary$pos.lfc>0)
colnames(SPH3_paired.gene_summary)
neg_data <- SPH3_paired.gene_summary[neg, c(1,4, 8)]
head(neg_data)
pos_data <- SPH3_paired.gene_summary[pos, c(1,10, 14)]
head(pos_data)
colnames(pos_data) <-  c("id" ,"p.value" ,"lfc"   )
colnames(neg_data) <-  c("id" ,"p.value" ,"lfc"   )
pos_data$log10pval <- -log10(pos_data$p.value)
neg_data$log10pval <- -log10(neg_data$p.value)
neg_data$log10pval <- neg_data$log10pval*-1
plot_data <-rbind(neg_data, pos_data)
#View(plot_data)
plot_data$group <- ifelse(plot_data$lfc < 0, "neg", "pos")
reduced_plot_data <- plot_data %>%
  group_by(group) %>%
  top_n(n = 20, wt = abs(log10pval))
reduced_plot_data

ggplot(data=reduced_plot_data, aes(x=reorder(id, log10pval),y= log10pval,  fill=group)) +
  geom_bar(stat="identity")+
  theme_classic() + scale_fill_manual(values=c('#81C341','#D12026'))+ coord_flip()#+
  #ggsave("/media/Coco/MOSAIC LIVER/Manuscript/Graphs/SPH2.pdf", width = 5, height = 4)

gdata3 = ReadRRA(SPH3_paired.gene_summary)
View(gdata3)
sdata3 = ReadsgRRA(read.delim("/media/Coco/MOSAIC LIVER/Experiments/Screen1/Fastqs/MAGECK/Bam_files/SPH3_paired.sgrna_summary.txt"))
sgRankView(sdata3, top = 2, gene = c("Psen1", "Plxnb2", "Gpc4", "App" )) +ggsave("/media/Coco/MOSAIC LIVER/Manuscript/Graphs/sgRNAplotSPH3.pdf", width = 5, height = 4)
gdata3$LogFDR = -log10(gdata3$FDR)
ScatterView(gdata3, x = "Score", y = "LogFDR", label = "id", 
            model = "volcano", top = 10, y_cut=0.25, x_cut = 0.05)
VolcanoView(gdata3, x = "Score", y = "FDR", Label = "id", x_cut = 0.05, y_cut = 0.965, alpha=1)+
  theme_classic() + scale_fill_manual(values=c('#81C341', "grey", '#D12026')) +#+ ylim(0,3)+ xlim(-2,2) #+
  ggsave("/media/Coco/MOSAIC LIVER/Manuscript/Graphs/VolcanoSPH3_paired.pdf", width = 5, height = 4)


#########PAIRED ANALYSIS OF ALL###########
SPH1_2_3 <- list(`3174.count_normalized` ,`3068.count_normalized`, `3039.count_normalized`, 
                 `3379.count_normalized`, `3425.count_normalized`, `3434.count_normalized`, 
                 `3436.count_normalized`) %>% 
  lapply(tibble::rownames_to_column) %>% purrr::reduce(full_join, by="sgRNA") %>% 
  # NB: below we are IMPUTING the genes where they miss in some samples, with value 0:
  # many ways possible, I choose one simple enough but also ok for big data sets:
  # https://stackoverflow.com/questions/8161836/how-do-i-replace-na-values-with-zeros-in-an-r-dataframe
  mutate_all(~replace(., is.na(.), 0))
colnames(SPH1_2_3)
SPH1_2_3_counts <- SPH1_2_3[, c(2,3, 4,5,8,9,12,13,16,17,20,21,24,25,28,29)]
colnames(SPH1_2_3_counts) <- c( "sgRNA" ,"Gene", "3174d","3174p","3068d","3068p","3039d","3039p",
                                "3379d", "3379p","3425d", "3425p", "3434d", "3434p", "3436d", "3436p")

head(SPH1_2_3_counts)
write.table(SPH1_2_3_counts, "/media/Coco/MOSAIC LIVER/Experiments/Screen1/Fastqs/MAGECK/Bam_files/SPH1_2_3_counts.txt",append = FALSE, sep= "\t",  row.names = FALSE, quote = FALSE, col.names = TRUE)
mageck test -k SPH1_2_3_counts.txt -t 3174p,3068p,3039p,3379p,3425p,3434p,3436p -c 3174d,3068d,3039d,3379d,3425d,3434d,3436d -n paired_proximal_distal --paired

#mageck test -k SPH1_2_3_counts.txt -t 3174p,3068p,3039p,3379p,3425p,3434p,3436p -c 3174d,3068d,3039d,3379d,3425d,3434d,3436d -n paired_proximal_distal --paired --gene-lfc-method alphamedian #alphamedian changes a lot
#mageck test -k SPH1_2_3_counts.txt -t 3174p,3068p,3039p,3379p,3425p,3434p,3436p -c 3174d,3068d,3039d,3379d,3425d,3434d,3436d -n paired_proximal_distal --paired --gene-lfc-method alphamedian --remove-zero both
#mageck test -k SPH1_2_3_counts.txt -t 3174p,3068p,3039p,3379p,3425p,3434p,3436p -c 3174d,3068d,3039d,3379d,3425d,3434d,3436d -n paired_proximal_distal --paired --adjust-method pounds --gene-lfc-method alphamedian

paired_proximal_distal.sgrna_summary <- read.delim("/media/Coco/MOSAIC LIVER/Experiments/Screen1/Fastqs/MAGECK/Bam_files/paired_proximal_distal.sgrna_summary.txt")
sdata = ReadsgRRA(paired_proximal_distal.sgrna_summary)
View(sdata)
sgRankView(sdata)
sgRankView(sdata, gene = c("Psen1","Nrp2",  "Plxnb2", "Ltb", "App", "Saa1"))+ggsave("/media/Coco/MOSAIC LIVER/Manuscript/Graphs/sgRNAplotpaired.pdf", width = 5, height = 4)

paired_proximal_distal.gene_summary <- read.delim("/media/Coco/MOSAIC LIVER/Experiments/Screen1/Fastqs/MAGECK/Bam_files/paired_proximal_distal.gene_summary.txt")
View(paired_proximal_distal.gene_summary)
gdata = ReadRRA(paired_proximal_distal.gene_summary)
View(gdata)
gdata$LogFDR = -log10(gdata$FDR)
ScatterView(gdata, x = "Score", y = "LogFDR", label = "id", 
            model = "volcano", top = 10, y_cut=0.25, x_cut = 0.05)
VolcanoView(gdata, x = "Score", y = "FDR", Label = "id", top = 10, x_cut = 0.05, y_cut = 0.8, alpha=1)+
  theme_classic() + scale_fill_manual(values=c('#81C341', "grey", '#D12026'))+ ylim(0,3)+ xlim(-1,1) +
 ggsave("/media/Coco/MOSAIC LIVER/Manuscript/Graphs/VolcanoSPH123_paired_all.pdf", width = 5, height = 4)



ggplot(data=reduced_plot_data, aes(x=reorder(id, Score),y= Score,  fill=group)) +
  geom_bar(stat="identity")+
  theme_classic() + scale_fill_manual(values=c('#81C341','#D12026'))+ coord_flip()#+
#ggsave("/media/Coco/MOSAIC LIVER/Manuscript/Graphs/SPH123.pdf", width = 5, height = 4)


ggplot(data=paired_proximal_distal.gene_summary, aes(x=neg.lfc, y=-log10(neg.p.value), label=id)) + 
  geom_vline(xintercept=c(-0.5, 0.5), col="red") + geom_text()+
  geom_hline(yintercept=-log10(0.05), col="red")+
  geom_point() + 
  theme_classic()
ggplot(data=paired_proximal_distal.gene_summary, aes(x=pos.lfc, y=-log10(pos.p.value), label=id)) + 
  geom_vline(xintercept=c(-0.5, 0.5), col="red") + geom_text()+
  geom_hline(yintercept=-log10(0.05), col="red")+
  geom_point() + 
  theme_classic()

#barplot of pvalue
#has both negative and positive pvalue > take positive p val for positive logfc and neg for neg
neg<-which(paired_proximal_distal.gene_summary$neg.lfc<0)
pos<-which(paired_proximal_distal.gene_summary$pos.lfc>0)
colnames(paired_proximal_distal.gene_summary)
neg_data <- paired_proximal_distal.gene_summary[neg, c(1,4, 8)]
head(neg_data)
pos_data <- paired_proximal_distal.gene_summary[pos, c(1,10, 14)]
head(pos_data)
colnames(pos_data) <-  c("id" ,"p.value" ,"lfc"   )
colnames(neg_data) <-  c("id" ,"p.value" ,"lfc"   )
pos_data$log10pval <- -log10(pos_data$p.value)
neg_data$log10pval <- -log10(neg_data$p.value)
neg_data$log10pval <- neg_data$log10pval*-1
plot_data <-rbind(neg_data, pos_data)
#View(plot_data)
plot_data$group <- ifelse(plot_data$lfc < 0, "neg", "pos")
reduced_plot_data <- plot_data %>%
  group_by(group) %>%
  top_n(n = 20, wt = abs(log10pval))
reduced_plot_data

ggplot(data=reduced_plot_data, aes(x=reorder(id, log10pval),y= log10pval,  fill=group)) +
  geom_bar(stat="identity")+
  theme_classic() + scale_fill_manual(values=c('#81C341','#D12026'))+ coord_flip()#+ggsave("/media/Coco/MOSAIC LIVER/Manuscript/Graphs/SPH2.pdf", width = 5, height = 4)

#barplot by logfc 
neg<-which(paired_proximal_distal.gene_summary$neg.lfc<0)
pos<-which(paired_proximal_distal.gene_summary$pos.lfc>0)
colnames(paired_proximal_distal.gene_summary)
neg_data <- paired_proximal_distal.gene_summary[neg, c(1,8)]
head(neg_data)
pos_data <- paired_proximal_distal.gene_summary[pos, c(1,14)]
head(pos_data)
colnames(pos_data) <-  c("id" ,"lfc"   )
colnames(neg_data) <-  c("id" ,"lfc"   )
plot_data <-rbind(neg_data, pos_data)
#View(plot_data)
plot_data$group <- ifelse(plot_data$lfc < 0, "neg", "pos")
reduced_plot_data <- plot_data %>%
  group_by(group) %>%
  top_n(n = 10, wt = abs(lfc))
reduced_plot_data

ggplot(data=reduced_plot_data, aes(x=reorder(id, lfc),y= lfc,  fill=group)) +
  geom_bar(stat="identity")+
  theme_classic() + scale_fill_manual(values=c('#81C341','#D12026'))+ coord_flip()+
  ggsave("/media/Coco/MOSAIC LIVER/Manuscript/Graphs/SPH123_allpaired_lgfc.pdf", width = 5, height = 4)


####GSEA#####
library(fgsea)
library(msigdbr)
library(data.table)
library(stringr)
library(dplyr)
library(ggplot2)

###GSEA ON DIFFERENTIAL EXPRESSION ##
#select species and set
msigdbr_show_species()
m_df<- msigdbr(species = "Mus musculus", category = "H")
Hallmarks <- m_df %>% split(x = .$gene_symbol, f = .$gs_name)
head(Hallmarks)
m_df<- msigdbr(species = "Mus musculus", category = "C2", subcategory = "CGP")
CGP <- m_df %>% split(x = .$gene_symbol, f = .$gs_name)
head(CGP)
m_df<- msigdbr(species = "Mus musculus", category = "C2", subcategory = "CP:KEGG")
KEGG <- m_df %>% split(x = .$gene_symbol, f = .$gs_name)
head(KEGG)
m_df<- msigdbr(species = "Mus musculus", category = "C2", subcategory = "CP:REACTOME")
REACTOME <- m_df %>% split(x = .$gene_symbol, f = .$gs_name)
head(REACTOME)
m_df<- msigdbr(species = "Mus musculus", category = "C3", subcategory = "TFT")
TFT <- m_df %>% split(x = .$gene_symbol, f = .$gs_name)
head(TFT)
m_df<- msigdbr(species = "Mus musculus", category = "C5", subcategory = "BP")
GOBP <- m_df %>% split(x = .$gene_symbol, f = .$gs_name)
head(GOBP)

head(gdata)
ranks <- gdata %>% 
  na.omit()%>%
  mutate(ranking=-log10(FDR)/sign(Score))
ranks <- ranks$ranking
names(ranks) <- gdata$id
head(ranks, 10)

BP<- fgsea(pathways = GOBP, 
               stats = ranks,
               minSize=10,
               maxSize=500,
               nperm=1000000)
View(BP)

BP$pathway<-gsub("GOBP_","",BP$pathway)
BP$pathway<-gsub("_"," ",BP$pathway)
BP$pathway <- tolower(BP$pathway)
head(BP)

ggplot(BP %>% filter(abs(NES)>1) %>% head(n= 20), aes(reorder(pathway, NES), NES)) +
  #geom_bar(stat="identity")+
  geom_point(aes(size=size, colour=NES)) +
  scale_size_area(max_size = 8)+
  scale_colour_gradient( low = '#81C341',
    high = '#D12026',   aesthetics = "colour")+ coord_flip() +
 # ylim(-3.5, 0)+
    #scale_x_reverse()+
 # labs(x=" ", y="Normalized Enrichment Score",
#       title=" ", cols="black") + 
  theme_classic()  +
  ggsave("/media/Coco/MOSAIC LIVER/Manuscript/Graphs/SPH123_allpaired_GOBP.pdf", width = 10, height = 4)
 # labs(y="")+   
 # theme(axis.text.y=element_text(size=15), axis.text.x=element_text(size=15))




############CORRELATIONS############
SPH1_paired.sgRNA_summary <- read.delim("/media/Coco/MOSAIC LIVER/Experiments/Screen1/Fastqs/MAGECK/Bam_files/SPH1_paired.sgrna_summary.txt")
View(SPH1_paired.sgRNA_summary)
SPH1_paired.sgRNA_summary$ratio <- SPH1_paired.sgRNA_summary$control_count/SPH1_paired.sgRNA_summary$treatment_count
SPH2_paired.sgRNA_summary <- read.delim("/media/Coco/MOSAIC LIVER/Experiments/Screen1/Fastqs/MAGECK/Bam_files/SPH2_paired.sgrna_summary.txt")
View(SPH2_paired.sgRNA_summary)
SPH2_paired.sgRNA_summary$ratio <- SPH2_paired.sgRNA_summary$control_count/SPH2_paired.sgRNA_summary$treatment_count
SPH3_paired.sgRNA_summary <- read.delim("/media/Coco/MOSAIC LIVER/Experiments/Screen1/Fastqs/MAGECK/Bam_files/SPH3_paired.sgrna_summary.txt")
colnames(SPH3_paired.sgRNA_summary)
SPH3_paired.sgRNA_summary$ratio <- SPH3_paired.sgRNA_summary$control_count/SPH3_paired.sgRNA_summary$treatment_count

SPH12<-merge(SPH1_paired.sgRNA_summary[, c(1,16)], SPH2_paired.sgRNA_summary[, c(1,16)], by="sgrna")
SPH123 <-merge(SPH12, SPH3_paired.sgRNA_summary[, c(1,16)], by="sgrna")
head(SPH123)
colnames(SPH123) <- c("gene","SPH1", "SPH2", "SPH3")
plot(x=SPH123$SPH1 ,y= SPH123$SPH2, xlim= c(0,10))
plot(x=SPH123$SPH1 ,y= SPH123$SPH3,  xlim= c(0,10))
plot(x=SPH123$SPH2 ,y= SPH123$SPH3,   xlim= c(0,10))


########## NO CRE ANALYSIS ###########
#sum all mice for SPH1, 2 and 3
mageck count -l library2.csv -n SPH1nocre_sum --sample-label "distal,proximal" --fastq 3070distal.bam 3070proximal.bam  
mageck count -l library2.csv -n SPH2nocre_sum --sample-label "distal,proximal" --fastq 3378distal.bam,3381distal.bam 3378proximal.bam,3381proximal.bam 
mageck count -l library2.csv -n SPH3nocre_sum --sample-label "distal,proximal" --fastq 3428distal.bam,3431proximal.bam 3428proximal.bam,3431distal.bam  

`SPH1nocre_sum.count_normalized` <- read.delim("/media/Coco/MOSAIC LIVER/Experiments/Screen1/Fastqs/MAGECK/Bam_files/SPH1nocre_sum.count_normalized.txt")
`SPH2nocre_sum.count_normalized` <- read.delim("/media/Coco/MOSAIC LIVER/Experiments/Screen1/Fastqs/MAGECK/Bam_files/SPH2nocre_sum.count_normalized.txt")
`SPH3nocre_sum.count_normalized` <- read.delim("/media/Coco/MOSAIC LIVER/Experiments/Screen1/Fastqs/MAGECK/Bam_files/SPH3nocre_sum.count_normalized.txt") 

#paired analysis
SPHnocre <- list(`SPH1nocre_sum.count_normalized` ,`SPH2nocre_sum.count_normalized`,`SPH3nocre_sum.count_normalized` ) %>% 
  lapply(tibble::rownames_to_column) %>% purrr::reduce(full_join, by="sgRNA") %>% 
  # NB: below we are IMPUTING the genes where they miss in some samples, with value 0:
  # many ways possible, I choose one simple enough but also ok for big data sets:
  # https://stackoverflow.com/questions/8161836/how-do-i-replace-na-values-with-zeros-in-an-r-dataframe
  mutate_all(~replace(., is.na(.), 0))
head(SPHnocre)
SPHnocre_counts <- SPHnocre[, c(2,3,4,5,8,9,12,13)]
head(SPHnocre_counts)
colnames(SPHnocre_counts) <- c( "sgRNA" ,"Gene", "SPH1d","SPH1p","SPH2d","SPH2p", "SPH3d","SPH3p")
View(SPHnocre_counts)
write.table(SPHnocre_counts, "/media/Coco/MOSAIC LIVER/Experiments/Screen1/Fastqs/MAGECK/Bam_files/SPHnocre_counts.txt",append = FALSE, sep= "\t",  row.names = FALSE, quote = FALSE, col.names = TRUE)

mageck test -k SPHnocre_counts.txt -t SPH1p,SPH2p,SPH3p -c SPH1d,SPH2d,SPH3d -n SPH123nocre_paired --paired
SPH123nocre_paired.gene_summary <- read.delim("/media/Coco/MOSAIC LIVER/Experiments/Screen1/Fastqs/MAGECK/Bam_files/SPH123nocre_paired.gene_summary.txt")
View(SPH123nocre_paired.gene_summary)


#barplot by logfc > has nrp2
#has both negative and positive pvalue > take positive p val for positive logfc and neg for neg
neg<-which(SPH123nocre_paired.gene_summary$neg.lfc<0)
pos<-which(SPH123nocre_paired.gene_summary$pos.lfc>0)
colnames(SPH123nocre_paired.gene_summary)
neg_data <- SPH123nocre_paired.gene_summary[neg, c(1,8)]
head(neg_data)
pos_data <- SPH123nocre_paired.gene_summary[pos, c(1,14)]
head(pos_data)
colnames(pos_data) <-  c("id" ,"lfc"   )
colnames(neg_data) <-  c("id" ,"lfc"   )
plot_data <-rbind(neg_data, pos_data)
#View(plot_data)
plot_data$group <- ifelse(plot_data$lfc < 0, "neg", "pos")
reduced_plot_data <- plot_data %>%
  group_by(group) %>%
  top_n(n = 10, wt = abs(lfc))
reduced_plot_data

ggplot(data=reduced_plot_data, aes(x=reorder(id, lfc),y= lfc,  fill=group)) +
  geom_bar(stat="identity")+
  theme_classic() + scale_fill_manual(values=c('#81C341','#D12026'))+ coord_flip()#+
  #ggsave("/media/Coco/MOSAIC LIVER/Manuscript/Graphs/SPH123lgfc.pdf", width = 5, height = 4)


gdata = ReadRRA(SPH123nocre_paired.gene_summary)
View(gdata)
gdata<- gdata[-which(gdata$id =="ctrl"),]

sdata = ReadsgRRA(SPH123_paired.sgRNA_summary)
View(sdata)
sgRankView(sdata, gene = c("Psen1","Nrp2",  "Plxnb2", "Gpc4", "App", "Saa1"))+ggsave("/media/Coco/MOSAIC LIVER/Manuscript/Graphs/sgRNAplot.pdf", width = 5, height = 4)


gdata$LogFDR = -log10(gdata$FDR)
ScatterView(gdata, x = "Score", y = "LogFDR", label = "id", 
            model = "volcano", top = 10, y_cut=0.25, x_cut = 0.05)
VolcanoView(gdata, x = "Score", y = "FDR", Label = "id", top = 2, x_cut = 0.05, y_cut = 0.8, alpha=1)+
  theme_classic() + scale_fill_manual(values=c('#81C341', "grey", '#D12026'))+ ylim(0,3)+ xlim(-2,2) +
  ggsave("/media/Coco/MOSAIC LIVER/Manuscript/Graphs/Volcanonocre.pdf", width = 5, height = 4)



######### SUMMED ANALYSIS: SPH1 and SPH2 and SPH3##########
mageck count -l library2.csv -n SPH_sum --sample-label "distal,proximal" --fastq 3174distal.bam,3068distal.bam,3039distal.bam,3379distal.bam,3425distal.bam,3434distal.bam,3436proximal.bam 3174proximal.bam,3068proximal.bam,3039proximal.bam,3379proximal.bam,3425proximal.bam,3434proximal.bam,3436distal.bam  
SPH_sum.count_normalized <- read.delim("~/NAS/Coco/MOSAIC LIVER/Experiments/Screen1/Fastqs/MAGECK/Bam_files/SPH_sum.count_normalized.txt", header=T)
View(SPH_sum.count_normalized)
SPH_sum.count_normalized$ratio <- SPH_sum.count_normalized$proximal/SPH_sum.count_normalized$distal

SPH_sum.countsummary <- read.delim("/media/Coco/MOSAIC LIVER/Experiments/Screen1/Fastqs/MAGECK/Bam_files/SPH_sum.countsummary.txt")
View(SPH_sum.countsummary)

a<- ggplot(data=SPH_sum.countsummary, aes(x=File, y=Percentage, fill=Label)) +
  geom_bar(stat="identity", color="white", position=position_dodge())+
  theme_classic()+ scale_fill_manual(values=c('#81C341','#D12026'))+ ylim(0,1)+
  scale_x_discrete(guide = guide_axis(angle = 45))

b <-ggplot(data=SPH_sum.countsummary, aes(x=File, y=Zerocounts, fill=Label)) +
  geom_bar(stat="identity", color="white", position=position_dodge())+ ylim(0,50)+
  theme_classic()+ scale_fill_manual(values=c('#81C341','#D12026'))+
  scale_x_discrete(guide = guide_axis(angle = 45))

c <- ggplot(data=SPH_sum.countsummary, aes(x=File, y=GiniIndex, fill=Label)) +
  geom_bar(stat="identity", color="white", position=position_dodge())+ ylim(0,1)+
  theme_classic()+ scale_fill_manual(values=c('#81C341','#D12026'))+
  scale_x_discrete(guide = guide_axis(angle = 45))

ggarrange(a,b,c, ncol = 1, nrow = 3) +ggsave("/media/Coco/MOSAIC LIVER/Manuscript/Graphs/stats.pdf", width = 5, height = 11)

SPH_sum.count_normalized <- read.delim("/media/Coco/MOSAIC LIVER/Experiments/Screen1/Fastqs/MAGECK/Bam_files/SPH_sum.count_normalized.txt")
View(SPH_sum.count_normalized)
ratio <- SPH_sum.count_normalized
ratio$ratio <-SPH_sum.count_normalized$proximal/SPH_sum.count_normalized$distal
View(ratio) #this is a sgRNA level

#use mageck to perform robust rank aggregation and look at gene level
mageck test -k SPH_sum.count_normalized.txt -t proximal -c distal -n SPH_sum
SPH_sum.gene_summary <- read.delim("/media/Coco/MOSAIC LIVER/Experiments/Screen1/Fastqs/MAGECK/Bam_files/SPH_sum.gene_summary.txt")
View(SPH_sum.gene_summary)


ggplot(data=SPH_sum.gene_summary, aes(x=neg.lfc, y=-log10(neg.p.value), label=id)) + 
  geom_vline(xintercept=c(-0.5, 0.5), col="red") + geom_text()+
  geom_hline(yintercept=-log10(0.05), col="red")+
  geom_point() + 
  theme_classic()
ggplot(data=SPH_sum.gene_summary, aes(x=pos.lfc, y=-log10(pos.p.value), label=id)) + 
  geom_vline(xintercept=c(-0.5, 0.5), col="red") + geom_text()+
  geom_hline(yintercept=-log10(0.05), col="red")+
  geom_point() + 
  theme_classic()

#barplot of pvalue
#has both negative and positive pvalue > take positive p val for positive logfc and neg for neg
neg<-which(SPH2_paired.gene_summary$neg.lfc<0)
pos<-which(SPH2_paired.gene_summary$pos.lfc>0)
colnames(SPH2_paired.gene_summary)
neg_data <- SPH2_paired.gene_summary[neg, c(1,4, 8)]
head(neg_data)
pos_data <- SPH2_paired.gene_summary[pos, c(1,10, 14)]
head(pos_data)
colnames(pos_data) <-  c("id" ,"p.value" ,"lfc"   )
colnames(neg_data) <-  c("id" ,"p.value" ,"lfc"   )
pos_data$log10pval <- -log10(pos_data$p.value)
neg_data$log10pval <- -log10(neg_data$p.value)
neg_data$log10pval <- neg_data$log10pval*-1
plot_data <-rbind(neg_data, pos_data)
View(plot_data)
plot_data$group <- ifelse(plot_data$lfc < 0, "neg", "pos")
reduced_plot_data <- plot_data %>%
  group_by(group) %>%
  top_n(n = 15, wt = abs(log10pval))
reduced_plot_data

ggplot(data=reduced_plot_data, aes(x=reorder(id, log10pval),y= log10pval,  fill=group)) +
  geom_bar(stat="identity")+
  theme_classic() + scale_fill_manual(values=c('#81C341','#D12026'))+ coord_flip()+
  ggsave("/media/Coco/MOSAIC LIVER/Manuscript/Graphs/SPH2.pdf", width = 5, height = 4)

######### SUM ACROSS BATCHED AND PAIRED FOR FINAL PLOT fig 3#####
#sum all mice for SPH1, 2 and 3
mageck count -l library2.csv -n SPH1_sum --sample-label "distal,proximal" --fastq 3174distal.bam,3068distal.bam 3174proximal.bam,3068proximal.bam  
mageck count -l library2.csv -n SPH2_sum --sample-label "distal,proximal" --fastq 3039distal.bam,3379distal.bam,3425distal.bam 3039proximal.bam,3379proximal.bam,3425proximal.bam  
mageck count -l library2.csv -n SPH3_sum --sample-label "distal,proximal" --fastq 3434distal.bam,3436proximal.bam 3434proximal.bam,3436distal.bam  #3436 was switched in lib prep so here reverse

`SPH1_sum.count_normalized` <- read.delim("/media/Coco/MOSAIC LIVER/Experiments/Screen1/Fastqs/MAGECK/Bam_files/SPH1_sum.count_normalized.txt")
`SPH2_sum.count_normalized` <- read.delim("/media/Coco/MOSAIC LIVER/Experiments/Screen1/Fastqs/MAGECK/Bam_files/SPH2_sum.count_normalized.txt")
`SPH3_sum.count_normalized` <- read.delim("/media/Coco/MOSAIC LIVER/Experiments/Screen1/Fastqs/MAGECK/Bam_files/SPH3_sum.count_normalized.txt") 

#paired analysis
SPH <- list(`SPH1_sum.count_normalized` ,`SPH2_sum.count_normalized`,`SPH3_sum.count_normalized` ) %>% 
  lapply(tibble::rownames_to_column) %>% purrr::reduce(full_join, by="sgRNA") %>% 
  # NB: below we are IMPUTING the genes where they miss in some samples, with value 0:
  # many ways possible, I choose one simple enough but also ok for big data sets:
  # https://stackoverflow.com/questions/8161836/how-do-i-replace-na-values-with-zeros-in-an-r-dataframe
  mutate_all(~replace(., is.na(.), 0))
head(SPH)
SPH_counts <- SPH[, c(2,3,4,5,8,9,12,13)]
head(SPH_counts)
colnames(SPH_counts) <- c( "sgRNA" ,"Gene", "SPH1d","SPH1p","SPH2d","SPH2p", "SPH3d","SPH3p")
head(SPH_counts)
write.table(SPH_counts, "/media/Coco/MOSAIC LIVER/Experiments/Screen1/Fastqs/MAGECK/Bam_files/SPH_counts.txt",append = FALSE, sep= "\t",  row.names = FALSE, quote = FALSE, col.names = TRUE)

mageck test -k SPH_counts.txt -t SPH1p,SPH2p,SPH3p -c SPH1d,SPH2d,SPH3d -n SPH123_paired --paired
SPH123_paired.gene_summary <- read.delim("/media/Coco/MOSAIC LIVER/Experiments/Screen1/Fastqs/MAGECK/Bam_files/SPH123_paired.gene_summary.txt")
View(SPH123_paired.gene_summary)

ggplot(data=SPH123_paired.gene_summary, aes(x=neg.lfc, y=-log10(neg.p.value), label=id)) + 
  geom_vline(xintercept=c(-0.5, 0.5), col="red") + geom_text()+
  geom_hline(yintercept=-log10(0.05), col="red")+
  geom_point() + 
  theme_classic()

SPH123_paired.gene_summary$genelabels <- ""
SPH123_paired.gene_summary$genelabels <- ifelse(SPH123_paired.gene_summary$pos.rank < 11, TRUE,FALSE)

ggplot(data=SPH123_paired.gene_summary, 
       aes(x=id, y=pos.lfc, label = genelabels, color=-log10(pos.p.value))) + geom_point()+ #size=-log10(pos.p.value)
        geom_text_repel(label = ifelse(SPH123_paired.gene_summary$genelabels == TRUE, as.character(SPH123_paired.gene_summary$id),""), size=4, max.overlaps = Inf, segment.size = 0.2, color="black")+
        scale_color_gradientn(colours = c( "grey", '#D12026')) 
 
       label = genelabels, color=log10(neg.p.value)) +
  geom_point(shape = 20) + scale_y_reverse(limits = c(0,-5)) + theme_bw() + scale_size_area(max_size = 5)+
  scale_color_gradientn(colours = c("#FF0000", "black")) +
  theme(axis.text.x=element_blank(), axis.ticks.x=element_blank())+
  scale_fill_gradientn(colours = c("#FF0000", "black"))+
  geom_label_repel(label = ifelse(DEGs$genelabels == TRUE, as.character(DEGs$id),""), size=4, max.overlaps = Inf, segment.size = 0.2, color="black", fill="white")


ggplot(data=SPH123_paired.gene_summary, 
       aes(x=id, y=pos.lfc, label = genelabels, color=-log10(pos.p.value))) + geom_point()+ #size=-log10(pos.p.value)
  geom_text_repel(label = ifelse(SPH123_paired.gene_summary$genelabels == TRUE, as.character(SPH123_paired.gene_summary$id),""), size=4, max.overlaps = Inf, segment.size = 0.2, color="black")+
  scale_color_gradientn(colours = c( "grey", '#D12026')) 


ggplot(data=SPH123_paired.gene_summary, aes(x=pos.lfc, y=-log10(pos.p.value), label=id)) + 
  geom_vline(xintercept=c(-0.5, 0.5), col="red") + geom_text()+
  geom_hline(yintercept=-log10(0.05), col="red")+
  geom_point() + 
  theme_classic()

ggplot(data=SPH123_paired.gene_summary, aes(x=pos.rank, y=pos.p.value))+geom_point()
head(SPH123_paired.gene_summary)
       
#barplot of score > has Taz
#has both negative and positive pvalue > take positive p val for positive logfc and neg for neg
SPH123_paired.gene_summary <- SPH123_paired.gene_summary[-which(SPH123_paired.gene_summary$id=="ctrl"),]
neg<-which(SPH123_paired.gene_summary$neg.rank<11)
pos<-which(SPH123_paired.gene_summary$pos.rank<11)
colnames(SPH123_paired.gene_summary)
neg_data <- SPH123_paired.gene_summary[neg, c(1,3)]
head(neg_data)
pos_data <- SPH123_paired.gene_summary[pos, c(1,9)]
head(pos_data)
colnames(pos_data) <-  c("id" ,"score")
colnames(neg_data) <-  c("id" ,"score")
pos_data$log10score <- -log10(pos_data$score)
pos_data$group <- "pos"
neg_data$log10score <- -log10(neg_data$score)
neg_data$log10score <- neg_data$log10score*-1
neg_data$group <- "neg"
plot_data <-rbind(neg_data, pos_data)
head(plot_data)
ggplot(data=plot_data, aes(x=reorder(id, log10score),y= log10score,  fill=group)) +
  geom_bar(stat="identity")+
  theme_classic() + scale_fill_manual(values=c('#81C341','#D12026'))+ coord_flip()+
  ggsave("/media/Coco/MOSAIC LIVER/Manuscript/Graphs/SPH123.pdf", width = 5, height = 4)

#barplot by logfc > has nrp2
#has both negative and positive pvalue > take positive p val for positive logfc and neg for neg
neg<-which(SPH123_paired.gene_summary$neg.lfc<0)
pos<-which(SPH123_paired.gene_summary$pos.lfc>0)
colnames(SPH123_paired.gene_summary)
neg_data <- SPH123_paired.gene_summary[neg, c(1,8)]
head(neg_data)
pos_data <- SPH123_paired.gene_summary[pos, c(1,14)]
head(pos_data)
colnames(pos_data) <-  c("id" ,"lfc"   )
colnames(neg_data) <-  c("id" ,"lfc"   )
plot_data <-rbind(neg_data, pos_data)
#View(plot_data)
plot_data$group <- ifelse(plot_data$lfc < 0, "neg", "pos")
reduced_plot_data <- plot_data %>%
  group_by(group) %>%
  top_n(n = 10, wt = abs(lfc))
reduced_plot_data

ggplot(data=reduced_plot_data, aes(x=reorder(id, lfc),y= lfc,  fill=group)) +
  geom_bar(stat="identity")+
  theme_classic() + scale_fill_manual(values=c('#81C341','#D12026'))+ coord_flip()+
  ggsave("/media/Coco/MOSAIC LIVER/Manuscript/Graphs/SPH123lgfc.pdf", width = 5, height = 4)

####MAGeCK FLUTE####
library(MAGeCKFlute)
library(clusterProfiler)
library(ggplot2)
SPH123_paired.gene_summary <- read.delim("/media/Coco/MOSAIC LIVER/Experiments/Screen1/Fastqs/MAGECK/Bam_files/SPH123_paired.gene_summary.txt")
View(SPH123_paired.gene_summary)
SPH123_paired.sgRNA_summary <- read.delim("/media/Coco/MOSAIC LIVER/Experiments/Screen1/Fastqs/MAGECK/Bam_files/SPH123_paired.sgrna_summary.txt")
View(SPH123_paired.sgRNA_summary)

# Run FluteRRA with both gene summary file and sgRNA summary file
FluteRRA(SPH123_paired.gene_summary, SPH123_paired.sgRNA_summary, proj="SPH", organism="mmu", outdir = "./")
#Error in download.file(entrezfile, tmpfile, quiet = TRUE) : cannot open URL 'ftp://ftp.ensembl.org/pub/release-109/xml/tsv/mus_musculus/Mus_musculus.GRCm38.109/xml.entrez.tsv.gz'

gdata = ReadRRA(SPH123_paired.gene_summary)
View(gdata)
gdata<- gdata[-which(gdata$id =="ctrl"),]

sdata = ReadsgRRA(SPH123_paired.sgRNA_summary)
View(sdata)
sgRankView(sdata, gene = c("Psen1","Nrp2",  "Plxnb2", "Gpc4", "App", "Saa1"))+ggsave("/media/Coco/MOSAIC LIVER/Manuscript/Graphs/sgRNAplot.pdf", width = 5, height = 4)


gdata$LogFDR = -log10(gdata$FDR)
ScatterView(gdata, x = "Score", y = "LogFDR", label = "id", 
                 model = "volcano", top = 10, y_cut=0.25, x_cut = 0.05)
VolcanoView(gdata, x = "Score", y = "FDR", Label = "id", top = 2, x_cut = 0.05, y_cut = 0.8, alpha=1)+
theme_classic() + scale_fill_manual(values=c('#81C341', "grey", '#D12026'))+ ylim(0,3)+ xlim(-2,2) #+
# ggsave("/media/Coco/MOSAIC LIVER/Manuscript/Graphs/Volcano.pdf", width = 5, height = 4)

gdata$Rank = rank(gdata$Score)
ScatterView(gdata, x = "Rank", y = "Score", label = "id", 
            top = 15, auto_cut_y = TRUE, ylab = "Log2FC", 
            groups = c("top", "bottom"), max.overlaps = Inf)

View(gdata)
plot_data<- gdata
plot_data$group <- ifelse(plot_data$Score < 0, "neg", "pos")
View(plot_data)
neg <- plot_data %>% top_n(n = -15, wt = Score)
neg
pos <- plot_data %>% top_n(n = 15, wt = Score)
pos
reduced_plot_data <- rbind(neg, pos)
reduced_plot_data

ggplot(data=reduced_plot_data, aes(x=reorder(id, Score),y= Score,  fill=group)) +
  geom_bar(stat="identity")+
  theme_classic() + scale_fill_manual(values=c('#81C341','#D12026'))+ coord_flip()#+
  #ggsave("/media/Coco/MOSAIC LIVER/Manuscript/Graphs/SPH123.pdf", width = 5, height = 4)
