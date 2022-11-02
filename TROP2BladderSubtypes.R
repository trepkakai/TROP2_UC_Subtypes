####################################################################################################
####################################################################################################
# Kai Trepka
# 2021-08-26
####################################################################################################
####################################################################################################
setwd('/Users/kaitrepka/Documents/LabWork/FengLab_Analyses/BladderCancer')
# Libraries
library(consensusMIBC)
library(rms)
library(stringr)
library(Biobase)
library(ggplot2)
library(survival)
library(data.table)
library(org.Hs.eg.db)
library(dplyr)
library(Cairo) 

####################################################################################################
####################################################################################################
# Directories and filenames
dir_results <- "~/Documents/LabWork/FengLab_Analyses/BladderCancer/results"
dir_data_in_helper <- "~/Documents/LabWork/FengLab_Analyses/BladderCancer"

fn_sjodahl2017 <- paste0(dir_data_in_helper,"/sjodahl2017_eset")
fn_sjodahl2017_suppl <- "~/Documents/LabWork/FengLab_Analyses/BladderCancer/GSE83586_Processed_Data_Supplement.txt"
fn_sjodahl2012 <- paste0(dir_data_in_helper,"/sjodahl2012_eset")
fn_seiler <- paste0(dir_data_in_helper,"/seiler_eset")
fn_tcga <- paste0(dir_data_in_helper,"/tcga_cell_2017_data")

fn_sbt <- "~/Documents/LabWork/FengLab_Analyses/BladderCancer/1-s2.0-S0302283819306955-mmc9.txt"
####################################################################################################
# Colors subtypes
sbt6_colors <- c("orangered2","royalblue2","olivedrab2","goldenrod2","darkorchid2","springgreen2")
names(sbt6_colors) <- c("Ba/Sq","LumNS","LumP","LumU","NE-like","Stroma-rich")
colors_low_high <- c("dodgerblue3", "red3")

####################################################################################################
####################################################################################################
# Define custom functions
####################################################################################################
####################################################################################################
# KM survival plot for quartiles of gene expression, add p-val and HR (per quartile step)
plotsurvQuartiles <- function(gex_vector,
                              col_time, 
                              col_event,
                              x_lab="Years",
                              y_lab="", 
                              t_lab="",
                              text_x = 1,
                              text_y = 0.2,
                              plot_at_risk = TRUE,
                              colors=c("olivedrab3", "royalblue2","goldenrod3","firebrick2"),
                              ltypes=c(1,1,1,1),lwid=c(2,2,2,2),
                              x_max=NULL,
                              y_min=NULL,
                              time_inc=24,
                              savePlot = F,
                              filename = NULL,
                              fileH = 7,
                              fileW = 10) {
  
  stopifnot(length(col_event)==length(gex_vector))
  
  # Create quartile grouping
  quartiles <- rep(NA, length(gex_vector))
  quartiles[gex_vector < quantile(gex_vector,0.25)] <- 1
  quartiles[gex_vector >= quantile(gex_vector,0.25)] <- 2
  quartiles[gex_vector >= quantile(gex_vector,0.5)] <- 3
  quartiles[gex_vector >= quantile(gex_vector,0.75)] <- 4
  
  survival <- Surv(col_time, col_event)
  
  cox_cat <- coxph(survival~quartiles)
  hr <- coef(summary(cox_cat))[1,2]
  ci <- exp(confint(cox_cat))
  pval <- summary(cox_cat)$coefficients[,'Pr(>|z|)']
  
  # Factoriation is done after cox to get continous cox results
  labls <- c("Q1","Q2","Q3","Q4")
  quartiles <- factor(quartiles,levels = c(1:4), labels = labls)
  
  fit <- npsurv(survival ~ quartiles)
  
  if(is.null(y_min)){
    yLim <- c(0,1)
  }else{
    yLim <- c(y_min,1)
  }
  
  if(is.null(x_max)){
    xLim <- c(0,pretty(fit$maxtime)[2])
  }else{
    xLim <- c(0,x_max)
  }
  
  par(mar=c(11,4.1,4.1,2.2))
  survplot(fit, xlab="", ylab=y_lab, lwd=lwid, col=colors, n.risk = plot_at_risk, adj.n.risk=1, sep.n.risk=0.05, y.n.risk = -0.5, time.inc=time_inc, lty=ltypes,
           label.curves=list(keys = "lines",cex=1), levels.only=TRUE, conf="none", xlim=xLim, ylim=yLim,
           cex.n.risk=1,cex.xlab=1,cex.ylab=1)
  
  numstring <- paste0("HR = ", signif(hr,2), " (95%CI: ", signif(ci[1],2), "-", signif(ci[2],2), ")\nP = ", signif(pval, 2), sep="")
  
  text(text_x, text_y, numstring, pos=4, cex=1)
  mtext(text = x_lab, side = 1, line = 2, at = (x_max+1)/2)
  mtext(text = "Number at risk", side = 1, line = 3, at = -0.3, adj = 0)
  title(main=paste(t_lab),cex.main=1.7)
  
  if(savePlot){
    for(fileType in c("png","pdf","svg")){
      if(fileType == "pdf"){
        pdf(file = paste0(filename,".",fileType), width = fileW, height = fileH)
      }else if(fileType == "png"){
        png(file = paste0(filename,".",fileType), width = fileW, height = fileH, units = "in", res = 300)
      }else{
        CairoSVG(file = paste0(filename,".",fileType), width = fileW, height = fileH)
      }
      par(mar=c(11,4.1,4.1,2.2))
      survplot(fit, xlab="", ylab=y_lab, lwd=lwid, col=colors,n.risk = plot_at_risk, adj.n.risk=1, sep.n.risk=0.05, y.n.risk = -0.5, time.inc=time_inc, lty=ltypes,
               label.curves=list(keys = "lines",cex=1), levels.only=TRUE, conf="none", xlim=xLim, ylim=yLim,
               cex.n.risk=1,cex.xlab=1,cex.ylab=1)
      
      numstring <- paste0("HR = ", signif(hr,2), " (95%CI: ", signif(ci[1],2), "-", signif(ci[2],2), ")\nP = ", signif(pval, 2), sep="")
      
      text(text_x, text_y, numstring, pos=4, cex=1)
      mtext(text = x_lab, side = 1, line = 2, at = (x_max+1)/2)
      mtext(text = "Number at risk", side = 1, line = 3, at = -0.3, adj = 0)
      title(main=paste(t_lab),cex.main=1.7)
      
      dev.off()
    }
  }
  
}

# Plot correlation of gene expression of 2 genes, add p and Rho values, possibility to color by group vector, and and optional ab-line
plot2GenesGG <- function(gene1,
                         gene2,
                         rnaData,
                         x_text = NULL,
                         y_text = NULL,
                         Title = NULL,
                         sub = "", 
                         xLab = NULL,
                         yLab = NULL,
                         Colors = NULL,
                         ColorGroup = NULL,
                         abline = FALSE,
                         LOG2 = FALSE,
                         base_fn = NULL,
                         fileH = 5,
                         fileW = 5,
                         savePlot = F){
  
  expr1 <- rnaData[gene1,]
  expr2 <- rnaData[gene2,]
  
  if(LOG2){
    expr1 <- log2(expr1 + 1)
    expr2 <- log2(expr2 + 1)
  }
  
  df <- data.frame(x = expr1,
                   y = expr2)
  
  if(!is.null(ColorGroup)){
    df$Color <- as.factor(ColorGroup)
  }
  
  C <- signif(cor(expr1, expr2, use = "complete.obs", method = "spearman"),2)
  
  if(cor.test(expr1, expr2, use = "complete.obs", method = "spearman")$p.value < 0.0001){
    textString <- paste0("Spearman's Rho: ",C,
                         "\n", "P < 0.0001")
  }else{
    P <- signif(cor.test(expr1, expr2, use = "complete.obs", method = "spearman")$p.value,2)
    textString <- paste0("Spearman's Rho: ",C,
                         "\n", "P = ",P)
  }
  
  if(is.null(x_text)){
    x_text <- min(expr1) + ((max(expr1) - min(expr2))*0.8) 
  }
  if(is.null(y_text)){
    y_text <- min(expr2) + ((max(expr2) - min(expr2))*0.2) 
  }
  
  if(is.null(Title)){
    Title = paste0(gene1, " vs. ", gene2)
  }
  
  if(is.null(xLab)){
    xLab <- paste0(gene1," expression")
  }
  
  if(is.null(yLab)){
    yLab <- paste0(gene2," expression")
  }
  
  p <- ggplot(df, aes(x = x, y = y)) +
    geom_point(aes(color = ColorGroup)) +
    labs(title = Title, subtitle = sub, x = xLab, y = yLab) +
    theme(axis.text.x = element_text(size = 12),
          axis.text.y = element_text(size = 12),
          plot.title = element_text(hjust = 0.5),
          plot.subtitle = element_text(hjust = 0.5),
          panel.background = element_rect(fill = "white"),
          panel.grid.major = element_line(size = 0.5, linetype = 'solid',
                                          colour = "gray"), 
          panel.grid.minor = element_line(size = 0.25, linetype = 'solid',
                                          colour = "white"),
          axis.line = element_line()) +
    annotate("text", x = x_text, y = y_text, label = textString)
  
  if(abline){
    p <- p + geom_abline(linetype = "dashed")
  }
  
  if(!is.null(ColorGroup)){
    p <- p + scale_color_manual(values=Colors) +
      labs(color = "Subtype") +
      theme(legend.position = c(2,2), # Originally (0.1, 0.8) - ALTER THIS TO CHANGE LEGEND 
            legend.background = element_rect(fill = alpha("white",1)),
            legend.key = element_rect(fill = alpha("white",0)))
    
  }
  
  print(p)
  
  if(savePlot == T){
    for(fileType in c("png","pdf")){
      if(fileType == "pdf"){
        pdf(file = paste0(base_fn,".",fileType), width = fileW, height = fileH)
      }else{
        png(file = paste0(base_fn,".",fileType), width = fileW, height = fileH, units = "in", res = 300)
      }
      print(p)
      dev.off()
    }
  }
}

# Violin plot to plot expression values per subtype
plotViolin <- function(gene,
                       rnaData, 
                       groupVariable, 
                       groupNames = NULL, 
                       yLab = "",
                       y_text = 7, 
                       x_text = 1.5,
                       Colors = NULL, 
                       Title = NULL,
                       sub = "Dataset",
                       base_fn = NULL,
                       savePlot = F,
                       fileW = 10,
                       fileH = 7){
  
  if(!is.factor(groupVariable)){
    groupVariable <- as.factor(groupVariable)
  }
  
  if(is.null(groupNames)){
    groupNames = levels(groupVariable)
  }
  
  
  #create data frame for ggplot2 functionality
  theData <- data.frame(expression = rnaData[gene,],
                        group = groupVariable)
  
  AnovaPval <- summary(aov(rnaData[gene,] ~ groupVariable))[[1]]["groupVariable","Pr(>F)"]
  KruskalWallisPval <- kruskal.test(rnaData[gene,] ~ groupVariable)$p.value
  pvals <- pairwise.wilcox.test(rnaData[gene,],groupVariable, p.adjust.method = "none")$p.value
  
  sink(file = paste0(base_fn,"_pvals.txt"))
  cat("ANOVA p-value: ")
  cat(AnovaPval)
  cat("\n")
  cat("Kruskal-Wallis p-value: ")
  cat(KruskalWallisPval)
  cat("\n")
  cat("P-values for pairwise Wilcoxon rank-sum test")
  cat("\n")
  print(pvals)
  cat("\n")
  cat("\n")
  cat("Descriptive statistics")
  cat("\n")
  print(theData %>%
          group_by(group) %>% 
          summarize(median = median(expression),
                    sd = sd(expression),
                    Q1 =  quantile(expression, 0.25),
                    Q3 =  quantile(expression, 0.75),
                    min = min(expression),
                    max = max(expression))
  )
  
  sink()
  
  if(!is.null(Colors)){
    stopifnot(names(Colors)==levels(groupVariable))
  }else{
    Colors = 1:length(levels(theData$group))
  }
  
  if(is.null(Title)){
    Title = gene
  }
  
  p = ggplot(theData, aes(x = group, y = expression, fill = group)) +
    geom_violin() + 
    #geom_dotplot(binaxis='y', stackdir='center', dotsize=0.5) #+
    geom_boxplot(width=0.05, fill = "white") +
    geom_jitter(width = 0.2) +
    scale_fill_manual(values=Colors) +
    scale_x_discrete(labels = groupNames)+
    labs(title = Title, subtitle = sub, x = "", y = yLab) +
    theme(axis.text.x = element_text(size = 12),
          plot.title = element_text(hjust = 0.5),
          plot.subtitle = element_text(hjust = 0.5),
          panel.background = element_rect(fill = "white"),
          panel.grid.major = element_line(size = 0.5, linetype = 'solid',
                                          colour = "gray"), 
          panel.grid.minor = element_line(size = 0.25, linetype = 'solid',
                                          colour = "white"),
          legend.position = "none",
          axis.line = element_line()) 
  
  
  print(p)
  
  if(savePlot == T){
    for(fileType in c("png","pdf","svg")){
      if(fileType == "pdf"){
        pdf(file = paste0(base_fn, ".",fileType), width = fileW, height = fileH)
      }else if(fileType == "png"){
        png(file = paste0(base_fn,".",fileType), width = fileW, height = fileH, units = "in", res = 300)
      }else{
        CairoSVG(file = paste0(base_fn,".",fileType), width = fileW, height = fileH)
      }
      
      print(p)
      
      dev.off()
    }
  }
}

####################################################################################################
####################################################################################################
# Load data
####################################################################################################
####################################################################################################
# Read in data
sjodahl17_eset <- readRDS(fn_sjodahl2017)
sjodahl12_eset <- readRDS(fn_sjodahl2012)
seiler_eset <- readRDS(fn_seiler)

tcga_list <- readRDS(fn_tcga)
tcga.rna <- tcga_list[["tcga.rna"]]
tcga.rna.mapping <- tcga_list[["tcga.rna.mapping"]]
tcga.clin <- tcga_list[["tcga.clin"]]

sbt_master <- read.table(file = fn_sbt, header = T, sep = "\t", stringsAsFactors = F)

####################################################################################################
# Format and add subtype
# Sjodahl 2017
# Data is  RMA normalized but not median centered. And not batch effect corrected?

# Add subtype from paper
pData(sjodahl17_eset)$ID <- paste0("Sjo.",pData(sjodahl17_eset)$title)
table(pData(sjodahl17_eset)$ID %in% sbt_master$ID) # 243 matches with the number in the paper
pData(sjodahl17_eset)$sbt <- NA
pData(sjodahl17_eset)$sbt <- sbt_master$Consensus_class[match(pData(sjodahl17_eset)$ID,sbt_master$ID)]

symbol <- str_split_fixed(fData(sjodahl17_eset)$gene_assignment, "//", n = 3)[,2]
symbol <- trimws(symbol)
fData(sjodahl17_eset)$gene <- symbol
entrez <- ifelse(grepl(pattern = "// [0-9]+$",  fData(sjodahl17_eset)$gene_assignment), str_extract(fData(sjodahl17_eset)$gene_assignment, pattern = "[0-9]+$"),NA)
fData(sjodahl17_eset)$Entrez_Gene_Id <- entrez 

# Subset for clinical characteristics:
sjodahl17_eset$stage <- str_split_fixed(sjodahl17_eset$`stage and grade:ch1`," ", n = 2)[,1] 
sjodahl17_eset$grade <- trimws(str_split_fixed(sjodahl17_eset$`stage and grade:ch1`," ", n = 2)[,2])

sjodahl17 <- sjodahl17_eset[,sjodahl17_eset$stage %in% c("pT2","pT3","pT4")] # gives 243, that are the same that had sbt called in the paper


# For comparision with the normalized sample table provided
sjodahl17_suppl_expr <- read.table(file = fn_sjodahl2017_suppl, stringsAsFactors = F, header = T, sep = "\t", skip = 10)
rownames(sjodahl17_suppl_expr) <- sjodahl17_suppl_expr$SYMBOL
sjodahl17_suppl_expr <- sjodahl17_suppl_expr[,3:ncol(sjodahl17_suppl_expr)]

sjodahl17_suppl_expr_pdata <- read.table(file = fn_sjodahl2017_suppl, stringsAsFactors = F, header = T, sep = "\t", nrows = 10)
sjodahl17_suppl_expr_pdata["ID",] <- sjodahl17_suppl_expr_pdata[1,]
sjodahl17_suppl_expr_pdata["ID",] <- paste0("Sjo.BLCA_Cx_",gsub("Sample ","",sjodahl17_suppl_expr_pdata["ID",])) 
sjodahl17_suppl_expr_pdata <- sjodahl17_suppl_expr_pdata[,2:ncol(sjodahl17_suppl_expr_pdata)]

sjodahl17_suppl_expr_pdata <- t(sjodahl17_suppl_expr_pdata)
colnames(sjodahl17_suppl_expr_pdata) <- sjodahl17_suppl_expr_pdata[1,] 
sjodahl17_suppl_expr_pdata <- sjodahl17_suppl_expr_pdata[-1,]
colnames(sjodahl17_suppl_expr_pdata)[ncol(sjodahl17_suppl_expr_pdata)] <- "ID"
rownames(sjodahl17_suppl_expr_pdata) <- sjodahl17_suppl_expr_pdata[,"ID"]
sjodahl17_suppl_expr_pdata <- as.data.frame(sjodahl17_suppl_expr_pdata)

colnames(sjodahl17_suppl_expr) <- sjodahl17_suppl_expr_pdata$ID
sjodahl17_suppl_expr <- as.matrix(sjodahl17_suppl_expr)

# Add subtype from paper

table(sjodahl17_suppl_expr_pdata[,"ID"] %in% sbt_master$ID) # 243 matches with the number in the paper
sjodahl17_suppl_expr_pdata$sbt <- NA
sjodahl17_suppl_expr_pdata$sbt <- sbt_master$Consensus_class[match(sjodahl17_suppl_expr_pdata$ID,sbt_master$ID)]
#table(sjodahl17_suppl_expr_pdata$sbt,sjodahl17_suppl_expr_pdata$ConsensusClusterPlus_Clusters)
#table(sjodahl17_suppl_expr_pdata[sjodahl17$ID,"sbt"], sbt_master[sjodahl17$ID,"Consensus_class"])

#sync
all(sjodahl17$ID %in% sjodahl17_suppl_expr_pdata$ID)
all(colnames(sjodahl17_suppl_expr)==rownames(sjodahl17_suppl_expr_pdata))

sjodahl17_suppl_expr_pdata <- sjodahl17_suppl_expr_pdata[sjodahl17$ID,]
sjodahl17_suppl_expr <- sjodahl17_suppl_expr[,sjodahl17$ID]

# sanity plots
#theGene <- "PVRL4"
#plot(exprs(sjodahl17)[fData(sjodahl17)$gene == theGene,],sjodahl17_suppl_expr[theGene,], 
#     col = as.factor(sjodahl17_suppl_expr_pdata$Labeling_Batch)) # probably the batch effect correctsion??

#table(sjodahl17_suppl_expr_pdata$Labeling_Batch == sjodahl17$`labeling batch:ch1`)
####################################################################################################
# Sjodahl 2012
# Data is processed for Illumina bead arrays, corrected for batch, background, QN, log2, mean centered and intensity filtered. 
# Add subtype from paper
sjodahl12_eset$ID <- gsub("\\.a1\\.lbe1","",sjodahl12_eset$title)
sjodahl12_eset$ID <- gsub("\\.RNA\\.e1","",sjodahl12_eset$ID)
sjodahl12_eset$ID <- gsub("\\.RNA\\.e2","",sjodahl12_eset$ID)
sjodahl12_eset$ID <- gsub("_h2","",sjodahl12_eset$ID)
sjodahl12_eset$ID <- paste0("Sjo.",sjodahl12_eset$ID)
table(sjodahl12_eset$ID %in% sbt_master$ID) # the supposed 93 do match

sjodahl12_eset$sbt <- NA
sjodahl12_eset$sbt <- sbt_master$Consensus_class[match(pData(sjodahl12_eset)$ID,sbt_master$ID)]

# Remove non MIBC
sjodahl12_eset$stage <- sjodahl12_eset$`tumor_stage:ch1`
sjodahl12_eset$grade <- sjodahl12_eset$`tumor_grade:ch1`
sjodahl12 <- sjodahl12_eset[,!sjodahl12_eset$stage %in% c("T1","Ta","Tx")] # gets the desired 93 with sbt data

####################################################################################################
# Seiler
# Human exon 1.0 ST, SCAN normalized
# Add subtype from paper
table(seiler_eset$geo_accession %in% sbt_master$ID) # the supposed 305 do match

seiler_eset$sbt <- NA
seiler_eset$sbt <- sbt_master$Consensus_class[match(seiler_eset$geo_accession,sbt_master$ID)]

seiler <- seiler_eset

####################################################################################################
# TCGA
# Add subtype from paper
dashToDot <- function(x){
  gsub("-","\\.",x)
}
sbt_master$ID_old <- sbt_master$ID
sbt_master$ID <- dashToDot(sbt_master$ID)
sbt_master$ID <- gsub("[A-Z]{1}$","",sbt_master$ID) #remove A at the end, one sample has B at the end

table(tcga.clin$SAMPLE_ID %in% sbt_master$ID) # the supposed 406 do match

# Probably most consistent to just go the 406 overlapping with the consensus paper. The original Cell paper states that it should be 406 MIBC with complete clinical data.
# However, there are inconsistensies regarding histological diagnsosis etc. Not really up to me to re-evaluate the TCGA dataset? And in any case, for this project, is has 0
# impact on the results and conclusions

tcga.clin$sbt <- NA
tcga.clin$sbt <- sbt_master$Consensus_class[match(tcga.clin$SAMPLE_ID,sbt_master$ID)]

tcga.clin <- tcga.clin[!is.na(tcga.clin$sbt),]

#Subset/check clinical characteristics
#tcga.clin <- tcga.clin[tcga.clin$SAMPLE_TYPE != "Metastasis",] # 1 met, 1 unknown, remove the met. Keep unknown since is still MIBC as diagnosis
table(tcga.clin$SAMPLE_TYPE, useNA = "ifany") # 1 unknown...

#tcga.clin <- tcga.clin[tcga.clin$HISTOLOGICAL_DIAGNOSIS == "Muscle invasive urothelial carcinoma (pT2 or above)" & !is.na(tcga.clin$HISTOLOGICAL_DIAGNOSIS),] # 4 NAs
table(tcga.clin$HISTOLOGICAL_DIAGNOSIS) # 3 missing, keep to be consistent with the consensuspaper
table(duplicated(tcga.clin$PATIENT_ID)) # no patients with >1 sample

table(tcga.clin$SAMPLE_ID %in% colnames(tcga.rna)) #all match
#tcga.clin <- tcga.clin[!is.na(tcga.clin$sbt),]

# Make clinical and rna data match
theSamples <- intersect(colnames(tcga.rna),tcga.clin$SAMPLE_ID)

tcga.rna <- tcga.rna[,theSamples]
rownames(tcga.clin) <- tcga.clin$SAMPLE_ID
tcga.clin <- tcga.clin[theSamples,]

stopifnot(all(colnames(tcga.rna)==rownames(tcga.clin)))


# The below was for the first round of analyses. Now I use the same 406 as the consensuspaper.

# The survival data is for 404 paitents BUT 2 of these are missing...


# 3 samples have missing FU data, 400 remain for survival analyses. I think the discrepancy to the consensus paper (where 2 is missing FU data) is that I exlude 1 pat with negative OS time, which I think makes sense to assume is an error, or at least should be excluded.

# I end up with 403 samples, of which 400 have survival data.

# The discrepancy to the 406 is because I exclude the 3 that to not have MIBC as diagnosis. I think this makes sense and I can certainly explain why we did this.Somehow TCGA data needs to be filtered anyway, since there are mets, mutliple samples from 1 patients etc.


####################################################################################################
####################################################################################################
# Plotting
####################################################################################################
####################################################################################################
#Violin plots

# TROP2 (TACSTD2)
dir_results_1 <- "~/Documents/LabWork/FengLab_Analyses/BladderCancer/TROP2_Figures"

fn = paste0(dir_results_1, "/TROP2_violinplot_seiler")
plotViolin(gene = "TACSTD2", rnaData = exprs(seiler), groupVariable =  pData(seiler)$sbt, yLab = "Gene expression (SCAN)", Colors = sbt6_colors,
           Title = "TROP2", sub = "Seiler data", base_fn = fn, fileW = 8, fileH = 4, savePlot = T)

fn = paste0(dir_results_1, "/TROP2_violinplot_sjodahl2017")
plotViolin(gene = which(fData(sjodahl17)$gene=="TACSTD2"),
           rnaData = exprs(sjodahl17), groupVariable =  sjodahl17$sbt, yLab = "Gene expression (RMA normalized)", Colors = sbt6_colors,
           Title = "TROP2", sub = "Sjodahl 2017 data", base_fn = fn, fileW = 8, fileH = 4, savePlot = T)

fn = paste0(dir_results_1, "/TROP2_violinplot_sjodahl2017_SUPPL_DATA")
plotViolin(gene = "TACSTD2",
           rnaData = as.matrix(sjodahl17_suppl_expr), groupVariable = sjodahl17_suppl_expr_pdata$sbt, yLab = "Gene expression (RMA normalized)", Colors = sbt6_colors,
           Title = "TROP2", sub = "Sjodahl 2017 SUPPLEMENTAL DATA", base_fn = fn, fileW = 8, fileH = 4, savePlot = T)

fn = paste0(dir_results_1,"/TROP2_violinplot_tcga")
plotViolin(gene = "TACSTD2",
           rnaData = tcga.rna, groupVariable =  tcga.clin$sbt, yLab = "Gene expression log2(RSEM+1)", Colors = sbt6_colors,
           Title = "TROP2", sub = "TCGA Cell 2017 data", base_fn = fn, fileW = 8, fileH = 4, savePlot = T)

fn = paste0(dir_results_1, "/TROP2_violinplot_sjodahl2012")
plotViolin(gene = which(fData(sjodahl12)$Symbol=="TACSTD2"), rnaData = exprs(sjodahl12), groupVariable =  sjodahl12$sbt,
           yLab = "Gene expression (QN, log2, mean centered)", Colors = sbt6_colors,Title = "TROP2", sub = "Sjodahl 2012 data",
           base_fn = fn, fileW = 8, fileH = 4, savePlot = T)


####################################################################################################
# Correlation with other genes
# TROP2 vs NECTIN 4
plot2GenesGG(gene1 = "TACSTD2",
             gene2 = "PVRL4",rnaData = tcga.rna,sub = "TCGA Cell 2017 data", ColorGroup = tcga.clin$sbt, Colors = sbt6_colors,abline = F,
             Title = "TROP2 vs. NECTIN4", base_fn = paste0(dir_results, "/cdcp1/correlation/corr_TROP2_NECTIN4_tcga"), savePlot = T, 
             x_text = 8, y_text = 14, xLab = "TROP2 expression")

plot2GenesGG(gene1 = "TACSTD2",
             gene2 = "PVRL4",rnaData = exprs(seiler),sub = "Seiler data", ColorGroup = seiler$sbt, Colors = sbt6_colors,abline = F,
             Title = "TROP2 vs. NECTIN4", base_fn = paste0(dir_results, "/cdcp1/correlation/corr_TROP2_NECTIN4_seiler"), 
             savePlot = T, x_text = 0.5, y_text = 2, xLab = "TROP2 expression")

plot2GenesGG(gene1 = which(fData(sjodahl17)$gene=="TACSTD2"), gene2 = which(fData(sjodahl17)$gene=="PVRL4"), rnaData = exprs(sjodahl17),
             sub = "Sjodahl 2017 data", ColorGroup = sjodahl17$sbt, Colors = sbt6_colors, abline = F,
             Title = "TROP2 vs. NECTIN4", base_fn = paste0(dir_results, "/cdcp1/correlation/corr_TROP2_NECTIN4_sjodahl2017"), 
             savePlot = T, x_text = 7, y_text = 8, xLab = "TROP2 expression", yLab = "NECTIN4 expression")

plot2GenesGG(gene1 = "TACSTD2",
             gene2 = "PVRL4", rnaData = sjodahl17_suppl_expr, sub = "Sjodahl 2017 SUPPLEMENTAL DATA", ColorGroup = sjodahl17_suppl_expr_pdata$sbt, 
             Colors = sbt6_colors,abline = F, Title = "TROP2 vs. NECTIN4", base_fn = paste0(dir_results, "/cdcp1/correlation/corr_TROP2_NECTIN4_sjodahl17_SUPPL_DATA"), 
             savePlot = T, x_text = -3, y_text = 1, xLab = "TROP2 expression")

plot2GenesGG(gene1 = which(fData(sjodahl12)$Symbol=="TACSTD2"), gene2 = which(fData(sjodahl12)$Symbol=="PVRL4"), rnaData = exprs(sjodahl12),
             sub = "Sjodahl 2012 data", ColorGroup = sjodahl12$sbt, Colors = sbt6_colors, abline = F,
             Title = "TROP2 vs NECTIN4", base_fn = paste0(dir_results, "/cdcp1/correlation/corr_TROP2_NECTIN4_sjodahl2012"), 
             savePlot = T, x_text = -2, y_text = 1, xLab = "TROP2 expression", yLab = "NECTIN4 expression")

# TROP2 vs PD1
plot2GenesGG(gene1 = "TACSTD2",
             gene2 = "PDCD1",rnaData = tcga.rna,sub = "TCGA Cell 2017 data", ColorGroup = tcga.clin$sbt, Colors = sbt6_colors,abline = F,
             Title = "TROP2 vs. PD1", base_fn = paste0(dir_results_1, "/revision_1/corr_TROP2_PD1_tcga"), savePlot = T, 
             x_text = 14, y_text = 12, xLab = "TROP2 expression")

plot2GenesGG(gene1 = "TACSTD2",
             gene2 = "PDCD1",rnaData = exprs(seiler),sub = "Seiler data", ColorGroup = seiler$sbt, Colors = sbt6_colors,abline = F,
             Title = "TROP2 vs. PD1", base_fn = paste0(dir_results_1, "/revision_1/corr_TROP2_PD1_seiler"), 
             savePlot = T, x_text = 2, y_text = 0.45, xLab = "TROP2 expression")

plot2GenesGG(gene1 = "TACSTD2",
             gene2 = "PDCD1", rnaData = sjodahl17_suppl_expr, sub = "Sjodahl 2017 SUPPLEMENTAL DATA", ColorGroup = sjodahl17_suppl_expr_pdata$sbt, 
             Colors = sbt6_colors,abline = F, Title = "TROP2 vs. PD1", base_fn = paste0(dir_results_1, "/revision_1/corr_TROP2_PD1_sjodahl17_SUPPL_DATA"), 
             savePlot = T, x_text = 0, y_text = 1.2, xLab = "TROP2 expression")


plot2GenesGG(gene1 = which(fData(sjodahl17)$gene=="TACSTD2"), gene2 = which(fData(sjodahl17)$gene=="PDCD1"), rnaData = exprs(sjodahl17),
             sub = "Sjodahl 2017 data", ColorGroup = sjodahl17$sbt, Colors = sbt6_colors, abline = F,
             Title = "TROP2 vs. PDCD1", base_fn = paste0(dir_results, "/cdcp1/correlation/corr_TROP2_PDCD1_sjodahl2017"), 
             savePlot = T, x_text = 7, y_text = 7.5, xLab = "TROP2 expression", yLab = "PDCD1 expression")


# TROP2 vs PD-L1
plot2GenesGG(gene1 = "TACSTD2",
             gene2 = "CD274",rnaData = tcga.rna,sub = "TCGA Cell 2017 data", ColorGroup = tcga.clin$sbt, Colors = sbt6_colors,abline = F,
             Title = "TROP2 vs. PD-L1", base_fn = paste0(dir_results_1, "/revision_1/corr_TROP2_PDL1_tcga"), savePlot = T, 
             x_text = 10, y_text = 10.5, xLab = "TROP2 expression")

plot2GenesGG(gene1 = "TACSTD2",
             gene2 = "CD274",rnaData = exprs(seiler),sub = "Seiler data", ColorGroup = seiler$sbt, Colors = sbt6_colors,abline = F,
             Title = "TROP2 vs. PD-L1", base_fn = paste0(dir_results_1, "/revision_1/corr_TROP2_PDL1_seiler"), 
             savePlot = T, x_text = 2, y_text = 1.2, xLab = "TROP2 expression")

##################################################################################################
# Other suggested revisions figures
# Sex. Sex not available for Seiler, Sjodahl 2017 datasets
fn = paste0(dir_results_1,"/revision_1/sex_TROP2_violinplot_tcga")
plotViolin(gene = "TACSTD2",
           rnaData = tcga.rna, groupVariable =  tcga.clin$SEX, yLab = "Gene expression log2(RSEM+1)", Colors = colors_low_high,
           Title = "TROP2", sub = "TCGA Cell 2017 data", base_fn = fn, fileW = 6, fileH = 4, savePlot = T)

# Stage
#tcga.clin$AJCC_PATHOLOGIC_TUMOR_STAGE # Pathologic stage
#seiler$`cstage:ch1` # Clinical stage, since pathologic stage not readily available
#sjodahl17$stage # Pathologic stage
fn = paste0(dir_results_1, "/revision_1/stagec_TROP2_violinplot_seiler")
plotViolin(gene = "TACSTD2", rnaData = exprs(seiler), groupVariable =  substr(str_replace_all(pData(seiler)$cstage, "[^[:alnum:]]", ""), 1, 2), yLab = "Gene expression (SCAN)", Colors = as.character(sbt6_colors[1-5]),
           Title = "TROP2", sub = "Seiler data", base_fn = fn, fileW = 7, fileH = 4, savePlot = T)

fn = paste0(dir_results_1, "/revision_1/stagep_TROP2_violinplot_sjodahl2017")
plotViolin(gene = which(fData(sjodahl17)$gene=="TACSTD2"),
           rnaData = exprs(sjodahl17), groupVariable =  sjodahl17$stage, yLab = "Gene expression (RMA normalized)", Colors = as.character(sbt6_colors[1-3]), 
           Title = "TROP2", sub = "Sjodahl 2017 data", base_fn = fn, fileW = 4, fileH = 4, savePlot = T)

fn = paste0(dir_results_1,"/revision_1/stagep_TROP2_violinplot_tcga")
plotViolin(gene = "TACSTD2",
           rnaData = tcga.rna, groupVariable =  tcga.clin$AJCC_PATHOLOGIC_TUMOR_STAGE, yLab = "Gene expression log2(RSEM+1)", Colors = as.character(sbt6_colors[1-5]),
           Title = "TROP2", sub = "TCGA Cell 2017 data", base_fn = fn, fileW = 7, fileH = 4, savePlot = T)


##############################################################################
# Neoadjuvant chemotherapy
library(GEOquery)
gse_nac <- getGEO(filename= "/Users/kaitrepka/Documents/LabWork/FengLab_Analyses/BladderCancer/TROP2_rawdata/GSE124305_series_matrix.txt")
# This is the dataset for post-NAC treatment (pre-treatment: Seiler data)

# Get TROP2 expression
sei_trop2 <- exprs(seiler)['TACSTD2',]
nac_trop2 <- exprs(gse_nac)['TACSTD2',]

# Make new matrix
expression_test <- c(as.numeric(nac_trop2), as.numeric(sei_trop2))
expression_test <- t(as.matrix(expression_test))
rownames(expression_test) <- c('TACSTD2')
group_test <- c(rep("Post-NAC", length(nac_trop2)), rep("Naive", length(sei_trop2)))

# Plot
fn = paste0(dir_results_1,"/revision_1/NAC_TROP2_violinplot")
plotViolin(gene = "TACSTD2", rnaData = expression_test, groupVariable = group_test, yLab = "Gene expression (SCAN)", Colors = colors_low_high,
           Title = "TROP2", sub = "NAC data", base_fn = fn, fileW = 4, fileH = 4, savePlot = T)

# This was unmatched samples! What happens if we use matched samples only?
gse_naive <- getGEO(filename= "/Users/kaitrepka/Documents/LabWork/FengLab_Analyses/BladderCancer/TROP2_rawdata/GSE87304_series_matrix.txt")
naive_trop2 <- exprs(gse_naive)['TACSTD2',]
library(tidyverse)
library(readxl)
nac_labels <- read_excel("TROP2_rawdata/ClinicalData_Pre_Post_NAC_with_GEO_and_WES-ID_v4.xlsx")
nac_labels <- nac_labels %>% select('case_id', 'pre_nac_geo_id_gse87304', 'post_nac_geo_id_gse124305') %>% # Select only columns of interest
  drop_na() # Drop all samples that aren't matched

# Filter the Seiler and post-NAC datasets so the samples are matched
naive_trop2_matched <- naive_trop2[names(naive_trop2) %in% nac_labels['pre_nac_geo_id_gse87304'][[1]]]
nac_trop2_matched <- nac_trop2[names(nac_trop2) %in% nac_labels['post_nac_geo_id_gse124305'][[1]]]

# Make new matrix
expression_matched <- c(as.numeric(nac_trop2_matched), as.numeric(naive_trop2_matched))
expression_matched <- t(as.matrix(expression_matched))
rownames(expression_matched) <- c('TACSTD2')
group_test_matched <- c(rep("Post-NAC", length(nac_trop2_matched)), rep("Naive", length(naive_trop2_matched)))

# Plot
fn = paste0(dir_results_1,"/revision_1/NAC_TROP2_matched_violinplot")
plotViolin(gene = "TACSTD2", rnaData = expression_matched, groupVariable = group_test_matched, yLab = "Gene expression (SCAN)", Colors = colors_low_high,
           Title = "TROP2", sub = "NAC data", base_fn = fn, fileW = 4, fileH = 4, savePlot = T)

# Paired Wilcoxon test
nac_paired_result <- wilcox.test(as.numeric(nac_trop2_matched), as.numeric(naive_trop2_matched), paired = TRUE, alternative = "two.sided")


# Get TROP2 and NECTIN4 expression values by patient ID and export
naive_jon <- t(exprs(gse_naive)[c('TACSTD2','PVRL4'),])
nac_jon <- t(exprs(gse_nac)[c('TACSTD2','PVRL4'),])
write.csv(naive_jon, 'TROP2_rawdata/patient_IDs_with_expression_for_Jon/gse87304_ids_trop2_nectin4.csv')
write.csv(nac_jon, 'TROP2_rawdata/patient_IDs_with_expression_for_Jon/gse124305_ids_trop2_nectin4.csv')

#################################################################################################
# IMVigor Data - installed from http://research-pub.gene.com/IMvigor210CoreBiologies/#downloading-the-imvigor210corebiologies-package
# Package did not compile, so I opened their "cds" datafile on a computer with R 3.4.3 and then saved to csv
IMV_cts <- read.csv("TROP2_rawdata/counts_IMVigor.csv") # Load counts matrix
IMV_fts <- read.csv("TROP2_rawdata/features_IMVigor.csv") # Load features matrix
IMV_pts <- read.csv("TROP2_rawdata/patients_IMVigor.csv") # Load patients + characteristics matrix
IMV_cts_int <- mapply('/', IMV_cts, as.numeric(colSums(IMV_cts)))*1000000 # TPM normalize
rownames(IMV_cts_int) <- rownames(IMV_cts) # Fix row names
IMV_cts <- data.frame(IMV_cts_int)
IMV_cts$gene_symbol <- IMV_fts[,'symbol'] # Add gene symbol to the counts dataframe
IMV_cts <- IMV_cts %>% drop_na() # Drop rows with NA gene symbols
IMV_cts = IMV_cts[!duplicated(IMV_cts$gene_symbol),] # Remove duplicated rows
rownames(IMV_cts) <- IMV_cts$gene_symbol # Make the rownames into the gene names
IMV_cts <- IMV_cts[ , !(colnames(IMV_cts) %in% c("X","gene_symbol"))]# Drop the extra gene symbol and "X" columns
IMV_cts <- as.matrix(IMV_cts) # Coerce to matrix
IMV_cts <- log2(IMV_cts + 1) # Transform to tpm

# Plot BCG Status
fn = paste0(dir_results_1,"/revision_1/bcg_violinplot")
plotViolin(gene = "TACSTD2", rnaData = IMV_cts, groupVariable = IMV_pts$Intravesical.BCG.administered, yLab = "Gene expression log2(TPM + 1)", Colors = colors_low_high,
           Title = "TROP2", sub = "Intravesical BCG administered", base_fn = fn, fileW = 4, fileH = 4, savePlot = T)

# Plot metastasis status
# First, define a metastasis vector to plot
metastasis_status <- IMV_pts$Met.Disease.Status
metastasis_status[metastasis_status %in% c('LN Only')] <- 'Locally advanced' # "NA" = "no metastasis"
metastasis_status[metastasis_status %in% c('Visceral', 'Liver')] <- 'Metastatic' 
# Modified dataset + metastasis status to drop na values
IMV_cts_mod <- IMV_cts[,!is.na(metastasis_status)]
metastasis_status <- metastasis_status[!is.na(metastasis_status)]
fn = paste0(dir_results_1,"/revision_1/metastasis_violinplot")
plotViolin(gene = "TACSTD2", rnaData = IMV_cts_mod, groupVariable = metastasis_status, yLab = "Gene expression log2(TPM + 1)", Colors = colors_low_high,
           Title = "TROP2", sub = "Metastasis Status", base_fn = fn, fileW = 5, fileH = 4, savePlot = T)

# Plot by tissue
fn = paste0(dir_results_1,"/revision_1/metastasis_tissue_violinplot")
tissue_status <- IMV_pts$Tissue
#tissue_status[tissue_status %in% c('bladder',"kidney","ureter")] <- 'Primary urothelial' # "NA" = "no metastasis"
#tissue_status[!(tissue_status %in% c('Bladder'))] <- 'Other (liver, lung, LN)' 
plotViolin(gene = "TACSTD2", rnaData = IMV_cts, groupVariable = tissue_status, yLab = "Gene expression log2(TPM + 1)", 
           Title = "TROP2", sub = "Tissue", base_fn = fn, fileW = 7, fileH = 4, savePlot = T)



###################################################################################################
# Survival analysis
all(colnames(tcga.rna)==rownames(tcga.clin))

fn = paste0(dir_results_1, "/revision_1/survival_tcga_os_quartiles_trop2")
plotsurvQuartiles(gex_vector = tcga.rna["TACSTD2",],
                  col_time = tcga.clin$os_time,
                  col_event = tcga.clin$os_stat, 
                  x_lab = "Years", y_lab = "Overall survival", t_lab = "TROP2, TCGA Cell 2017 data", text_x = 1, text_y = 0.2, plot_at_risk = T, 
                  x_max = 10, time_inc = 1, savePlot = T, filename = fn)

clin.sj12 <- pData(sjodahl12)
all(colnames(exprs(sjodahl12))==rownames(clin.sj12))

clin.sj12$os <- NA
clin.sj12$os[clin.sj12$`dod_event_(yes/no):ch1`=="no"] <- 0
clin.sj12$os[clin.sj12$`dod_event_(yes/no):ch1`=="yes"] <- 1

clin.sj12$os_time <- as.numeric(clin.sj12$`time_to_dod_(months):ch1`)/12

fn = paste0(dir_results_1, "/survival_sjodahl12_os_quartiles_trop2")
plotsurvQuartiles(gex_vector = exprs(sjodahl12)[which(fData(sjodahl12)$Symbol=="TACSTD2"),],
                  col_time = clin.sj12$os_time, 
                  col_event = clin.sj12$os,
                  x_lab = "Years", y_lab = "Overall survival", t_lab = "TROP2, Sjodahl 2012 data", text_x = 1, text_y = 0.2, plot_at_risk = T, 
                  x_max = 6, time_inc = 1, savePlot = T, filename = fn)





