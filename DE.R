### 13.DE_master.R ###

# This script takes a list of dfs of different classes and types of repeat counts for
# TNB and control RNA-seq data sets and performs DE analysis:


### 0. Define variables/paths ###

# load packages needed:
library(tibble)
library(dplyr)
library(RColorBrewer)
library(RUVSeq)
library(ggplot2)
library(ggrepel)
library(preprocessCore)
library(edgeR)
library(org.Hs.eg.db)

# define starting variables:
project <- "hgsoc_repeats"
expName <- "triple_neg_breast_cancer"
Type <- "custom3"


################################################################################
### Options ###
################################################################################

################################################################################
### htseq_EdgeR_TNB_vs_normal ###
# sTypes <- c("normal", "TNB")
# sGroups <- list(c("SRR2148259", "SRR2148260", "SRR2148261"), 
#   c("SRR2148262", "SRR2148263", "SRR2148264"))
# names(sGroups) <- sTypes
# descrip <- "htseq_EdgeR_TNB_vs_normal"
################################################################################

################################################################################
# ### htseq_EdgeR_TNB_vs_normal ###
# sTypes <- c("TNB", "TNB_ant")
# sGroups <- list(c("SRR2148262", "SRR2148263", "SRR2148264"),
#                 c("SRR2148265", "SRR2148266", "SRR2148267"))
# names(sGroups) <- sTypes
# descrip <- "htseq_EdgeR_primary_TNB_vs_TNB_ant"
################################################################################

################################################################################
### htseq_EdgeR_non-normals_vs_normal ###
sTypes <- c("normal", "TNB", "TNB_ant")
sGroups <- list(c("SRR2148259", "SRR2148260", "SRR2148261"),
  c("SRR2148262", "SRR2148263", "SRR2148264"),
  c("SRR2148265", "SRR2148266", "SRR2148267"))
names(sGroups) <- sTypes
descrip <- "htseq_EdgeR_non-normals_vs_normal_with_controls"
################################################################################

################################################################################

# specify what combination of repeat genes (repeats) and other genes,
# (all, both, other) should contribute to the results:
resultTypes <- c("repeats")

# define sample group to use as control:
ctl <- "normal"

# specify what FDR and log2 fold change thresholds to use:
FDRthresh <- 0.3
FCthresh <- 0

# specify control genes to include:
posGeneIDs <- c("ENSG00000123500", "ENSG00000073282")
posGeneNames <- c("S100P", "TP63")
negGeneIDs <- c("ENSG00000163993", "ENSG00000089157")
negGeneNames <- c("beta-actin", " RPLP0")

# specify other genes to include if necessary:
#otherIDs <- c("ENSG00000130816", "ENSG00000119772", "ENSG00000088305", 
#	"ENSG00000276043", "ENSG00000138336", "ENSG00000168769", "ENSG00000187605", 
# "ENSG00000101945")

#otherSym <- c("DNMT1", "DNMT3A", "DNMT3B", "UHRF1", "TET1", "TET2", "TET3", 
#	"SUV39H1")

# specify whether to include controls:
include_ctls <- TRUE

# define directories:
homeDir <- "/Users/jamestorpy/clusterHome/"
#homeDir <- "/share/ScratchGeneral/jamtor/"
projectDir <- paste0(homeDir, "/projects/", project)
rawDir <- paste0(projectDir, 
                 "/RNA-seq/raw_files/")
resultsDir <- paste0(projectDir, "/RNA-seq/", expName, "/results")
RobjectDir <- paste0(projectDir, "/RNA-seq/",
                     expName, "/Robjects/")
newRobjectDir <- paste0(projectDir, "/RNA-seq/",
                        expName, "/Robjects/", descrip, "/")
plotDir <- paste0(resultsDir, "/R/plots/DEplots/", descrip, "/")

system(paste0("mkdir -p ", plotDir))
system(paste0("mkdir -p ", newRobjectDir))


################################################################################
### 1. Load in all counts ###
################################################################################

# define functions for this section:
counts_bind <- function(counts1, counts2) {
  # append counts1 to counts2:
  counts_all <- rbind(custom3Counts, gcCounts)
  
  # make rownames gene_id, get rid of latter column and change
  # storage mode from factor to integer:
  rownames(counts_all) <- counts_all$gene_id
  return(subset(counts_all, select=-gene_id))
}

if ( !file.exists(paste0(RobjectDir, "/", Type, "_counts.rds")) ) {
  
  writeLines("\n")
  print("Counts data frame does not exist, creating now...")
  
  custom3Counts <- readRDS(paste0(RobjectDir, "/", Type, 
                                  "_allcounts.htseq.rds"))
  gcCounts <- readRDS(paste0(RobjectDir, "/gc_allcounts.htseq.rds"))
  
  # append gcCounts to custom3Counts:
  Counts <- counts_bind(custom3Counts, gcCounts)
  
  saveRDS(Counts, paste0(RobjectDir, "/", Type, "_counts.rds"))
} else {
  
  print("Loading counts data frame...")
  Counts <- readRDS(paste0(RobjectDir, "/", Type, "_counts.rds"))
}

#  remove samples not needed:
Counts <- Counts[,colnames(Counts) %in% unlist(sGroups)]
  
# change sample names to reflect groups:
for ( n in 1:length(sTypes) ) {
  colnames(Counts)[colnames(Counts) %in% sGroups[[n]]] <- paste0(sTypes[n], "_", 
    colnames(Counts)[colnames(Counts) %in% sGroups[[n]]])
}

   
# eliminate lowly expressed genes (rows where there are less than 3 counts 
# where df > 4):
print(paste0("No. rows before filtering is: ", nrow(Counts)))
Counts <- Counts %>%
  rownames_to_column('gene_id') %>%
  dplyr::filter(rowSums(Counts > 2) >= (2)) %>%
  column_to_rownames('gene_id')
print(paste0("No. rows after  filtering: ", nrow(Counts)))
    
    
############################################################################
### 2. Perform pre-normalisation PCA and RLE plots ###
############################################################################
  
# create pre-normalised PCA plot from counts and plot:
if (ncol(Counts) > nrow(Counts)) {
  pca <- prcomp(Counts)
} else {
  pca <- princomp(Counts)	  
}
  
if (file.exists(paste0(plotDir, "/", Type, "_pcaCompsPrenormGC.pdf"))) {
  print(paste0(plotDir, "/", Type, "_pcaCompsPrenormGC.pdf already exists,
               no need to create"))
} else {
  print(paste0("Creating ", plotDir, "/", Type, "_pcaCompsPrenormGC.pdf"))
  pdf(file = paste0(plotDir, "/", Type, "_pcaCompsPrenormGC.pdf"))
  plot(pca)
  dev.off()
}
    
# change the order of columns of Counts to alphabetical order:
Counts <- Counts[,order(
  gsub(
    "_SRR.*$", "", colnames(Counts)
  )
)]
  
# define sample groups:
splt <- unlist(
  lapply(
    split(
      colnames(Counts), gsub(
        "_SRR.*$", "", colnames(Counts)
      )
    ), length
  )
)
  
for (i in 1:length(splt)) {
  if (i==1) {
    typeF <- c(rep(names(splt)[i], splt[i]))
  } else {
    typeF <- c(typeF, rep(names(splt)[i], splt[i]))
  }
}
levels(typeF) <- sTypes
  
sampleNos <- unlist(
  lapply(
    split(
      colnames(Counts), gsub(
        "\\.1", "", gsub(
          "_.*$", "", colnames(Counts)
        )
      )
    ), length
  )
)

saveRDS(sampleNos, file = paste0(newRobjectDir, "/sample_no_per_cat.rds"))
  
# convert Counts into SeqExpressionSet - elements need to be delisted and 
# changed to integers first:
Counts <- apply(Counts, 2, unlist)
storage.mode(Counts) <- "integer"
set <- newSeqExpressionSet(Counts, phenoData = data.frame(typeF, 
                                                          row.names=colnames(Counts)))

# create pre-norm RLE plot:
if (file.exists(paste0(plotDir, "/", Type, "_RLEPrenormGC.pdf"))) {
  print(paste0(plotDir, "/", Type, "_RLEPrenormGC.pdf already exists, 
             no need to create"))
} else {
  print(paste0("Creating ", plotDir, "/", Type, "_RLEPrenormGC.pdf"))
  par(mar=c(1,1,1,1))
  pdf(file = paste0(plotDir, "/", Type, "_RLEPrenormGC.pdf"))
  plotRLE(set)
  dev.off()
}

# create RUVseq pre-norm PCA:
if (file.exists(paste0(plotDir, "/", Type, "_pcaPrenormGC.pdf"))) {
  print(paste0(plotDir, "/", Type, "_pcaPrenormGC.pdf already exists, 
             no need to create"))
} else {
  print(paste0("Creating ", plotDir, "/", Type, "_pcaPrenormGC.pdf"))
  pdf(file = paste0(plotDir, "/", Type, "_pcaPrenormGC.pdf"), height = 10, 
      width = 12)
  plotPCA(set, cex=0.7)
  dev.off()
}
    
##############################################################################
### 3. perform normalisation and DE on counts:
##############################################################################

# perform between lane full normalisation:
y <- DGEList(counts = Counts, group = typeF)

# normalise for library size:
y <- calcNormFactors(y)

# create an MDS plot to show relative similarities of the samples and save to Dir:
if (file.exists(paste0(plotDir, "/edger_MDS.pdf"))) {
  paste0(plotDir, "/edger_MDS.pdf already exists")
  pdf(paste0(plotDir, "/edger_MDS.pdf"),width=16,height=12)
  plotMDS(y)
} else {
  paste0("Generating ", plotDir, "/edger_MDS.pdf")
  pdf(paste0(plotDir, "/edger_MDS.pdf"),width=16,height=12)
  plotMDS(y)
  dev.off()
}

# design matrix labelling all sample types:
design <- model.matrix(~0+typeF)

# estimate dispersion:
disp <- estimateDisp(y, design=design)

# adjust values using dispersion:
fit <- glmFit(disp, design=design, robust=TRUE)

saveRDS(fit, file=paste0(newRobjectDir, "/", Type, "DEfit.rds"))
#save.image(paste0(newRobjectDir, "/", Type, "DEdone.rds"))


# determine which column has control:
ctlInd <- grep(ctl, colnames(design))
con <- c(rep(0, (ctlInd - 1) ), -1, rep(0, (ncol(design) - ctlInd)))

# put sTypes in alphabetical order:
sTypes <- sTypes[order(sTypes)]

# check parameters and give the user the option to continue or not:
writeLines("\n")
print("Contrast is: ")
print(con)
print("Column names of design are: ")
print(colnames(design))
print("Design matrix is: ")
print(design)

for (i in 1:ncol(design)) {
  print(i)
  if ( i!=ctlInd ) {
    comp <- paste0(sTypes[i], "_vs_", ctl)
    
    # perform likelihood ratio test:
    con[i] <- 1
    lrt <- glmLRT(fit, contrast = con)
    
    # determine the top DE genes:
    topTags(lrt)
    
    if (file.exists(paste0(newRobjectDir, sTypes[i], "_vs_", ctl, "_lrt.rds"))) {
      print(paste0(newRobjectDir, sTypes[i], "_vs_", ctl, "_lrt.rds already exists, 
                   no need to create"))
    } else {
      print(paste0("Creating ", newRobjectDir, sTypes[i], "_vs_", ctl, "_lrt.rds"))
      saveRDS(lrt, file = paste0(newRobjectDir, sTypes[i], "_vs_", ctl, "_lrt.rds"))
    }
    
    # fetch summary of differentially expressed genes (those  with FDR =< 0.05:
    DEs <- summary(result <- decideTestsDGE((lrt)))
    
    # fetch all gene DE info, 
    allGenes <- as.data.frame(topTags(lrt, n=Inf))
    
    # annotate allGenes with entrez ids and symbols in separate columns:
    egENSEMBL <- toTable(org.Hs.egENSEMBL)
    egSYMBOL <- toTable(org.Hs.egSYMBOL)
    
    allGenes$gene_id <- egENSEMBL$gene_id[match(rownames(allGenes), egENSEMBL$ensembl_id)]
    allGenes$symbol <- egSYMBOL$symbol[match(allGenes$gene_id, egSYMBOL$gene_id)]

    if (!(ctlInd==1)) {
      if (i==1) {
        allGenesList <- list(allGenes)
      } else {
        allGenesList[[i]] <- allGenes
      }
      
    } else {
      if (i==2) {
        allGenesList <- list(allGenes)
      } else {
        allGenesList[[i]] <- allGenes
      }
    }
    
    
    ##############################################################################
    ### 4. Create DE data frames for repeats:
    ##############################################################################
    
    # define repeat and sig DE repeat dfs:
    repGenes <- allGenes[grep("ENS",  rownames(allGenes), invert = T),]
    print(repGenes)
    
    if ( is.na(FCthresh) ) {
      sigGenes <- filter(repGenes, FDR < FDRthresh)
      repGenes$threshold <- as.factor(repGenes$FDR < FDRthresh)
    } else if ( is.na(FDRthresh) ) {
      sigGenes <- repGenes[
        (repGenes$logFC > FCthresh)|(repGenes$logFC < -(FCthresh)), ]
      repGenes$threshold <- as.factor( 
        (repGenes$logFC > FCthresh)|(repGenes$logFC < -(FCthresh))
      )
    } else {
      sigGenes <- filter(
        repGenes, 
        (FDR < FDRthresh & logFC < -(FCthresh))|(
          FDR < FDRthresh & logFC > FCthresh
        )
      )
      repGenes$threshold <- as.factor(
        (repGenes$FDR < FDRthresh & repGenes$logFC < -(FCthresh))|(
          repGenes$FDR <  FDRthresh & repGenes$logFC > FCthresh
        )
        )
    }
    
    repGenes$type <- "repeat"
    sig_rep <- subset(repGenes, threshold == T)
    
    # include the control genes for labelling:
    # for (j in 1:length(posGeneIDs)) {
    #   if (j==1) {
    #     posGenes <- allGenes[ posGeneIDs[j],]
    #   } else {
    #     posGenes <- rbind(posGenes, allGenes[posGeneIDs[j],])
    #   }
    # }
    # rownames(posGenes) <- posGeneNames
    # 
    # for (j in 1:length(negGeneIDs)) {
    #   if (j==1) {
    #     negGenes <- allGenes[ negGeneIDs[j],]
    #   } else {
    #     negGenes <- rbind(negGenes,   allGenes[negGeneIDs[j],])
    #   }
    # }
    # rownames(negGenes) <- negGeneNames
    # 
    # # set default threshold statuses for control genes:
    # posGenes$threshold <- "POSITIVE"
    # if (nrow(posGenes[posGenes$FDR< FDRthresh,])>0) {
    #   posGenes[posGenes$FDR<  FDRthresh,]$threshold <- "POSSIG"
    # }
    # 
    # negGenes$threshold = "NEGATIVE"
    # if (nrow(negGenes[negGenes$FDR< FDRthresh,])>0) {
    #   negGenes[negGenes$FDR<  FDRthresh,]$threshold <-  "NEGSIG"
    # }
    
    if (!(ctlInd==1)) {
      if (i==1) {
        allReps <- list(repGenes)
      } else {
        allReps[[i]] <- repGenes
      }
      
      if (i==1) {
        sigReps <- list(sig_rep)
      } else {
        sigReps[[i]] <- sig_rep
      }
    } else {
      if (i==2) {
        allReps <- list(repGenes)
      } else {
        allReps[[i]] <- repGenes
      }
      
      if (i==2) {
        sigReps <- list(sig_rep)
      } else {
        sigReps[[i]] <- sig_rep
      }
    }
    
    
    ##############################################################################
    ### 5. Create DE data frames for gencode genes:
    ##############################################################################
    
    if (length(FCthresh) == 0) {
      sigGenes <- filter(allGenes, FDR < FDRthresh)
      allGenes$threshold <- as.factor(allGenes$FDR < FDRthresh)
    } else {
      sigGenes <- filter(allGenes, (FDR < FDRthresh & logFC < -(FCthresh))|(FDR < FDRthresh & logFC > FCthresh))
      allGenes$threshold <- as.factor((allGenes$FDR < FDRthresh & allGenes$logFC < -(FCthresh))|(allGenes$FDR <  FDRthresh & allGenes$logFC > FCthresh))
    }
    
    sig_gc <- subset(allGenes, threshold == T)
    
    if (!(ctlInd==1)) {
      if (i==1) {
        sig_gc_GenesList <- list(sig_gc)
      } else {
        sig_gc_GenesList[[i]] <- sig_gc
      }
    } else {
      if (i==2) {
        sig_gc_GenesList <- list(sig_gc)
      } else {
        sig_gc_GenesList[[i]] <- sig_gc
      }
    }
    
    # include the control genes for labelling:
    for (j in 1:length(posGeneIDs)) {
      if (j==1) {
        posGenes <- allGenes[ posGeneIDs[j],]
      } else {
        posGenes <- rbind(posGenes, allGenes[posGeneIDs[j],])
      }
    }
    rownames(posGenes) <- posGeneNames
    
    for (j in 1:length(negGeneIDs)) {
      if (j==1) {
        negGenes <- allGenes[ negGeneIDs[j],]
      } else {
        negGenes <- rbind(negGenes,   allGenes[negGeneIDs[j],])
      }
    }
    rownames(negGenes) <- negGeneNames
    
    # set default threshold statuses for control genes:
    posGenes$threshold <- "POSITIVE"
    if (nrow(posGenes[posGenes$FDR< FDRthresh,])>0) {
      posGenes[posGenes$FDR<  FDRthresh,]$threshold <- "POSSIG"
    }
    
    negGenes$threshold = "NEGATIVE"
    if (nrow(negGenes[negGenes$FDR< FDRthresh,])>0) {
      negGenes[negGenes$FDR<  FDRthresh,]$threshold <-  "NEGSIG"
    }
    
    ctlGenes <- rbind(posGenes, negGenes)
    
    # add 'type' identifier column:
    ctlGenes$type <- "ctl"
    
    ##############################################################################
    ### 6. Create volcano plots:
    ##############################################################################
    
    if ("repeats" %in% resultTypes) {
      
      if (include_ctls) {
        
        lab <- rbind(sig_rep, ctlGenes)
        lab$genes <- rownames(lab)
        repGenes <- rbind(repGenes, ctlGenes)
        
        # combine 'threshold' and 'type' columns:
        repGenes$type_thresh <- paste0(repGenes$type, "_", repGenes$threshold)
        lab$type_thresh <- paste0(lab$type, "_", lab$threshold)
        
      } else {
        
        lab <- sig_rep
        lab$genes <- rownames(lab)
        
        # combine 'threshold' and 'type' columns:
        repGenes$type_thresh <- paste0(repGenes$type, "_", repGenes$threshold)
        repGenes$type_thresh <- factor(repGenes$type_thresh, levels = c("ctl_NEGATIVE", 
          "ctl_POSITIVE", "ctl_POSSIG", "repeat_FALSE", "repeat_TRUE"))
        lab$type_thresh <- paste0(lab$type, "_", lab$threshold)
        lab$type_thresh <- factor(lab$type_thresh, levels = c("ctl_NEGATIVE", 
          "ctl_POSITIVE", "ctl_POSSIG", "repeat_FALSE", "repeat_TRUE"))
      }
      
      
      # plot on volcano plot:
      p <- ggplot(data=repGenes, aes( x=logFC, y=-log10(FDR), color=type_thresh))
      p <- p + geom_point(data=repGenes)
      p <- p + geom_text_repel(data=lab, aes(label=genes))
      p <- p + theme(legend.position =  "none")
      p <- p + labs(x="log2 fold change vs FT control", y="-log10 FDR")
      # key for colours = c("neg_ctls", "pos_ctls", "neg_gc", "pos_gc")
      #          p <- p + scale_colour_manual(values = c("#114477", "firebrick4", 
      #            "dodgerblue1", "firebrick1"))
      p <- p + scale_colour_manual(values = c("#205C89", "#C98912", "#C98912", "#A1D495",
                                              "#169132"))
      p <- p +  xlim(c(-4, 4))
      if (length(FCthresh) == 0) {
        if (file.exists(paste0(plotDir,   "/", Type,  "_volcano_FDR_",   FDRthresh, "_", comp, "_reps.pdf"))) {
          print(paste0(plotDir, "/",  Type,  "_volcano_FDR_",   FDRthresh, "_", comp, "_reps.pdf"))
          p
        } else {
          print(paste0("Creating  ",plotDir, "/", Type,    "_volcano_FDR_", FDRthresh, "_", comp, "_reps.pdf"))
          pdf(file = paste0(plotDir, "/",   Type,  "_volcano_FDR_",  FDRthresh, "_", comp, "_reps.pdf"))
          print(p)
          dev.off()
        }
      } else {
        if (file.exists(paste0(plotDir, "/",  Type,  "_volcano_FDR",   FDRthresh, "_FC", FCthresh, "_", comp, "_reps.pdf"))) {
          print(paste0(plotDir, "/",  Type,  "_volcano_FDR",   FDRthresh, "_FC", FCthresh, "_", comp, "_reps.pdf already exists"))
          p
        } else {
          print(paste0("Creating  ", plotDir, "/",  Type,  "_volcano_FDR", FDRthresh, "_FC", FCthresh, "_", comp, "_reps.pdf"))
          pdf(file = paste0(plotDir, "/",  Type,  "_volcano_FDR",   FDRthresh, "_FC", FCthresh, "_", comp, "_reps.pdf"))
          print(p)
          dev.off()
        }
      } 
    } 
    con <- c(rep(0, (ctlInd - 1) ), -1, rep(0, (ncol(design) - ctlInd)))
  }
}

# remove the NULL list dfs created when avoiding clt vs ctl:
if (length(sTypes)>2) {
  allReps <- allReps[-ctlInd]
  sigReps <- sigReps[-ctlInd]
  # name the list elements:
  names(allReps) <- paste0(sTypes[-ctlInd], "_vs_", ctl)
  names(sigReps) <- paste0(sTypes[-ctlInd], "_vs_", ctl)
}