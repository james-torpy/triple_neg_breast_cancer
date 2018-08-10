### 10.prepareCounts.R ###

# breaks gcCounts up into dataframe for each class and
# saves each as RData objects #

# define starting variables:
project <- "hgsoc_repeats"
expName <- "triple_neg_breast_cancer"

# define directories:
#homeDir <- "/share/ScratchGeneral/jamtor/"
homeDir <- "/Users/jamestorpy/clusterHome"
projectDir <- paste0(homeDir, "/projects/", 
  project, "/RNA-seq/")
resultsDir <- paste0(projectDir, "/", expName, "/results/")
RobjectDir <- paste0(projectDir, "/", expName, "/Robjects/")
inDir <- paste0(resultsDir, "/htseq/")
refDir <- paste0(projectDir, 
  "/", expName, "/refs/")
starGC_dir <- paste0(resultsDir, "/star/GC/", expName)
starRibo_dir <- paste0(resultsDir, "/star/ribo/", 
  expName)

system(paste0("mkdir -p ", RobjectDir))


### 1. Check for and discard bad inputs ###

# fetch counts file ids:
compFiles <- unique(gsub("\\..*", "", list.files(inDir)))


# check gc files have the same number of lines and put
# these in a vector:
for (i in 1:length(compFiles)) {
  if (i==1) {
    lineNo1 <- strsplit(system(paste0("wc -l ", inDir, 
      "/", compFiles[i], ".gc.htseq.txt"), intern=T), 
    " ")[[1]][4]
    print(paste0(compFiles[i], " has ", lineNo1, 
      " lines"))
    print(paste0("Appending ", compFiles[i], 
      " to goodFilesGC"))
    goodFilesGC <- c(compFiles[i])
    badFilesGC <- c()
  } else {
    lineNo <- strsplit(system(paste0("wc -l ", inDir, 
      "/", compFiles[i], 
      ".gc.htseq.txt"), intern=T), " ")[[1]][4]
    print(paste0(compFiles[i], " has ", lineNo, 
      " lines"))
    if (lineNo == lineNo1) {
      print(paste0("Appending ", compFiles[i], 
        " to goodFilesGC"))
      goodFilesGC <- append(goodFilesGC, 
        compFiles[i])
    } else {
      print(paste0("Appending ", compFiles[i], 
        " to \
        badFilesGC"))
      badFilesGC <- append(badFilesGC, 
        compFiles[i])
    }
  }
}

# check custom3files corresponding with ids of goodFilesGC have the same number of lines and put this into a vector:
for (i in 1:length(goodFilesGC)) {
  if (i==1) {
    lineNo1 <- strsplit(system(paste0("wc -l ", inDir, 
      "/", goodFilesGC[i], 
      ".custom3.htseq.txt"), intern=T), " ")[[1]][6]
    print(paste0(goodFilesGC[i], " has ", lineNo1, 
      " lines"))
    print(paste0("Appending ", goodFilesGC[i], 
      " to goodFilesCustom3"))
    goodFilesCustom3 <- c(goodFilesGC[i])
    badFilesCustom3 <- c()
  } else {
    lineNo <- strsplit(system(paste0("wc -l ", inDir, 
      "/", goodFilesGC[i], 
      ".custom3.htseq.txt"), intern=T), " ")[[1]][6]
    print(paste0(goodFilesGC[i], " has ", lineNo, 
      " lines"))
    if (lineNo == lineNo1) {
      print(paste0("Appending ", goodFilesGC[i], 
        " to goodFilesCustom3"))
      goodFilesCustom3 <- append(goodFilesCustom3, 
        goodFilesGC[i])
    } else {
      print(paste0("Appending ", goodFilesGC[i], 
        " to badFilesCustom3"))
      badFilesCustom3 <- append(badFilesCustom3, 
        goodFilesGC[i])
    }
  }
}

goodFiles <- goodFilesCustom3

saveRDS(goodFiles, paste0(RobjectDir, "/goodFiles.rds"))


### 2. Load in GC inputs ###

# load in gc_countFiles into a df:
for (i in 1:length(goodFiles)) {
  if (i==1) {
    print(goodFiles[i])
    CountsGC <- data.frame(read.table(file=paste0(inDir, "/", goodFiles[i], ".gc.htseq.txt")))
  } else {
    print(goodFiles[i])
    CountsGC <- cbind(CountsGC, read.table(file=paste0(inDir, "/", goodFiles[i], ".gc.htseq.txt"))[,2])
  }
}
colnames(CountsGC) <- c("gene_id", goodFiles)

######
saveRDS(CountsGC, paste0(RobjectDir, "gcCounts_temp.rds"))
######

# aggregate multiple types of same gene:
CountsGC$gene_id <- gsub("\\..*$", "", CountsGC$gene_id)
CountsGC <- aggregate(.~gene_id, CountsGC, mean)
# remove duplicate rows:
CountsGC <- unique(CountsGC)
# remove specs lines:
CountsGC <- CountsGC[grep("__", CountsGC$gene_id, invert=T),]

# save the counts:
if (!file.exists(paste0(RobjectDir, "/gc_allcounts.htseq.rds"))) {
  saveRDS(CountsGC, file=paste0(RobjectDir, "/gc_allcounts.htseq.rds"))
}


### 3. Load in 'custom3' inputs ###

# load in custom3_countFiles into a df:
for (i in 1:length(goodFiles)) {
  if (i==1) {
    print(goodFiles[i])
    CountsCustom3 <- data.frame(read.table(file=paste0(inDir, "/", goodFiles[i], ".custom3.htseq.txt")))
  } else {
    print(goodFiles[i])
    CountsCustom3 <- cbind(CountsCustom3, read.table(file=paste0(inDir, "/", goodFiles[i], ".custom3.htseq.txt"))[,2])
  }
}
colnames(CountsCustom3) <- c("gene_id", goodFiles)

# remove duplicate rows:
CountsCustom3 <- unique(CountsCustom3)
# remove specs lines:
CountsCustom3 <- CountsCustom3[grep("__", CountsCustom3$gene_id, invert=T),]

# save the counts:
if (!file.exists(paste0(RobjectDir, "/custom3_allcounts.htseq.rds"))) {
  saveRDS(CountsCustom3, file=paste0(RobjectDir, "/custom3_allcounts.htseq.rds"))
}
