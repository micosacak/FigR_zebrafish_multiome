{  
  args <- commandArgs(trailingOnly = TRUE)
  
  trtid = as.numeric(gsub("value=", "", args[1]))
  ncpus = as.numeric(gsub("ncpus=", "", args[2]))
  nrams = as.numeric(gsub("nrams=", "", args[3]))
  nidxs = gsub("nidxs=","",args[4])
  mxram = as.numeric(gsub("mxram=", "", args[5]))
  runAll = as.numeric(gsub("runall=", "", args[6]))
    
  if(runAll == 0){
      runAll = TRUE
  }else{
      runAll = FALSE
  }

  library(future)
  options(future.globals.maxSize = mxram*1024^3)

  print(trtid)

  source("/group/crtd_becker/Data/00_biocluster4/MIC_scRNA_Seq/77_CommonFunctionS5.R")
  ldFnx(nCPUs = nCPUs, mRAM = mRAM)
  # change the directory

  setRAM(nrams)
  nCPUs = ncpus

  print(list(
    "trtid" = trtid,
    "ncpus" = ncpus,
    "nidxs" = nidxs,
    "nrams" = nrams,
    "mxram" = mxram,
    "runAll" = runAll
    ))

  setDir(theFolder = "LAB5255newCompletedMultiERG")

  the_prefix = "20250909_normRNA_erg_LAB5255newCompletedMultiERG_main3allcellsA"
  crtD(the_prefix)


  if(!file.exists("pbmc3_annotations_small.rds")){
    seuobj = readRDS("pbmc3_annotations.rds")
    seuobj@assays$GAC = NULLx
    seuobj@assays$SCT = NULL
    seuobj@assays$ATAC = NULL
    seuobj@assays$GEX$scale.data = NULL
    saveRDS(seuobj, file = "pbmc3_annotations_small.rds")
  }else{
    print("Loading ...")
    seuobj = readRDS("pbmc3_annotations_small.rds")
  }
  
  nw.idx = seuobj$main.cells
  # edit the cluster names and pass it to seuobj$the_clusters
  
  seuobj$the_clusters = nw.idx
  Idents(seuobj) = seuobj$the_clusters

  the_samples = gsub("_R[0-9]+","",seuobj$samples)
  seuobj$the_samples = the_samples
  all_samples = unique(seuobj$the_samples)

  trt = all_samples[trtid]
  all_samples = c(all_samples[trtid])
  all_samples = sort(all_samples)
  print(all_samples)

  the_uniqs = unique(seuobj$the_samples)
  the_samples = the_uniqs
  names(the_samples) = the_uniqs

  if(nidxs == ""){
    nidxs = 1
    Idents(seuobj) = rep("allcells", ncol(seuobj))
  }else{
    if(length(unlist(strsplit(nidxs, split = ":"))) == 2){
      tmp_split = unlist(strsplit(nidxs, split = ":"))
      nidxs = seq(as.integer(tmp_split[1]),as.integer(tmp_split[2]))
      Idents(seuobj) = seuobj$the_clusters
    }else{
      Idents(seuobj) = seuobj$the_clusters
      nidxs = as.integer(nidxs)
    }
  }
  seuobj$ident2use = Idents(seuobj)
  #the_clusters = unique(Idents(seuobj))[nidxs]
  the_clusters = levels(seuobj$ident2use)[nidxs]
  print(the_clusters)

  
  ### SUBSET all_samples
  curDir = getwd()
  seuobj$samples = seuobj$the_samples

  if(runAll){
    seuobj$samples = rep("ALL",ncol(seuobj))
    seuobj$ident2use = rep("ALL", ncol(seuobj))
    trt = "ALL"
    the_clusters = c("ALL")
  }

  qwq = seuobj

  # edit for the correct path for gtf file.
  # here we used the gtf file used for multiome alignment.
  
  gtf <- rtracklayer::import("/group/crtd_becker/Data/00_biocluster4/MIC_scRNA_Seq/001_Alignments/Genomes/dre/22_Multiome/dre_GRCz11_atacARC_custom/genes/genes.gtf.gz")
  genes <- gtf[gtf$type == "gene"]
  grcZ11TSSRanges = genes
  
  #some default values
  ndims = 30
  minNgeneCutoff = 1
  nBG = 100
  nWindowsSize = 50000
  
  # for loop here
  #for(trt in all_samples){
  if(!runAll){
    seuobj_trt = subset(qwq, samples %in% the_uniqs[grep(trt,the_samples)])
  }
  newDir = paste0(the_prefix,"/",trt)
  crtDfull(newDir)
  
  expected_files = c("seuobj_dorcGenes.rds","seuobj_cisCorr.filt.rds",
      "seuobj_numDorcs.rds","seuobj_dorcMat.rds",
      "seuobj_dorcMat.s.rds","seuobj_RNAmat.s.rds","seuobj_dorcMat.rds") 

for(the_cluster in the_clusters){
    if(!runAll){
      seuobj = subset(seuobj_trt, ident2use %in% the_cluster)
    }
    if(ncol(seuobj) < 9){
      print("Skipping ..")
      next
    }
    print(table(seuobj$samples))
    newDir = paste0(the_prefix,"/",trt,"/",gsub("/","_",the_cluster))
    #
    getwd()
    setwd(curDir)
    
    #
    crtD(newDir, setDir = T)
    atac_se = SummarizedExperiment(assays = list(counts = GetAssayData(seuobj, assay = "peaks", layer = "data")), rowRanges = granges(seuobj[["peaks"]]))
    rna_mat = GetAssayData(seuobj, assay = "GEX", layer = "data")
    if(ncol(rna_mat) < 30) ndims = ncol(rna_mat)-1
    cellkNN = get.knn(seuobj@reductions$umap.gex@cell.embeddings, k = ndims)$nn.index
    rownames(cellkNN) = colnames(seuobj)
    colnames(cellkNN) = colnames(seuobj@reductions$umap.gex@cell.embeddings)[1:ndims]
    # Visualize on pre-computed UMAP
    umap.d <- as.data.frame(seuobj@reductions$umap.gex@cell.embeddings)
    colnames(umap.d) = c("UMAP1","UMAP2")
    saveRDS(umap.d, file = "seuobj_umap.d.rds")
    
    #
    if(!file.exists("seuobj_footprints_cisCor.rds")){
      cisCorr <- runGenePeakcorr_dre(ATAC.se = atac_se, 
                                     RNAmat = rna_mat, 
                                     genome = "GRCz11", 
                                     nCores = nCPUs, 
                                     p.cut = NULL, 
                                     normalizeATACmat = FALSE, 
                                     keepPosCorOnly = FALSE,
                                     keepMultiMappingPeaks = TRUE,
                                     n_bg = nBG, 
                                     windowPadSize = nWindowsSize)
      saveRDS(cisCorr, file = "seuobj_footprints_cisCor.rds")
    }else{
      print("Loading ...  seuobj_footprints_cisCor.rds")
      cisCorr = readRDS("seuobj_footprints_cisCor.rds")
    }
    
    if(!all(file.exists(expected_files))){
      oo(8,8)
      cisCorr.filt <- cisCorr %>% filter(pvalZ <= 0.1)
      dorcGenes <- dorcJPlot(dorcTab = cisCorr.filt, cutoff = minNgeneCutoff, labelTop = 20, returnGeneList = TRUE, force=2)
      saveRDS(dorcGenes, file = "seuobj_dorcGenes.rds")
      saveRDS(cisCorr.filt, file = "seuobj_cisCorr.filt.rds")
      # Unfiltered
      numDorcs <- cisCorr.filt %>% group_by(Gene) %>% tally() %>% arrange(desc(n))
      saveRDS(numDorcs, file = "seuobj_numDorcs.rds")
      dorcMat <- getDORCScores(ATAC.se = atac_se, dorcTab = cisCorr.filt, geneList = dorcGenes, nCores = nCPUs)
      saveRDS(dorcMat, file = "seuobj_dorcMat.rds")
      # Smooth dorc scores using cell KNNs (k=30)
      dorcMat.s <- smoothScoresNN(NNmat = cellkNN[,1:ndims], mat = dorcMat, nCores = nCPUs)
      # Smooth RNA using cell KNNs
      # This takes longer since it's all genes
      RNAmat.s <- smoothScoresNN(NNmat = cellkNN[,1:ndims], mat = rna_mat, nCores = nCPUs)
      saveRDS(dorcMat.s, file = "seuobj_dorcMat.s.rds")
      saveRDS(RNAmat.s, file = "seuobj_RNAmat.s.rds")
      saveRDS(dorcMat, file = "seuobj_dorcMat.rds")
    }else{
      print("Loading ...    seuobj_dorcGenes.rds")
      dorcGenes = readRDS(file = "seuobj_dorcGenes.rds")
      cisCorr.filt = readRDS(file = "seuobj_cisCorr.filt.rds")
      numDorcs = readRDS(file = "seuobj_numDorcs.rds")
      dorcMat = readRDS(file = "seuobj_dorcMat.rds")
      dorcMat.s = readRDS(file = "seuobj_dorcMat.s.rds")
      RNAmat.s = readRDS(file = "seuobj_RNAmat.s.rds")
      dorcMat = readRDS(file = "seuobj_dorcMat.rds") 
    }
    setwd(curDir)
  }
  #}
  
  # for loop here
  #for(trt in all_samples){
  if(!runAll){
    seuobj_trt = subset(qwq, samples %in% the_uniqs[grep(trt,the_samples)])
  }
  newDir = paste0(the_prefix,"/",trt)
  crtDfull(newDir)
  
  for(the_cluster in the_clusters){
    if(!runAll){
      seuobj = subset(seuobj_trt, the_clusters %in% the_cluster)
    }
    if(ncol(seuobj) < 9){
      print("Skipping ..")
      next
    }

    print(table(seuobj$samples))
    newDir = paste0(the_prefix,"/",trt,"/",gsub("/","_",the_cluster))
    #
    getwd()
    setwd(curDir)
    
    #
    crtD(newDir, setDir = T)
    
    #if(!all(file.exists(expected_files))){
      print("Loading ... for TFenrich analyses")
      dorcGenes = readRDS(file = "seuobj_dorcGenes.rds")
      cisCorr.filt = readRDS(file = "seuobj_cisCorr.filt.rds")
      numDorcs = readRDS(file = "seuobj_numDorcs.rds")
      dorcMat = readRDS(file = "seuobj_dorcMat.rds")
      dorcMat.s = readRDS(file = "seuobj_dorcMat.s.rds")
      RNAmat.s = readRDS(file = "seuobj_RNAmat.s.rds")
      dorcMat = readRDS(file = "seuobj_dorcMat.rds") 

#run_FigGRN_dre_final this fucntion was edited in such a way to be used for zebrafish!
# TO DO: the motifs used are from JASPAR! However, figR requires motif names not motif IDs.
# edit the code in such a way it can use motif IDs.
    
      if(!file.exists("seuobj_TFenrich.d.rds")){
        TFenrich.d = run_FigGRN_dre_final(ATAC.se = atac_se, dorcMat = dorcMat.s,  
                                          rnaMat = RNAmat.s, dorcTab = cisCorr.filt,  
                                          #DORC.knn = cellkNN, 
                                          dorcGenes = NULL, dorcK = 30, n_bg = 100, 
                                          genome = "GRCz11", nCores = nCPUs)
        saveRDS(TFenrich.d, file = "seuobj_TFenrich.d.rds")
      }else{
        print("Loading ...")
        TFenrich.d = readRDS(file = "seuobj_TFenrich.d.rds")
      }
    #}
    
    setwd(curDir)
    #}
  }

  quit(save = "no")
}
