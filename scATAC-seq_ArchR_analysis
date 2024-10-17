#####时间：20241012
#####安装包：ArchR
#####网络问题一堆包下载failed,下载到本地直接安装
install.packages("E:/Google_download/annotate_1.76.0.zip")
install.packages("E:/Google_download/DirichletMultinomial_1.40.0.zip")
install.packages("E:/Google_download/CNEr_1.34.0.zip")
install.packages("E:/Google_download/Rhdf5lib_1.20.0.zip")
install.packages("E:/Google_download/TFBSTools_1.36.0.zip")
install.packages("E:/Google_download/rhdf5filters_1.10.1.zip")
install.packages("E:/Google_download/ComplexHeatmap_2.14.0.zip")
install.packages("E:/Google_download/chromVAR_1.20.2.zip")
install.packages("E:/Google_download/motifmatchr_1.20.0.zip")
install.packages("E:/Google_download/rhdf5_2.42.1.zip")
install.packages("E:/Google_download/BSgenome.Hsapiens.UCSC.hg38_1.4.5.tar.gz")

devtools::install_github("GreenleafLab/ArchR", ref="master", repos = BiocManager::repositories())

# 载入ArchR包和依赖包
library(ArchR)
library(ggplot2)
#ArchR::installExtraPackages()
# 设置随机数种子
set.seed(1)
# 读入输入数据(/path是存放fragment文件或者bam文件的路径，如果是bam文件，请将gz$改为bam$)
inputFiles <- list.files(path = "./human_brain_3k_atac_fragments/", pattern = "gz$") 
#为输入的文件添加名字
names <- unlist(lapply(inputFiles , function(x){
  name <- unlist(strsplit(x, "[.]"))
  name <- name[1]
}))
inputFiles <- unlist(lapply(inputFiles, function(x){
  x <- file.path("F:/Data_R/test/scATAC-seq/human_brain_3k_atac/human_brain_3k_atac_fragments/", x)
}))
names(inputFiles) <- names
inputFiles
# 设置线程数和默认参考基因组
addArchRThreads(threads = 16)
addArchRGenome("hg38") #ArchR natively supports hg19, hg38, mm9, and mm10.
# 创建arrow文件
ArrowFiles <- createArrowFiles(
  inputFiles = inputFiles,
  sampleNames = names(inputFiles),
  filterTSS = 4, #Dont set this too high because you can always increase later
  filterFrags = 1000, 
  addTileMat = TRUE,
  addGeneScoreMat = TRUE
)
ArrowFiles
# 推断双细胞
doubScores <- addDoubletScores(
  input = ArrowFiles,
  k = 10, #Refers to how many cells near a "pseudo-doublet" to count.
  knnMethod = "UMAP", #Refers to the embedding to use for nearest neighbor search.
  LSIMethod = 1
)
# 创建ArchRProject对象
proj <- ArchRProject(
  ArrowFiles = ArrowFiles, 
  outputDirectory = "Human_Brain",
  copyArrows = TRUE #This is recommened so that you maintain an unaltered copy for later usage.
)
# 查看有哪些可用矩阵
getAvailableMatrices(proj)
# 过滤双细胞
proj <- filterDoublets(ArchRProj = proj)
# 降维聚类
proj <- addIterativeLSI(ArchRProj = proj, useMatrix = "TileMatrix", name = "IterativeLSI")
# 调用Seurat包的graph clustering作为默认分群方法
proj <- addClusters(input = proj, reducedDims = "IterativeLSI")
# 可视化(2D UMAP嵌入)
proj <- addUMAP(ArchRProj = proj, reducedDims = "IterativeLSI")
p1 <- plotEmbedding(ArchRProj = proj, colorBy = "cellColData", name = "Sample", embedding = "UMAP")
p2 <- plotEmbedding(ArchRProj = proj, colorBy = "cellColData", name = "Clusters", embedding = "UMAP")
ggAlignPlots(p1, p2, type = "h")
# 保存为pdf
plotPDF(p1,p2, name = "Plot-UMAP-Sample-Clusters.pdf",
        ArchRProj = proj, addDOC = FALSE, width = 5, height = 5)
# 用marker gene标记scATAC群
markerGenes  <- c(
  "CD34",  #Early Progenitor
  "GATA1", #Erythroid
  "PAX5", "MS4A1", "MME", #B-Cell Trajectory
  "CD14", "MPO", #Monocytes
  "CD3D", "CD8A"#TCells
)
p <- plotEmbedding(
  ArchRProj = proj, 
  colorBy = "GeneScoreMatrix", 
  name = markerGenes, 
  embedding = "UMAP",
  imputeWeights = getImputeWeights(proj)
)
p$CD14
# 用grid整合多张图
p2 <- lapply(p, function(x){
  x + guides() + 
    theme_ArchR(baseSize = 6.5) +
    theme(plot.margin = unit(c(0, 0, 0, 0), "cm")) +
    theme(
      axis.text.x=element_blank(), 
      axis.ticks.x=element_blank(), 
      axis.text.y=element_blank(), 
      axis.ticks.y=element_blank()
    )
})
do.call(cowplot::plot_grid, c(list(ncol = 3),p2))
# 保存为pdf
plotPDF(plotList = p, 
        name = "Plot-UMAP-Marker-Genes-W-Imputation.pdf", 
        ArchRProj = proj, 
        addDOC = FALSE, width = 5, height = 5)
# 可视化Genome Browser Tracks基因组浏览器轨道
p <- plotBrowserTrack(
  ArchRProj = proj, 
  groupBy = "Clusters", 
  geneSymbol = markerGenes, 
  upstream = 50000,
  downstream = 50000
)
grid::grid.newpage()
grid::grid.draw(p$CD14)
# 保存为pdf
plotPDF(plotList = p, 
        name = "Plot-Tracks-Marker-Genes.pdf", 
        ArchRProj = proj, 
        addDOC = FALSE, width = 5, height = 5)
# 使用shiny app更进一步可视化！
ArchRBrowser(ArchRProj = proj)
# 保存ArchRProject对象到本地
proj <- saveArchRProject(ArchRProj = proj)
# 后续重新读入ArchRProject对象
proj <- loadArchRProject(path = "./HemeTutorial")
