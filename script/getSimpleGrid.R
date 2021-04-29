#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)
seqscope1st=args[1]
DGEdir=args[2]
spatial=args[3]
tiles=args[4]
nrow=as.numeric(args[5])
ncol=as.numeric(args[6])
sidesize=as.numeric(args[7])
outpath=args[8]
collapsePath=args[9]
tiles=as.numeric(unlist(strsplit(tiles,',')))
####################################################################################################3
#' This function collapse tiles into small grids and create Seurat objects
#' @param tile_df Dataframe including coordinates
#' @param i integer of tile number
#' @param binx binning size for x coordinates
#' @param biny binning size for y coordinates
#' @param m_tile sparse count matrix
#' @export
collapseTiles=function(tile_df,i,binx,biny,m_tile)
{

  
  tile_df_i=tile_df[tile_df$tile_miseq==i,]
  miny = min(tile_df_i$y_miseq)
  maxy = max(tile_df_i$y_miseq)
  minx = min(tile_df_i$x_miseq)
  maxx = max(tile_df_i$x_miseq)
  xlim= c(min(tile_df_i$x_miseq),max(tile_df_i$x_miseq))
  ylim = c(min(tile_df_i$y_miseq),max(tile_df_i$y_miseq))
  grd = make.grid(tile_df_i$x_miseq,tile_df_i$y_miseq,tile_df_i$UMI, binx,biny, xlim, ylim)
  grd=t(grd)
  #colom
  fn1=paste0('Temp_CollapsedHDMIsIndLength','.csv')
  fn2=paste0('Temp_CollapsedHDMIsInd','.txt')
  if (file.exists(fn1)) {
    #Delete file if it exists
    file.remove(fn1)
  }

  if (file.exists(fn2)) {
    #Delete file if it exists
    file.remove(fn2)
  }
  grd_re = make.grid(tile_df_i$x_miseq,tile_df_i$y_miseq,tile_df_i$tileHDMIind, binx, biny, xlim, ylim,function(x) {write(length(x), file=fn1,append = T)})
  grd_re = make.grid(tile_df_i$x_miseq,tile_df_i$y_miseq,tile_df_i$tileHDMIind, binx, biny, xlim, ylim,function(x) {cat(x,file=fn2,append=TRUE,sep='\n')})

  collapseLen = read.csv(fn1,header=F)
  collapseInd = read.table(fn2,header=F)


  collapseLen = cbind(collapseLen,cumsum(collapseLen$V1))
  colnames(collapseLen) =c("len","end")
  interv = c(0,collapseLen$end)


  #sourceCpp(collapsePath)
  print('Start Simple Square Gridding!')
  tic();collapseM = collapse(m_tile,collapseInd$V1,interv);toc()
  rownames(collapseM) = rownames(m_tile)
  colnames(collapseM) = paste0("Collapse_",i,'_',1:(length(interv)-1))
  sparse.gbm <- Matrix(collapseM , sparse = T )
  writeMM(obj = sparse.gbm, file=paste0('tile',i,"collapsedMatrix.mtx"))
  write.csv(rownames(collapseM),paste0('tile',i,'collapsedGenes.csv'))
  write.csv(colnames(collapseM),paste0('tile',i,'collapsedBarcodes.csv'))

  grd=t(grd)
  pos = which(!is.na(grd), TRUE)
  pos_coor = t(sapply(1:(dim(pos)[1]),function(x) {c(as.numeric(rownames(grd)[as.numeric(pos[x,1])]),as.numeric(colnames(grd)[as.numeric(pos[x,2])]) )})) #here the grd is transposed
  colnames(pos_coor) =c("X","Y")
  coord.df = data.frame("Y"=pos_coor[,2], "X"=pos_coor[,1],"tile"=i, stringsAsFactors=FALSE)
  #write.csv(coord.df,paste0('tile',i,'coor_df_stratgy1.csv'))


  obj1 = CreateSeuratObject(counts=collapseM,assay='Spatial')
  #obj1$status = "Original"
  obj1@meta.data$tile = coord.df$tile
  obj1@meta.data$X = coord.df$X
  obj1@meta.data$Y = coord.df$Y
  obj1@meta.data$tile = i


  if(i==tiles[1])
  {
    obj=obj1
  }
  else
  {
    obj = merge(obj,obj1)

  }

  return (obj)
}




####################################################################################################3
#' This function grids the data with user-defined binning size and collapsed HDMIs within each grid
#' @param seqscope1st Data flatform,
#' @param DGEdir folder that stores barcodes.tsv, features.tsv and matrix.mtx
#' @param spatial txt file stores spatial informaiton with four columns: 'HDMI','tile_miseq','x_miseq','y_miseq'
#' @param tiles a vector of tiles that the user is interested in collapsing
#' @param nrow an integer of how many rows to organize the tiles
#' @param ncol an integer of how many cols to organize the tiles
#' @param sidesize side size of the square grid 300 represents 10um
#' @param outpath path to store the output, you need to make sure the path exists before running the function
#' @import Seurat
#' @import Matrix
#' @import ggplot2
#' @export

getSimpleGrid = function(seqscope1st,DGEdir,spatial,tiles,nrow,ncol,sidesize,outpath,collapsePath)
{
  if(missing(seqscope1st))
  {
    seqscope1st="MiSeq"
  }
  if(missing(sidesize))
  {
    sidesize=300
  }

  if (!dir.exists(DGEdir)){
    stop("DGEdir does not exist")
  }

  if (!dir.exists(outpath)){
    stop("outpath does not exist")
  }


  #install required packages
  packages = c("Matrix", "tictoc", "ggplot2", "ggsci","Seurat","mapplots","rlist","cowplot","dplyr","Rcpp")
  ## add more packages to load if needed
  ## Now load or install&load all
  package.check <- lapply(
    packages,
    FUN = function(x) {
      if (!require(x, character.only = TRUE)) {
        install.packages(x, repos="https://cran.rstudio.com", dependencies = TRUE)
        library(x, character.only = TRUE)
      }
    }
  )

  #source collapsed.cpp
  sourceCpp(collapsePath)

  setwd(DGEdir)
  biny = binx = sidesize

  #read files
  print("Read files")
  bc = read.table("barcodes.tsv",header=F)$V1
  features = read.table('features.tsv',header=F)$V2
  m = readMM('matrix.mtx')
  if(any(c(length(features),length(bc)) != dim(m)))
  {
    stop('Dimension of matrix.mtx does not match with features or barcodes')
  }
  rownames(m) = features
  colnames(m) = bc
  m = m[,colSums(m)<=100&colSums(m)>0]  #remove outliers

  #get spatial info
  miseq_pos = read.table(spatial)
  colnames(miseq_pos) = c('HDMI','lane_miseq','tile_miseq','x_miseq','y_miseq')
  tiles = intersect(tiles,unique(miseq_pos$tile_miseq))
  bottom=miseq_pos[miseq_pos$tile_miseq %in% tiles,]

  if (seqscope1st=='MiSeq')
  {
    #bottom = miseq_pos[miseq_pos$tile>2100,]
    plotwidth = plotheight=3.5
  }
  else
  {
    #bottom = miseq_pos
    plotheight=3.5
    plotwidth=plotheight*3
  }

  print('merge')

  df = data.frame("HDMI" =colnames(m),"HDMIind" = 1:(dim(m)[2]))
  tile_df = merge(bottom,df,by = "HDMI")
 
  #aggregate all tils by expanding coord
  tile_df$aggrInd =  as.numeric(factor(tile_df$tile_miseq))-1
  #tile_df$aggrInd = tile_df$tile_miseq - min(tile_df$tile_miseq)
  m_tile = m[,tile_df$HDMIind]
  tile_df$UMI = colSums(m_tile)
  tile_df$tileHDMIind= match(tile_df$HDMI,colnames(m_tile))
  setwd(outpath)
  print('Start collapsing')
  #change the following to either cpp or apply
  obj=sapply(as.list(tiles),collapseTiles, tile_df=tile_df,binx=binx,biny=biny,m_tile=m_tile)
  obj=obj[[1]]
  
  #remodifying coordinates and super tile
  #for super tile
  tile_df = obj@meta.data
  addson_hori = max(tile_df$Y)
  addson_verti = max(tile_df$X)
  tile_df$aggrInd =  as.numeric(factor(tile_df$tile))-1

  tile_df$aggrInd2 = 0
  tile_df[tile_df$aggrInd<ncol,"aggrInd2"] = tile_df[tile_df$aggrInd<ncol,]$aggrInd+ncol
  tile_df[tile_df$aggrInd>=ncol,"aggrInd2"] = tile_df[tile_df$aggrInd>=ncol,]$aggrInd-ncol

  tile_df$y_miseq_expand = tile_df$Y +  addson_hori*(((tile_df$aggrInd2%%ncol)))
  tile_df$x_miseq_expand =   tile_df$X +  addson_verti*(floor((tile_df$aggrInd2/ncol)))
  #tile_df$orig.ident = rownames(tile_df)
  obj@meta.data$X_expand = tile_df$x_miseq_expand
  obj@meta.data$Y_expand = tile_df$y_miseq_expand

  obj@images$image = new(
    Class = 'SlideSeq',
    assay = "Spatial",
    key = "image_",
    coordinates = obj@meta.data[,c('Y_expand','X_expand')]
  )
  saveRDS(obj,'SimpleSqureGrids.RDS')
  print('Done!')


}



getSimpleGrid(seqscope1st,DGEdir,spatial,tiles,nrow,ncol,sidesize,outpath,collapsePath)






