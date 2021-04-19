
####################################################################################################3
#' This function grids the data with user-defined binning size and collapsed HDMIs within each grid
#' @param seqscope1st Data flatform, please use "MiSeq" for now
#' @param DGEdir folder that stores barcodes.tsv, features.tsv and matrix.mtx
#' @param spatial txt file stores spatial informaiton with four columns: 'HDMI','tile_miseq','x_miseq','y_miseq'
#' @param tile tile number. ONLY one tile number is allowed in current version. We are currently working on the scalability.
#' @param nrow an integer of how many rows to organize the tile
#' @param ncol an integer of how many cols to organize the tile
#' @param sidesize side size of the square grid 300 represents 10um
#' @param outpath path to store the output, you need to make sure the path exists before running the function
#' @param window size of sliding window
#' @param subXmin min value of X axis of the sub-field
#' @param subXmax max value of X axis of the sub-field
#' @param subYmin min value of Y axis of the sub-field
#' @param subYmax max value of Y axis of the sub-field
#' @export
getSlidingGrid = function(seqscope1st,DGEdir,spatial,tile,nrow,ncol,sidesize,outpath,window,subXmin,subXmax,subYmin,subYmax)
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

  biny = binx = sidesize
  setwd(DGEdir)
  #libraries
  # packages = c("Rcpp","Matrix", "tictoc", "ggplot2", "ggsci","Seurat","rlist","cowplot","dplyr","mapplots","tictoc")
  # ## add more packages to load if needed
  # ## Now load or install&load all
  # package.check <- lapply(
  #   packages,
  #   FUN = function(x) {
  #     if (!require(x, character.only = TRUE)) {
  #       install.packages(x, repos="https://cran.rstudio.com", dependencies = TRUE)
  #       library(x, character.only = TRUE)
  #     }
  #   }
  # )

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
  dim(m)
  #get spatial info
  miseq_pos = read.table(spatial)
  #colnames(miseq_pos) = c('HDMI','tile_miseq','x_miseq','y_miseq')
  colnames(miseq_pos) = c('HDMI','lane_miseq','tile_miseq','x_miseq','y_miseq')

  if (seqscope1st=='MiSeq')
  {
    bottom = miseq_pos[miseq_pos$tile>2100,]
    plotwidth = plotheight=3.5
  }
  else
  {
    bottom = miseq_pos
    plotheight=3.5
    plotwidth=plotheight*3
  }


  df = data.frame("HDMI" =colnames(m),"HDMIind" = 1:(dim(m)[2]))
  df_merge = merge(bottom,df,by = "HDMI")
  head(df_merge)

  tile_df = df_merge[df_merge$tile_miseq %in% tile,]
  setwd(outpath)
  m_tile = m[,tile_df$HDMIind]
  tile_df$UMI = colSums(m_tile)
  tile_df$tileHDMIind= match(tile_df$HDMI,colnames(m_tile))
  tile_df$aggrInd =  as.numeric(factor(tile_df$tile_miseq))-1


   #focus on given tile, sub-field
    tile_df_exact = tile_df[tile_df$tile_miseq==tile,]
    tile_df_sub = tile_df_exact[tile_df_exact$x_miseq> subXmin&tile_df_exact$x_miseq< subXmax & tile_df_exact$y_miseq< subYmax & tile_df_exact$y_miseq> subYmin,]

    miny = min(tile_df_sub$y_miseq)
    maxy = max(tile_df_sub$y_miseq)
    minx = min(tile_df_sub$x_miseq)
    maxx = max(tile_df_sub$x_miseq)
    slidestarts = seq(0,(binx/window-1),1)

    print("Start Sliding Square Grids!")
    #write these in Rcpp.
    for(j in slidestarts)
    {

      #print(j)

      for(t in slidestarts)
      {
        miny = min(tile_df_sub$y_miseq)
        maxy = max(tile_df_sub$y_miseq)
        minx = min(tile_df_sub$x_miseq)
        maxx = max(tile_df_sub$x_miseq)

        tile_df_sub_down = tile_df_sub[tile_df_sub$x_miseq>=(minx+window*j)&tile_df_sub$x_miseq<=(maxx-window*j),]
        miny = min(tile_df_sub_down$y_miseq)
        maxy = max(tile_df_sub_down$y_miseq)
        minx = min(tile_df_sub_down$x_miseq)
        maxx = max(tile_df_sub_down$x_miseq)

        #print(t)

        tile_df_sub_wind = tile_df_sub_down[tile_df_sub_down$y_miseq>=(miny+window*t) & tile_df_sub_down$y_miseq<=(maxy-window*t),]
        xlim2 = c(min(tile_df_sub_wind$x_miseq),max(tile_df_sub_wind$x_miseq))
        ylim2 = c(min(tile_df_sub_wind$y_miseq),max(tile_df_sub_wind$y_miseq))
        grd = make.grid(tile_df_sub_wind$x_miseq,tile_df_sub_wind$y_miseq,tile_df_sub_wind$UMI, binx,biny, xlim2, ylim2)
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
        grd_re = make.grid(tile_df_sub_wind$x_miseq,tile_df_sub_wind$y_miseq,tile_df_sub_wind$tileHDMIind, binx, biny, xlim2, ylim2,function(x) {write(length(x), file=fn1,append = T)})
        grd_re = make.grid(tile_df_sub_wind$x_miseq,tile_df_sub_wind$y_miseq,tile_df_sub_wind$tileHDMIind, binx, biny, xlim2, ylim2,function(x) {cat(x,file=fn2,append=TRUE,sep='\n')})


        #draw the grids centers:
       # write.csv(as.numeric(colnames(grd)),paste0('grd_col',j,t,'.csv'))
      #  write.csv(as.numeric(rownames(grd)),paste0('grd_row',j,t,'.csv'))

        collapseLen = read.csv(fn1,header=F)
        collapseInd = read.table(fn2,header=F)


        collapseLen = cbind(collapseLen,cumsum(collapseLen$V1))
        colnames(collapseLen) =c("len","end")
        interv = c(0,collapseLen$end)

        #create dataframe of the assignment of collapsed grids for each HDMI
        df=data.frame('HDMIind' = collapseInd$V1,"HDMI" = colnames(m_tile)[collapseInd$V1])
        assign=c()
        out=sapply(1:(dim(collapseLen)[1]),function(x) {nrep=collapseLen$len[x];return(c(assign,rep(paste0('Collapse2_',x),nrep)))})
        assign = unlist(out)
        df$assign = assign
        #fn3=paste0('HDMI_collapsing_assignment','_',j,'_',t,'.csv')

      #  write.csv(df,fn3)

        #sourceCpp(collapsePath)
        tic();collapseM = collapse(m_tile,collapseInd$V1,interv);toc()
        rownames(collapseM) = rownames(m_tile)
        colnames(collapseM) = paste0("Collase_",j,"_",t,1:(length(interv)-1))
        sparse.gbm <- Matrix(collapseM , sparse = T )
        # fn4=paste0('tile2112_collapsedMatrix_',j,'_',t,'.mtx')
        # fn5=paste0('tile2112_collapsedGenes_',j,'_',t,'.csv')
        # fn6=paste0('tile2112_collapsedBarcodes_',j,'_',t,'.csv')
        #
        # writeMM(obj = sparse.gbm, file=fn4)
        # write.csv(rownames(collapseM),fn5)
        # write.csv(colnames(collapseM),fn6)

        grd=t(grd)
        pos = which(!is.na(grd), TRUE)
        pos_coor = t(sapply(1:(dim(pos)[1]),function(x) {c(as.numeric(rownames(grd)[as.numeric(pos[x,1])]),as.numeric(colnames(grd)[as.numeric(pos[x,2])]) )})) #here the grd is transposed
        colnames(pos_coor) =c("X","Y")
        coord.df = data.frame("Y"=pos_coor[,2], "X"=pos_coor[,1],"tile"=tile, stringsAsFactors=FALSE)
        fn7=paste0('tile2112_coor_df_',j,'_',t,'.csv')
        #write.csv(coord.df,fn7)

        dge1 = collapseM
        spatial1 = coord.df
        obj1 = CreateSeuratObject(counts=dge1,assay='Spatial')
        #obj1$status = "Original"
        obj1@meta.data$tile = spatial1$tile
        obj1@meta.data$X = spatial1$X
        obj1@meta.data$Y = spatial1$Y
        obj1@meta.data$interationi=j
        obj1@meta.data$interationj=t

        if(j==0 & t==0)
        {
          obj=obj1
        }
        else
        {
          obj = merge(obj,obj1)

        }
      }
    }

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
    obj@meta.data$X_expand = tile_df$x_miseq_expand
    obj@meta.data$Y_expand = tile_df$y_miseq_expand

    obj@images$image = new(
      Class = 'SlideSeq',
      assay = "Spatial",
      key = "image_",
      coordinates = obj@meta.data[,c('Y_expand','X_expand')]
    )
    saveRDS(obj,'SlidingSquareGrids.RDS')
    print('Done!')
}



#
#
#   #setwd('~/scrna/leejun/ngst/fastqs/HiSeq/19129-9and11/ColonCombAll-Update-9and11/analysis/Seurat/grid_10um/fourtimesGrid/')
#
#
#
#   protein_coding = read.csv('~/TenXGenes_GMremoved.csv',row.names=1)
#   protein_coding=protein_coding$V3
#
#   dge1 = c()
#   spatial1 = c()
#
#   for(j in slidestarts)
#   {
#    print(j)
#     for (t in slidestarts)
#     {
#       fn4=paste0('tile2112_collapsedMatrix_',j,'_',t,'.mtx')
#       fn5=paste0('tile2112_collapsedBarcodes_',j,'_',t,'.csv')
#       fn6=paste0('tile2112_collapsedGenes_',j,'_',t,'.csv')
#       fn7=paste0('tile2112_coor_df_',j,'_',t,'.csv')
#
#       m1 = readMM(fn4)
#       bc1 = read.csv(fn5)
#       gene1=read.csv(fn6)
#
#       colnames(m1) =paste0(i,'_',bc1$x)
#       rownames(m1) =gene1$x
#       dge1 = m1
#       coor1 = read.csv(fn7,row.names=1)
#       spatial1 = coor1
#       ind = match(protein_coding,rownames(dge1))
#       ind = ind[complete.cases(ind)]
#       dge1 = dge1[ind,]
#       dge1 = dge1[!duplicated(rownames(dge1)),]
#
#       #protein_coding
#       obj1 = CreateSeuratObject(counts=dge1,assay='Spatial')
#       #obj1$status = "Original"
#       obj1@meta.data$tile = spatial1$tile
#       obj1@meta.data$X = spatial1$X
#       obj1@meta.data$Y = spatial1$Y
#       obj1@meta.data$interationi=j
#       obj1@meta.data$interationj=t
#
#       # obj1@images$image = new(
#       #   Class = 'SlideSeq',
#       #   assay = "Spatial",
#       #   key = "image_",
#       #   coordinates = spatial1
#       # )
#       if(j==0 & t==0)
#       {
#         obj=obj1
#       }
#       else
#       {
#         obj = merge(obj,obj1)
#
#       }
#
#
#     }
#
#   }
#
#
#   #for super tile
#   tile_df = obj@meta.data
#   addson_hori = max(tile_df$Y)
#   addson_verti = max(tile_df$X)
#   tile_df$aggrInd =  as.numeric(factor(tile_df$tile))-1
#   tile_df$aggrInd2 = 0
#   tile_df[tile_df$aggrInd<ncol,"aggrInd2"] = tile_df[tile_df$aggrInd<ncol,]$aggrInd+ncol
#   tile_df[tile_df$aggrInd>=ncol,"aggrInd2"] = tile_df[tile_df$aggrInd>=ncol,]$aggrInd-ncol
#   tile_df$y_miseq_expand = tile_df$Y +  addson_hori*(((tile_df$aggrInd2%%ncol)))
#   tile_df$x_miseq_expand =   tile_df$X +  addson_verti*(floor((tile_df$aggrInd2/ncol)))
#   obj@meta.data$X_expand = tile_df$x_miseq_expand
#   obj@meta.data$Y_expand = tile_df$y_miseq_expand
#
#   obj@images$image = new(
#     Class = 'SlideSeq',
#     assay = "Spatial",
#     key = "image_",
#     coordinates = obj@meta.data[,c('Y_expand','X_expand')]
#   )
#   saveRDS(obj,'SlidingSquareGrids.RDS')


 # }

#
#   png('spatialBeforeCutoff.png',height=3,width=6,res=300,units='in')
#   SpatialFeaturePlot(obj, features = "nCount_Spatial",images = 'image') + theme(legend.position = "right")
#   dev.off()
#
#   png('test1.png',height=6,width=6,res=300,units='in')
#   x = obj@meta.data
#   qplot(y=x$Y,x=x$X,color=as.factor(x$interationi))+geom_point(size=0.1,alpha=1);dev.off()
#
#
#
#   png('test1.png',height=6,width=6,res=300,units='in')
#   qplot(y=x$Y,x=x$X,color=x$interationi);dev.off()
#   #SpatialFeaturePlot(, features = "nCount_Spatial",images = 'image') + theme(legend.position = "right")
#   dev.off()
#
#   saveRDS(obj,'Liver_10Xlist_Gmremoved_noisoformdup_tile2118_middle_Sliding100Times.RDS')
#
#   #saveRDS(obj,'OriginalAndFourTimes.RDS')
#
#   #
#   # #remodifying coordinates and super tile
#   # tile_df = obj@meta.data
#   # tile_df$aggrInd =  as.numeric(factor(tile_df$tile))-1
#   # nrow = 2
#   # ncol=6
#   # tile_df$aggrInd2 = 0
#   # tile_df[tile_df$aggrInd<ncol,"aggrInd2"] = tile_df[tile_df$aggrInd<ncol,]$aggrInd+ncol
#   # tile_df[tile_df$aggrInd>=ncol,"aggrInd2"] = tile_df[tile_df$aggrInd>=ncol,]$aggrInd-ncol
#   #
#   # tile_df$y_miseq_expand = tile_df$Y + 30000*(((tile_df$aggrInd2%%ncol)))
#   #
#   # tile_df$x_miseq_expand =   tile_df$X+ 30000*(floor((tile_df$aggrInd2/ncol)))
#   #
#   # obj@meta.data$X_expand = tile_df$x_miseq_expand
#   # obj@meta.data$Y_expand = tile_df$y_miseq_expand
#   obj@images$image = new(
#     Class = 'SlideSeq',
#     assay = "Spatial",
#     key = "image_",
#     coordinates = obj@meta.data[,c('Y_expand','X_expand')]
#   )
#   saveRDS(obj,'OriginalAndFourTimes.RDS')
#
# }
#
#
#
# png('spatialBeforeCutoff.png',height=10*6,width=2*6,res=300,units='in')
# SpatialFeaturePlot(obj, features = "nCount_Spatial",images = 'image') + theme(legend.position = "right")
# dev.off()
#
#
#
# png('nFeature_Spatial.png',width = 5,height=4,res=300,units='in')
#
# VlnPlot(obj, features = "nFeature_Spatial", pt.size = 0,group.by='status') + NoLegend()
#
# dev.off()
#
#
#
# #ploting for miseq
# png('Original_gridding_tilei.png',res=300,width = plotwidth*ncol,height=plotheight*nrow,units='in')
# basemap(ylim, xlim, main = paste0("Super Tile"),xlab="Coord 1",ylab="Coord 2")
# p=draw.grid(grd,breaks)
# dev.off()
#
# #second grids
# miny = min(tile_df_sub_wind$y_miseq)
# maxy = max(tile_df_sub_wind$y_miseq)
# minx = min(tile_df_sub_wind$x_miseq)
# maxx = max(tile_df_sub_wind$x_miseq)
#
# xlim2 = c(minx,maxx)
# ylim2 = c(miny,maxy)
# grd = make.grid(tile_df_sub_wind$x_miseq,tile_df_sub_wind$y_miseq,tile_df_sub$UMI, binx,biny, xlim2, ylim2)
# grd = t(grd)
# breaks = breaks.grid(grd,zero=FALSE)
# #ploting for miseq
# png('new_gridding_tilei.png',res=300,width = plotwidth*ncol,height=plotheight*nrow,units='in')
# basemap(ylim2, xlim2, main = paste0("Super Tile"),xlab="Coord 1",ylab="Coord 2")
# p=drawgrids(grd,breaks)
# dev.off()
#
#
#
# #draw the grids centers:
# grd_y = as.numeric(colnames(grd))
# grd_x = as.numeric(rownames(grd))
#
# grd_y = as.numeric(colnames(grd))
# grd_x = as.numeric(rownames(grd))
#
# png('grd.png',width=7,height=8,res=300,units='in')
# plot()
# abline(v = grd_y, col="red", lwd=3, lty=2)
# dev.off()
#
#
# drawgrids=function (grd, breaks = NULL, col = NULL)
# {
#   if (missing(breaks))
#     breaks <- breaks.grid(grd)
#   ncol <- length(breaks) - 1
#   # if (missing(col))
#   #   col <- colorRampPalette(c("lightyellow", "yellow",
#   #                             "orange", "red", "brown4"))(ncol)
#   # if (length(col) != ncol)
#   #   stop("The number of breakpoints should be one more than the number of colours")
#   # grd <- ifelse(grd > max(breaks), max(breaks), grd)
#   grd <- ifelse(grd < min(breaks), min(breaks), grd)
#   x <- as.numeric(rownames(grd))
#   y <- as.numeric(colnames(grd))
#   image(x, y, grd, breaks = breaks)
#   box()
# }
#
# df2=c()
# for(i in 1:length(grd_x))
# {
#   t=data.frame('y' =grd_y,'x' = grd_x[i])
#   df2 = rbind(df2,t)
# }
#
#
# df1=c()
# for(i in 1:length(grd_x))
# {
#   t=data.frame('y' =grd_y,'x' = grd_x[i])
#   df1 = rbind(df1,t)
# }
#
#
#
# df1$group = "Red"
# df2$group= "Green"
# df = rbind(df1,df2)
# png('grd.png',width=7,height=8,res=300,units='in')
# ggplot(df, aes(x=x, y=y,color = group)) +
#   geom_tile(fill='transparent') +
#   scale_y_reverse() +
#   theme_classic() +
#   theme(axis.text  = element_blank(),
#         panel.grid = element_blank(),
#         axis.line  = element_blank(),
#         axis.ticks = element_blank(),
#         axis.title = element_blank())
# dev.off()
#
# library(ggplot2)
# library(scales)
# ggplot(df2, aes(x=x, y=y)) +
#   geom_tile(fill='transparent', colour = 'red') +
#   scale_y_reverse() +
#   theme_classic() +
#   theme(axis.text  = element_blank(),
#         panel.grid = element_blank(),
#         axis.line  = element_blank(),
#         axis.ticks = element_blank(),
#         axis.title = element_blank())
# par(new=TRUE)
# ggplot(df1, aes(x=x, y=y)) +
#   geom_tile(fill='transparent', colour = 'green') +
#   scale_y_reverse() +
#   theme_classic() +
#   theme(axis.text  = element_blank(),
#         panel.grid = element_blank(),
#         axis.line  = element_blank(),
#         axis.ticks = element_blank(),
#         axis.title = element_blank())
# dev.off()
#
# # qplot(y=grd_y,x =grd_x);
#
#
# tile_i=subset(obj,subset=(tile==2103))
# png('spatialAfterCutoff500.png',height=10*6,width=2*6,res=300,units='in')
# SpatialPlot(tile_i, features = "nCount_Spatial") + theme(legend.position = "right")
# dev.off()
#
#
# #cluster concordance:
# tile_2110=ngst@meta.data[ngst@meta.data$tile==2113,]
# tile_2110_four = tile_2110[tile_2110$status=='Sliding_4times',]
# tile_2110_Original = tile_2110[tile_2110$status=='Original',]
# item_2 = as.matrix(tile_2110_four[,c('X_expand','Y_expand')])
# i_vec = c()
# or_clus=c()
# new_clus=c()
# new_vec=c()
# for (i in 1:dim(tile_2110_Original)[1]){
#   print(i)
#   or = tile_2110_Original[i,]
#   four_i = tile_2110_four[tile_2110_four$X_expand==or$X_expand+150 & tile_2110_four$Y_expand==or$Y_expand-150,]
#   if(dim(four_i)[1]==0)
#   {next}
#   or_clus = c(or_clus,or$seurat_cluster)
#   new_vec=rbind(new_vec,four_i)
#   i_vec = c(i_vec,i)
#   new_clus = c(new_clus,four_i$seurat_cluster)
#   #dist_i = dist_euclidean(as.matrix(or[,c('X_expand','Y_expand')]),item_2)
#   #ind=sort(dist_i[dist_i!=0],decreasing = F,index.return=T)
#   #ind[[2]] [1:4]
# }
# out_vec = tile_2110_Original[i_vec,]
# out_vec$new_seurat_clusters = new_clus-1
#
#
# png('tile2113_mixture.png',width=7,height=6,res=300,units='in')
# ggplot(out_vec,aes(Y,X,color = new_seurat_clusters==seurat_clusters))+geom_point(size=1,alpha=0.9)+theme_bw()+theme(legend.position = 'bottom');
# dev.off()
#
# png('tile2110_1.png',width=7,height=6,res=300,units='in')
# ggplot(out_vec,aes(X,Y,color = seurat_clusters))+geom_point(size=1,alpha=0.8)+theme_bw()+theme(legend.position = 'bottom');
# dev.off()
# png('tile2110_2.png',width=7,height=6,res=300,units='in')
# ggplot(out_vec,aes(X,Y,color = as.factor(new_seurat_clusters)))+geom_point(size=1,alpha=0.8)+theme_bw()+theme(legend.position = 'bottom');
# dev.off()
#
# out_old_new = rbind(tile_2110_Original[i_vec,],new_vec)
#
#
# png('tile2110_merge.png',width=7,height=6,res=300,units='in')
# ggplot(out_old_new,aes(X,Y,color = as.factor(seurat_clusters)))+geom_point(size=1,alpha=0.8)+theme_bw()+theme(legend.position = 'bottom');
# dev.off()
#
#
# png('tile2113_merge.png',width=7,height=7,res=300,units='in')
# ggplot(out_old_new,aes(Y,X,color = as.factor(seurat_clusters)))+geom_point(size=1,alpha=0.9)+theme_bw()+theme(legend.position = 'bottom');
# dev.off()
#
#
# png('tile2113_1.png',width=7,height=7,res=300,units='in')
# ggplot(out_vec,aes(Y,X,color = seurat_clusters))+geom_point(size=1,alpha=0.9)+theme_bw()+theme(legend.position = 'bottom');
# dev.off()
#
# png('tile2110_dim_plot.png',width=7,height=6,res=300,units='in')
# SpatialDimPlot(out_old_new)
# dev.off()
#
#
# png('dim_plot.png',width=7,height=6,res=300,units='in')
# SpatialDimPlot(ngst,group)
# dev.off()
#
#
#
#
#
# #raster
# meta = tile2117@meta.data
# points <- SpatialPoints(meta[,c('X','Y')], meta[,c('nCount_Spatial')])
#
#
#
# s100 <- matrix(c(267573.9, 2633781, 213.29545, 262224.4, 2633781, 69.78261, 263742.7, 2633781, 51.21951, 259328.4, 2633781, 301.98413, 264109.8, 2633781, 141.72414, 255094.8, 2633781, 88.90244),  ncol=3,  byrow=TRUE)
# colnames(s100) <- c('X', 'Y', 'Z')
#
# library(raster)
# # set up an 'empty' raster, here via an extent object derived from your data
# e <- extent(s100[,1:2])
# e <- e + 1000 # add this as all y's are the same
#
# r <- raster(e, ncol=10, nrow=2)
# # or r <- raster(xmn=, xmx=,  ...
#
# # you need to provide a function 'fun' for when there are multiple points per cell
# x <- rasterize(s100[, 1:2], r, s100[,3], fun=mean)
# plot(x)
#
#
#
# s100 = meta[,c('X','Y','nCount_Spatial')]
# colnames(s100) <- c('x', 'y', 'z')
# e <- extent(s100[,1:2])
# e <- e + 1000
# r <- raster(e, ncol=92, nrow=81)
# x <- rasterize(s100[, 1:2], r, s100[,3], fun=mean)
# png('test2.png',width=7,height=6,units='in',res=300)
# plot(x)
# dev.off()
# plot(x)
#
#
#
# # Following is the map of ELSA, lower values represent higher local spatial autocorrelation
# png('tile2117_elsa_300dist.png',units='in',res=300,width=7,height=7)
# e <- elsa(x,d=320,categorical=FALSE)
# cl <- colorRampPalette(c('darkblue','yellow','red','black'))(100) # specifying a color scheme
# plot(e,col=cl, main='ELSA for the local distance of 10um')
# #plot(x, main="local Moran's I (Z.Ii)")
#
# dev.off()
#

# #example: Colon data tile 2112 top left
# #window = 2/10*300
# window = 2/10*300
#
# DGEdir='/net/fantasia/home/jyxi/scrna/leejun/ngst/fastqs/HiSeq/19129-9and11/ColonCombAll-Update-9and11/analysis/align/ColonCombAll-Update-9and11-STARsolo_trimA10_bottomSolo.out/GeneFull/raw/'
# #nrow = 2
# #ncol=6
# spatial = '~/scrna/leejun/ngst/fastqs/MiSeq-DraI-100pM-mbcore-RD4-revHDMIs-pos-uniq.txt'
# tile = 2110
#
# #binx=biny=150 #5um
# binx=biny=sidesiz=300 #5
# outpath="/net/fantasia/home/jyxi/scrna/leejun/ngst/fastqs/HiSeq/19129-9and11/ColonCombAll-Update-9and11/analysis/Seurat/slide2um/tile2110"
# collapsePath = '/net/fantasia/home/jyxi/scrna/leejun/ngst/fastqs/HiSeq/ColonCombAll/analysis/Seurat/collapse.cpp'
#
#
#
# ##########Liver data
# window = 1/10*300
# DGEdir='/net/fantasia/home/jyxi/scrna/leejun/ngst/fastqs/HiSeq/LiverCombAll/analysis/align/LiverCombAll-Starsolo-trimA10-bottomSolo.out/GeneFull/raw'
# nrow = 2
# ncol=6
# spatial = "~/scrna/leejun/ngst/fastqs/MiSeq-DraI-100pM-mbcore-RD2-revHDMIs-pos-uniq.txt"
# tile = 2118
# #tile_df_file = "/net/fantasia/home/jyxi/scrna/leejun/ngst/fastqs/HiSeq/LiverCombAll/analysis/Seurat/tile_df.csv"
#
# #binx=biny=150 #5um
# binx=biny=sidesize=300 #5
# outpath='/net/fantasia/home/jyxi/scrna/leejun/ngst/fastqs/HiSeq/LiverCombAll/analysis/Seurat/slide_1um/tile2118'
# collapsePath = '/net/fantasia/home/jyxi/scrna/leejun/ngst/fastqs/HiSeq/ColonCombAll/analysis/Seurat/collapse.cpp'


#add a parameter of isoform?
#getSlidingGrid(seqscope1st="MiSeq",DGEdir,spatial,tile,nrow,ncol,sidesize,outpath,window =150,subXmin,subXmax,subYmin,subYmax)

