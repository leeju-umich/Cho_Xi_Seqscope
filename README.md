# Cho_Xi_Seqscope
Seq-Scope is a spatial barcoding technology with a resolution almost comparable to an optical microscope. Seq-Scope is based on a solid-phase amplification of randomly barcoded single-molecule oligonucleotides using an Illumina sequencing-by-synthesis platform. The resulting clusters annotated with spatial coordinates are processed to expose RNA-capture moiety. These RNA-capturing barcoded clusters define the pixels of Seq-Scope that are approximately 0.5-1 μm apart from each other. For more information, please refer to the link''.

This github page includes the pipepline and codes for processing the data, including tissue boundary detection, alignment, and downstream analysis(gridding, collapsing, clustering analysis, etc.).

## Requirements
linux, STARSolo>=2.7.5c, seqtk, R, Python, perl...

## Data Access
The two datasets used in the paper (link) can be found as following:
- Liver Data
  - 1st Seq data: https://drive.google.com/file/d/1vv6Uy7Ovvw-Sx3WCnRw0tgwuyvrKdVfo/view?usp=sharing
  - 2nd Seq data: (too big)
- Colon Data
  - 1st Seq data: 
  - 2nd Seq data: (too big)
 
 The structures of the the data can be found in the original paper. 
## Tissue Boundary Estimation
To estimate the tissue boundary, the 2nd Seq data are joined into 1st Seq data according to their HDMI sequence. HDMI discovery plot is generated to visualize the density of HiSeq HDMI in a given XY space of each tile. Then density plots can be manually assigned to the corresponding H&E images. We use liver data as an example. 

- Step 1: Extract HDMI and spatial information from 1st Seq; Extract HDMI information from 2nd Seq.
  - Input files:xxx
  - Linux codes:xxx
    ``` zcat /net/fantasia/home/jyxi/scrna/leejun/ngst/fastqs/DraI-100pM-mbcore-RD2.fastq.gz | sed -n '1~4s/:/ /gp' | cut -d ' ' -f 5-7 >  /net/fantasia/home/jyxi/scrna/leejun/ngst/fastqs/pos-MiSeq-DraI-100pM-mbcore-RD2.txt
    zcat /net/fantasia/home/jyxi/scrna/leejun/ngst/fastqs/DraI-100pM-mbcore-RD2.fastq.gz  | perl -lane 'print $_ if ( $. % 4 == 2 )'  | cut -c 1-20  > /net/fantasia/home/jyxi/scrna/leejun/ngst/fastqs/HDMIs-MiSeq-DraI-100pM-mbcore-RD2.txt
     cat /net/fantasia/home/jyxi/scrna/leejun/ngst/fastqs/HDMIs-MiSeq-DraI-100pM-mbcore-RD2.txt | rev | tr ACGTN TGCAN > /net/fantasia/home/jyxi/scrna/leejun/ngst/fastqs/HDMIs-MiSeq-DraI-100pM-mbcore-RD2-rev.txt
     paste ~/scrna/leejun/ngst/fastqs//HDMIs-MiSeq-DraI-100pM-mbcore-RD2-rev.txt ~/scrna/leejun/ngst/fastqs/pos-MiSeq-DraI-100pM-mbcore-RD2.txt | column -s $'\t' -t > ~/scrna/leejun/ngst/fastqs/MiSeq-DraI-100pM-mbcore-RD2-revHDMIs-pos.txt
     zcat /net/fantasia/home/jyxi/scrna/leejun/ngst/fastqs/HiSeq/ColonCombAll/combR1.fastq.gz | perl -lane 'print $_ if ( $. % 4 == 2 )'  | cut -c 1-20 > /net/fantasia/home/jyxi/scrna/leejun/ngst/fastqs/HiSeq/ColonCombAll/HDMI-20mers-from-HiSeq-colon.txt

    ```
   - Output files:xxxxxx
- Step 2 : To join 2nd Seq data  into 1st Seq data and estimate the tissue boundary, please refer to the python codes below.
  - Input files:
  - Python codes:
     ```import pandas as pd
     miseq_pos = pd.read_csv('LiverCombo/MiSeq-DraI-100pM-mbcore-RD2-revHDMIs-pos-uniq.txt',delim_whitespace=True, header=None)
     miseq_pos.columns=['HDMI','tile','x','y']
     hiseq = pd.read_csv("All-20mers-from-HiSeq-liver.txt",delim_whitespace = True,header=None)
     hiseq.columns = ['HDMI']
     merge_df = pd.merge(hiseq,miseq_pos, on='HDMI',how="inner")
     
     
    import mpl_scatter_density # adds projection='scatter_density'
    from matplotlib.colors import LinearSegmentedColormap
    import pylab as plt

    white_viridis = LinearSegmentedColormap.from_list('white_viridis', [
        (0, '#ffffff'),
        (1e-20, '#440053'),
        (0.2, '#404388'),
        (0.4, '#2a788e'),
        (0.6, '#21a784'),
        (0.8, '#78d151'),
        (1, '#fde624'),
    ], N=256)

    def using_mpl_scatter_density(fig, x, y):
        ax = fig.add_subplot(1, 1, 1, projection='scatter_density')
        density = ax.scatter_density(x, y, cmap=white_viridis)
        #density = ax.scatter_density(x, y, cmap=white_viridis,vmax =200) # can change vmax here to eliminate outliers.
        fig.colorbar(density, label='Number of points per pixel')

    for i in range(2101,2120):
        x = merge_df[merge_df.tile.eq(i)]
        fig = plt.figure()
        using_mpl_scatter_density(fig, x['y'], x['x'])
        plt.gca().set_aspect('equal', adjustable='box')
        plt.savefig('tile'+str(i)+'.png',dpi=1000)
        plt.show()
     ```
   - Output files
## Alignment
STARsolo software is used to align the reads.
  - Step 1: Trimming Read 1 (this helps with the computational timing) and copy randomer to Read 1.
    - Input files
    - Linux Codes
    
    ```
    #Trimming R1
     ~/seqtks/seqtk trimfq -q 0 -l 29 /net/fantasia/home/jyxi/scrna/leejun/ngst/fastqs/HiSeq/19129FL0601/19129FL-06-01-01_S93_L008_R1_001.fastq.gz >            /net/fantasia/home/jyxi/scrna/leejun/ngst/fastqs/HiSeq/19129FL0601/19129FL-06-01-01_S93_L008_R1_001_trimmed.fastq
     

    pigz -p 8 
    /net/fantasia/home/jyxi/scrna/leejun/ngst/fastqs/HiSeq/19129FL0601/19129FL-06-01-01_S93_L008_R1_001_trimmed.fastq


    #Copy and paste randomer from R2 to R1
    paste <(zcat /net/fantasia/home/jyxi/scrna/leejun/ngst/fastqs/HiSeq/19129FL0601/19129FL-06-01-01_S93_L008_R1_001_trimmed.fastq.gz) <(zcat /net/fantasia/home/jyxi/scrna/leejun/ngst/fastqs/HiSeq/19129FL0601/19129FL-06-01-01_S93_L008_R2_001.fastq.gz) | perl -lane 'if ( $. % 4 == 1 ) { print "$F[0] $F[1]"; } elsif ( $. % 4 == 3 ) { print "+"; } else { print substr($F[0],0,20).substr($F[1],0,9).substr($F[0],50); }' > /net/fantasia/home/jyxi/scrna/leejun/ngst/fastqs/HiSeq/19129FL0601/19129FL-06-01-01_S93_L008_R1_001_final.fastq

    pigz -p 8 /net/fantasia/home/jyxi/scrna/leejun/ngst/fastqs/HiSeq/19129FL0601/19129FL-06-01-01_S93_L008_R1_001_final.fastq
    ```
    - Output files
  - Step 2: 
    - Input files
  
    - STARsolo alignment:
   ```
  ./STAR    --genomeDir  /net/fantasia/home/jyxi/scrna/dropseq/mouse/mm10/geneIndex \
          --readFilesIn  /net/fantasia/home/jyxi/scrna/leejun/ngst/fastqs/HiSeq/19129FL0601/19129FL-06-01-01_S93_L008_R2_001.fastq.gz /net/fantasia/home/jyxi/scrna/leejun/ngst/fastqs/HiSeq/19129FL0601/19129FL-06-01-01_S93_L008_R1_001_final.fastq.gz  \
          --outSAMtype BAM SortedByCoordinate  \
          --readFilesCommand zcat \
          --runDirPerm All_RWX \
          --outFileNamePrefix /net/fantasia/home/jyxi/scrna/leejun/ngst/fastqs/HiSeq/19129FL0601/analysis/align/N709-liver-Starsolo-trimA10-bottom \
          --soloType CB_UMI_Simple \
          --soloCBstart 1 --soloCBlen 20 \
          --soloUMIstart 21 --soloUMIlen 9 \
          --soloCBwhitelist 
/net/fantasia/home/jyxi/scrna/leejun/ngst/fastqs/HiSeq/19129-04/PsomagenAdmera/N709N710_bottom_tiles_concord_whitelists.txt \
          --runThreadN 6 \
          --soloBarcodeReadLength 0 \
          --outSAMattributes NH HI nM AS CR UR CB UB GX GN sS sQ sM \
          --outFilterScoreMinOverLread 0 \
          --outFilterMatchNminOverLread 0 \
          --clip3pAdapterSeq AAAAAAAAAA \
          --clip3pAdapterMMp 0.1 \
          --soloFeatures Gene GeneFull SJ Velocyto \
          --soloCellFilter None 

   ```
    - Output files

## Data gridding and collapsing (need to remodifying the funcition to overwrite the ind.txt and some default inputs,change Original.RDS to SeqScope.RDS?)
Square gridding with user-specified binning size is conducted on the raw digital expression matrix genearted by STARsolo alignment. 
  - Input files: raw digital expression matrix generated by STARsolo. Note that with the above alignment pipeline, raw digital expression matrix is generated under Gene and GeneFull folder(corresponding to 'Gene' and 'GeneFull' option). Either is fine to be used for collapsing and downstream analysis. 
  - Codes
  Please refer to the main function squareGrid.R and its documentaion in the repo. The example usage of the function is:
  ```
  workingdir = '~/scrna/leejun/ngst/fastqs/HiSeq/ColonCombAll/analysis/align/ColonCombAll_Starsolo_trimA10_bottomSolo.out/Gene/raw'
spatial = '~/scrna/leejun/ngst/fastqs/MiSeq-DraI-100pM-mbcore-RD4-revHDMIs-pos-uniq.txt'
#colon: 
tiles = c(2103:2106,2110:2114,2118:2119)
nrow = 2
ncol=6
#binx=biny=150 #5um
binx=biny=300 #10um

outpath = '/net/fantasia/home/jyxi/scrna/leejun/ngst/fastqs/HiSeq/ColonCombAll/analysis/Seurat/grid_10um_Gene/'
collapsePath = '/net/fantasia/home/jyxi/scrna/leejun/ngst/fastqs/HiSeq/ColonCombAll/analysis/Seurat/collapse.cpp'
  ```
  - The function squareGrid outputs intermediant and final files. The final file that's to be used for further analysis is Original.R. The spatial coordinates information is stored in the meta data; Also, you can call 'GetTissueCoordinates' Seurat function for spatial information.

## Clustering and find gene markers
We will use Seurat R package for downstream analysis including clustering and finding cluster specific gene markers. 
- Input files: Original.RDS from gridding and collapsing step
- Codes: The examples codes is as follows
```
library('ggplot2')
library('Seurat')
library('cowplot')
obj = readRDS('Original.RDS')
png('nCount_Spatial.png',width = 5,height=4,res=300,units='in')
VlnPlot(obj, features = "nCount_Spatial", pt.size = 0) + NoLegend()
dev.off()

png('nFeature_Spatial.png',width = 5,height=4,res=300,units='in')
VlnPlot(obj, features = "nFeature_Spatial", pt.size = 0) + NoLegend()
dev.off()

#remove non-tissue areas by setting a cutoff of nFeature from the feature plot
feature_cutoff = 100
ngst = subset(obj, subset = nFeature_Spatial > feature_cutoff) #clon comb
#spatial plot to visualize the tissue area after filtering out non tissue areas
SpatialFeaturePlot(ngst, features = "nCount_Spatial") + theme(legend.position = "right")

#find clusters and dimmension reduction
ngst = SCTransform(ngst, assay = "Spatial", verbose = FALSE)
ngst= RunPCA(ngst, assay = "SCT", verbose = FALSE)
ngst = FindNeighbors(ngst, reduction = "pca", dims = 1:30)
ngst = FindClusters(ngst, verbose = FALSE)
ngst = RunUMAP(ngst, reduction = "pca", dims = 1:30)
ngst = RunTSNE(ngst, reduction = "pca", dims = 1:30,check_duplicates = FALSE)
saveRDS(ngst,'Original_cutoff.RDS')
p1 = DimPlot(ngst, reduction = "tsne", label = F)
p2 = DimPlot(ngst, reduction = "umap", label = F)
plot_grid(p1,p2,ncol=2,nrow=1)


#find cluster specific gene markers:
ngst.markers = FindAllMarkers(ngst, only.pos = TRUE, min.pct = 0.1, logfc.threshold = 0.1)
nGene = 10
topn = ngst.markers %>% group_by(cluster) %>% top_n(n = nGene, wt = avg_logFC)
write.csv(topn,'FindMarkers.csv')
```

## Spatial plotting

