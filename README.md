# Seq-Scope Data Analysis Pipeline

Seq-Scope is a spatial barcoding technology with a resolution almost comparable to an optical microscope. Seq-Scope is based on a solid-phase amplification of randomly barcoded single-molecule oligonucleotides using an Illumina sequencing-by-synthesis platform. The resulting clusters annotated with spatial coordinates are processed to expose RNA-capture moiety. These RNA-capturing barcoded clusters define the pixels of Seq-Scope that are approximately 0.5-1 μm apart from each other. For more information, please refer to the link (https://doi.org/10.1101/2021.01.25.427807).

This is a brief tutorial that includes scripts used for the SeqScope paper. The bash scripts are used for preprocessing the data (tissue boundary detection, STARsolo alignment), python and R scipts to bin the data into square grids and conduct part of subcellular analysis. All scripts can be found under the script folder in this repository. 

**Note**: All scripts are writen based on the SeqScope sequencing structure as described in our SeqScope paper. The users will need to modify the scripts by themselves to make it compatible to their experimental design. A more flexible and user-friendly software tool is under active development. We will update this page when the new tool is ready. 

# Getting Started
## Required Sofware Tools

You need to install the following software tools before using this pipeline. Linux operating system is necessary. All scripts are running wiithin the shell.
* STARSolo>=2.7.5c
* seqtk
* R 
* Python >=3.0
* perl
* pigz

## Example Data
The raw dataset used for Seq-Scope paper will be available in GEO and SRA (GSE169706). The annotated RDS file and H&E images can be found at https://doi.org/10.7302/cjfe-wa35. Here we assume that you already have access to these example dataset. 

* 1st-seq data (typically from MiSeq)
  - abc_SeqScope_1st.fastq.gz
* 2nd-seq data (typically from NovaSeq or HiSeq X)
  - abc_SeqScope_2nd_R1.fastq.gz
  - abc_SeqScope_2nd_R2.fastq.gz
* Reference genome and STAR index
  - mm10 under GenomeDir folder in this repository
 
## Overall Workflow
This image shows the overall workflow for Seq-Scope data. We will introduce the implementations for each workflow section. 
<p>
    <img src="Workflow.png" width="1000" height="400" />
</p>

## Step1: Tissue Boundary Estimation
In this section, we process 1st-seq data to extract spatial coordinates and match the HDMIs from 1st-seq to HDMIs from 2nd-seq and to visualize the tissue boundary captured by Seq-Scope compared to H&E images. The bash script takes two file paths as arguments and outputs files in the current working directory. The tissueBoundaryPlot function visualize the tissue boundary.

### Step1a: Preprosessing
First, we need to process our data with bash script extractCoord.sh, which can be found under script folder in this repository. 

 * Input files:
  ```
  abc_SeqScope_1st.fastq.gz:  path of read file from 1st-Seq
  abc_SeqScope_2nd_R1.fastq.gz: path of read1 from 2nd-Seq
  hdmilength: An integer indicating the length of the HDMIs; For now, it can only take 20 or 30. In default, we assume if MiSeq is used for 1st-Seq, then hdmilenght=20; if HiSeq is used for 1st-Seq, then hdmilength=30.
  ```
 * Codes:
```
dos2unix extractCoord.sh
chmod +x extractCoord.sh
bash extractCoord.sh [abc_SeqScope_1st.fastq.gz] [abc_SeqScope_2nd_R1.fastq.gz] [20]
```
 * Output:

```
spatialcoordinates.txt: Five columns representing 1st-Seq HDMIs, lane, tile, X, Y 
whitelists.txt: This is the whitelists of HDMIs used for STARsolo alignment. If MiSeq is used for 1st-Seq, then whitelists are the reverse complementary of HDMIs in bottom tiles from 1st-Seq ; if HiSeq is used for 1st-Seq,  whitelists are the reverse complementary of HDMIs in all tiles in lane 2 from 1st-Seq.
HDMI_SeqScope_2nd.txt
```

### Step1b: Discovery plot of tissue boundary
To Visualize the spatial map of HDMI barcode and estimate the tissue boundary, please run estimateTissueBoundary function within the shell. Please install the following python modules before running the script.
* Required python modules
  *  os
  *  sys
  *  numpy
  *  pandas
  *  os
  *  mpl_scatter_density
  *  matplotlib
  *  pylab
  
The script estimateTissueBoundary.py can be found under script folder in this repository. Please download the script in your working directory and run it within the shell.

* Input:
```
[pos1stSeq]:  txt file with spatial information from 1st-Seq. The txt file have five columns representing 1st-Seq HDMIs, lane, tile, X, Y. We can use spatialcoordinates.txt from output of extractCoord.sh
[hdmi2ndSeq]: txt file with HDMIs from the 2nd-Seq. We can use HDMI_SeqScope_2nd.txt from extractCoord.sh
[maxScale]: vmax value for the colorbar; If not known, just put 0 as input. Sometimes outliers exist and make it hard to visulize the tissue boundary. maxScale of colorbar helps with a better visualization
[outpath]: path to output the plots
```
* Code
```
python3 estimateTissueBoundary.py [pos1stSeq] [hdmi2ndSeq] [maxScale] [outpath]
```
* Output
```
tile*.png: The discovery plot can be uesd to compare with H&E images
```

## Step2: STARsolo Alignment and Data Binning
In this subsection, we first preprocess the data and run alignment with reference genome using STARsolo. Then the digital expression matrix (DGE) is binned into square grids with user defined options.

### Step2a: Alignment
This step is to preprocess the fastq files and to align the data to reference genome.The bash script takes several user defined parameters and produces STARsolo summary statistics, and DGE in the current directory. Note: Here we assume you already have the reference genome that is needed for STARsolo alignment. If not please refer to https://hbctraining.github.io/Intro-to-rnaseq-hpc-O2/lessons/03_alignment.html. 
* Input
```
abc_SeqScope_2nd_R1.fastq.gz: Read1 from 2nd-Seq 
abc_SeqScope_2nd_R2.fastq.gz: Read2 from 2nd-Seq
hdmilength: An integer indicating the length of the HDMIs; For now, it can only take 20 or 30. In default, we assume if MiSeq is used for 1st-Seq, then hdmilenght=20; if HiSeq is used for 1st-Seq, then hdmilength=30.
whitelists.txt: whitelists of barcodes from extractCoord.sh
outprefix: prefix for STARsolo output
starpath: path for STAR software
seqtkpath: path for seqtk tool
geneIndex: reference genome directory

```
* Code
```
dos2unix align.sh
chmod +x align.sh
bash align.sh [abc_SeqScope_2nd_R1.fastq.gz] [abc_SeqScope_2nd_R2.fastq.gz] [hdmilength] [whitelists.txt] [outprefix] [starpath] [seqtkpath] [geneIndex]
```
* Output
```
#The STARsolo output that is used for downstream analysis (such as data binning, clustering, cell type mapping)
prefixSolo.out/GeneFull/raw/matrix.mtx
prefixSolo.out/GeneFull/raw/barcodes.mtx
prefixSolo.out/GeneFull/raw/features.mtx

#The STARsolo output that is used for subcellular analysis 
prefixSolo.out/Velocyto/raw/matrix.mtx
prefixSolo.out/Velocyto/raw/barcodes.mtx
prefixSolo.out/Velocyto/raw/features.mtx
```

### Step2b: Data Binning
DGE(prefixSolo.out/GeneFull/raw/) from STARsolo are binned into square grids. In our paper, we tried simple square binning and sliding window binning. Simple square binning generate a super tile with the tiles that the users are insterested in. For sliding window binning, currently it is only available for sub-field of one tile. We would improve this and make updates in the near future.

Please install the following R packages before running the script.
* Required R libraries
  * Rcpp
  * Matrix
  * tictoc
  * ggplot2
  * ggsci
  * Seurat
  * rlist
  * cowplot
  * dplyr
  * mapplots


#### (i) Simple Square Binning
Please download the the script getSimpleGrids.R and collapse.cpp to your working directory and run the command within your shell. 
* Input
```
[seqscope1st]: "MiSeq" or "HiSeq"
[DGEdir]: directory path for DGE from STARsolo output that stores matrix.mtx,features.tsv,barcodes.tsv
[spatial]: txt file that stores spatial coordinates. This file is 'spatialcoordinates.txt' from  extractCoord.sh
[tiles]:a vector of tile numbers that the user is interested in
[nrow]: number of rows for super tile
[ncol]: number of columns for super tile
[sidesize]: side size of squre grids
[outpath]: output directory path
```
* Code
```
Rscript getSimpleGrid.R [seqscope1st] [DGEdir] [spatial] [tiles],[nrow],[ncol],[sidesize],[outpath],[collapsePath])
```
* Output
```
simpleSquareGrids.RDS
(other temporary files)
```
The simpleSquareGrids.RD stores the count matrix, meta information and the spatial coordinates in images slot.

#### (ii) Sliding window binning (sub-field)
Sliding window grids are useful in doing high-resolution cell type mapping. Please download the script getSlidingGrids.R to your working directory and run the command within your shell.
This script collapses HDMIs within a square grids with user-defined grid side length and a sliding window size. In this version, the sliding window function can be only run on a small sub-field. We are currently working on a more flexbile software with the flexibility to make the sliding window function scalable. We will have an update on this page when the package is available.
* Input
```
[seqscope1st]: "MiSeq" or "HiSeq"
[DGEdir]: directory path for DGE from STARsolo output that stores matrix.mtx,features.tsv,barcodes.tsv
[spatial]: txt file that stores spatial coordinates. This file is 'spatialcoordinates.txt' from  extractCoord.sh
[tiles]:a vector of tile numbers that the user is interested in
[nrow]: number of rows for super tile
[ncol]: number of columns for super tile
[sidesize]: side size of squre grids
[outpath]: output directory path
[window]: size of sliding window
[subXmin]: min value of X axis of the sub-field
[subXmax]: max value of X axis of the sub-field
[subYmin]: min value of Y axis of the sub-field
[subYmax]: max value of Y axis of the sub-field

```
* Code
```
```
* Output
```
slidingSquareGrids.RDS
```

After running these steps, the final output files that you will find the most useful will be the following:
- DGE (matrix.mtx, barcodes.tsv, features.tsv) from STARsolo alignment under GeneFull and Velocyto folder
- simpleSquareGrids.RDS
- slidingSquareGrids.RDS

Downstream analysis (clustering,cell type mapping, etc) can be conducted using the three output files. Please refer to Seurat tutorials https://satijalab.org/seurat/archive/v3.2/spatial_vignette.html


### SubCellular Analysis
The DGE (matrix.mtx, barcodes.tsv, features.tsv) for subcellular analysis is under the folder of STARsolo alignment: DGE(xxx/Velocyto/raw/). In this step, we divide genes of spliced and unspliced into three groups and plot them in the sub-micrometer resolution to visulalize the spatial pattern. The script subCellularAna.py can be found under script folder in this repository. Please download subCellularAna.py in your working directory and run this within the shell.

* Required python modules
  *  os
  *  sys
  *  csv
  *  gzip
  *  scipy
  *  numpy
  *  pandas
  *  seaborn
  *  matplotlib
  *  scprep

* Input
```
[DGEdir]: directory of digital expression matrix from STARsolo alignment under Velocyto option (xxx/Velocyto/raw/)
[workingdir]: working directory
[spatial]: spatial coordinates 
[seqscope1st]: string to indicaate if the seqscope1st file is from MiSeq or HiSeq
[tiles]: array of tiles, for example: 2106,2107,2108. Please use comma to separate the multiple tile numbers
[alpha]: transpanrency, if unknown, please set to 0, the default alpha is 0.01
[vmin]: min value for color bar, if unknown, please set to 0, the default vmin is 0
[vmax]: max value for color bar, if unknown, please set to 0, the default vmax is 2
```
* Code
```

python3 subCellularAna [DGEdir],[workingdir],[spatial],[seqscope1st],[tiles],[alpha],[vmin],[vmax]

```
* Output
```
splice1_pos.csv
splice2_pos.csv
splice3_pos.csv
unsplice1_pos.csv
unsplice2_pos.csv
unsplice3_pos.csv'

splice_subset_1_tilexxxx.png
splice_subset_2_tilexxxx.png
splice_subset_3_tilexxxx.png
unsplice_subset_1_tilexxxx.png
unsplice_subset_2_tilexxxx.png
unsplice_subset_3_tilexxxx.png

```






