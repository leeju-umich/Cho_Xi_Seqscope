import sys
import os
import csv
import gzip
import scipy.io
import numpy as np
import pandas as pd 
import seaborn as sns
import matplotlib.pyplot as plt
import matplotlib.colors as colors
import scprep
from sys import argv

def subPlot(df_pos,tiles,alpha,vmin,vmax,title):
    """Function to plot unspliced or spliced genes on subcelullar resolution"""
    for i in tiles:
        #print(i)
        print(i)
        print(df_pos.head())
        z = df_pos[df_pos.tile_miseq==int(i)]
        print(z.head())
        maxV = max(z['umi'])
        minV = 0
        centerV = 0.5
        cmap = colors.LinearSegmentedColormap.from_list('', ['white','grey','black'])
        norm = colors.DivergingNorm(vmin=minV, vcenter=centerV, vmax=maxV)

        fig, ax = plt.subplots()
        fig.gca().set_aspect('equal', adjustable='box')

        points = ax.scatter(z['y_miseq'],z['x_miseq'], c=z['umi'], s=1 ,alpha=0.01,vmin=0,vmax=2,cmap = cmap,norm=norm)  #tune alpha and size parameter here
        cbar = fig.colorbar(points, ax=ax,extend = 'max')
        plt.title('tile' + str(i))
        #plt.show()
        plt.savefig(title+"tile"+str(i)+".png", dpi=500)
        plt.show()

def getSubset(unspliced,spliced,geneInd,barcode_df,bottom):
    """Function to subset unspliced or spliced data"""

    unspliced_geneset_df=unspliced.loc[unspliced['gene'].isin(geneInd)]
    spliced_geneset_df = spliced.loc[spliced['gene'].isin(geneInd)]
    unsplicedUMI_geneset = unspliced_geneset_df[['barcode','umi']].groupby(['barcode']).sum().reset_index()
    splicedUMI_geneset = spliced_geneset_df[['barcode','umi']].groupby(['barcode']).sum().reset_index()
    df_unspliced=barcode_df.merge(unsplicedUMI_geneset,how='inner',left_on="barcodeInd",right_on="barcode")
    df_spliced=barcode_df.merge(splicedUMI_geneset,how='inner',left_on="barcodeInd",right_on="barcode")
    unspliced_pos = pd.merge(bottom, df_unspliced, how="right",left_on='HDMI', right_on='BARCODE')
    spliced_pos = pd.merge(bottom, df_spliced, how="right",left_on='HDMI', right_on='BARCODE')
    return unspliced_pos,spliced_pos

def subCellularAna(DGEdir,workingdir,spatial,seqscope1st,tiles,alpha,vmin,vmax):
    """Function to plot unspliced or spliced genes on subcelullar resolution"""
    if os.path.isdir(DGEdir)==True:
        os.chdir(DGEdir)
    else:
        raise NameError('Not valide DGE directory')
    if os.path.isdir(workingdir)==False:
        raise NameError('Not valide working directory')
    if os.path.isfile(spatial)==False:
        raise NameError('Not valide spatial file')
    if alpha is None:
        alpha=0.01
    if vmin is None:
        vmin=0
    if vmax is None:
        vmax=2



    miseq_pos = pd.read_csv(spatial,delim_whitespace=True, header=None)
    miseq_pos.columns = ['HDMI','lane_miseq','tile_miseq','x_miseq','y_miseq']
    tiles_cat = '|'.join(tiles)
    print(tiles)
    print(tiles_cat)
#    tiles_cat ='|'.join([2106,2107])
    miseq_pos[miseq_pos['tile_miseq'].astype(str).str.contains(tiles_cat)]
    if seqscope1st=='MiSeq':
        bottom= miseq_pos[miseq_pos['tile_miseq']>=2000]    #This should be made more flexible
    if seqscope1st=='HiSeq':
        bottom=miseq_pos
        
    os.chdir(DGEdir) #Velocyto/raw/
    #Read featrures.tsv
    gene_names = [row[1] for row in csv.reader(open("features.tsv"), delimiter="\t")]
    gene_names = np.array(gene_names)
    #Read barcodes.tsv
    bc = []
    with open(r"barcodes.tsv") as f:
        reader = csv.reader(f, delimiter='\t', quotechar='"')
        for row in reader:
            if row:
                bc.append(row[0])
    #Read  matrix.mtx and create splice and unsplice matrix
    m = pd.read_csv(r"matrix.mtx", sep=' ',skiprows=[0,1,2],header=None)
    print(m.shape)
    splice=m[[0,1,2]]
    splice.columns = ['gene','barcode','umi']
    print(splice.shape)
    unsplice = m[[0,1,3]]
    unsplice.columns = ['gene','barcode','umi']
    print(unsplice.shape)

    #randomly divide the genes to three sets
    #set=[1,2,3]
    rng = np.random.default_rng()
    geneset_ind = np.array(np.arange(len(gene_names))+1)
    rng.shuffle(geneset_ind)
    print(geneset_ind.shape)
    geneset_ind = np.array_split(geneset_ind,3)  #instead of np.split()
    geneset1_ind = geneset_ind[0]
    geneset2_ind = geneset_ind[1]
    geneset3_ind = geneset_ind[2]
    barcode_df = pd.DataFrame({"BARCODE":bc,"barcodeInd":range(1,(len(bc)+1))})

    unspliced1_pos,spliced1_pos = getSubset(unsplice,splice,geneset1_ind,barcode_df,bottom)
    unspliced2_pos,spliced2_pos = getSubset(unsplice,splice,geneset1_ind,barcode_df,bottom)
    unspliced3_pos,spliced3_pos = getSubset(unsplice,splice,geneset1_ind,barcode_df,bottom)

    os.chdir(workingdir)
    for i in range(1,4):
        file1_name = 'unspliced'+str(i)+'_pos'
        file2_name = 'spliced'+str(i)+'_pos'
        out1_name = 'unsplice'+str(i)+'_pos.csv'
        out2_name = 'splice'+str(i)+'_pos.csv'
        plot1_name = 'unsplice_subset_'+str(i)
        plot2_name = 'splice_subset_'+str(i)
        eval( file1_name).to_csv(out1_name)
        eval(file2_name).to_csv(out2_name)
        subPlot(eval(file1_name),tiles,alpha,vmin,vmax,plot1_name)
        subPlot(eval(file2_name),tiles,alpha,vmin,vmax,plot2_name)




DGEdir=sys.argv[1] 
workingdir=sys.argv[2] 
spatial=sys.argv[3] 
seqscope1st=sys.argv[4] 
tiles = sys.argv[5].split(',')
alpha=sys.argv[6]
vmin=sys.argv[7]
vmax=sys.argv[8]
subCellularAna(DGEdir,workingdir,spatial,seqscope1st,tiles,alpha,vmin,vmax)

