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
        z = df_pos[df_pos.tile_miseq.eq(i)]
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


    # df_pos = unsplice1_pos
    # alpha=0.01,vmin=0,vmax=2
def subCellularAna(DGEdir,workingdir,spatial,seqscope1st,tiles,alpha,vmin,vmax):
    """Function to plot unspliced or spliced genes on subcelullar resolution"""
    if os.path.isdir(path)==True:
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
    if vmax=None:
        vmax=2



    miseq_pos = pd.read_csv(spatial,delim_whitespace=True, header=None)
    miseq_pos.columns = ['HDMI','tile_miseq','lane_miseq','x_miseq','y_miseq']

    if seqscope1st=='MiSeq':
        bottom= miseq_pos[miseq_pos['tile_miseq']>=2000]
    if seqscope1st=='HiSeq':
        bottom=miseq_pos

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
    rng = np.random.default_rng()
    geneset_ind = np.array(np.arange(len(gene_names))+1)
    rng.shuffle(geneset_ind)
    print(geneset_ind.shape)
    geneset_ind = np.split(geneset_ind,3) 
    geneset1_ind = geneset_ind[0]
    geneset2_ind = geneset_ind[1]
    geneset3_ind = geneset_ind[2]

    unspliced_geneset1_df = unsplice.loc[unsplice['gene'].isin(geneset1_ind)]
    unspliced_geneset2_df = unsplice.loc[unsplice['gene'].isin(geneset2_ind)]
    unspliced_geneset3_df = unsplice.loc[unsplice['gene'].isin(geneset3_ind)]
    spliced_geneset1_df = splice.loc[splice['gene'].isin(geneset1_ind)]
    spliced_geneset2_df = splice.loc[splice['gene'].isin(geneset2_ind)]
    spliced_geneset3_df = splice.loc[splice['gene'].isin(geneset3_ind)]

    #Total UMI per HDMI for splice and unspliced separately
    spliceUMI_geneset1 = spliced_geneset1_df[['barcode','umi']].groupby(['barcode']).sum().reset_index()
    spliceUMI_geneset2 = spliced_geneset2_df[['barcode','umi']].groupby(['barcode']).sum().reset_index()
    spliceUMI_geneset3 = spliced_geneset3_df[['barcode','umi']].groupby(['barcode']).sum().reset_index()
    unspliceUMI_geneset1 = unspliced_geneset1_df[['barcode','umi']].groupby(['barcode']).sum().reset_index()
    unspliceUMI_geneset2 = unspliced_geneset2_df[['barcode','umi']].groupby(['barcode']).sum().reset_index()
    unspliceUMI_geneset3 = unspliced_geneset3_df[['barcode','umi']].groupby(['barcode']).sum().reset_index()
    barcode_df = pd.DataFrame({"BARCODE":bc,"barcodeInd":range(1,(len(bc)+1))})

    df_splice1=barcode_df.merge(spliceUMI_geneset1,how='inner',left_on="barcodeInd",right_on="barcode")
    df_splice2=barcode_df.merge(spliceUMI_geneset2,how='inner',left_on="barcodeInd",right_on="barcode")
    df_splice3=barcode_df.merge(spliceUMI_geneset3,how='inner',left_on="barcodeInd",right_on="barcode")

    df_unsplice1=barcode_df.merge(unspliceUMI_geneset1,how='inner',left_on="barcodeInd",right_on="barcode")
    df_unsplice2=barcode_df.merge(unspliceUMI_geneset2,how='inner',left_on="barcodeInd",right_on="barcode")
    df_unsplice3=barcode_df.merge(unspliceUMI_geneset3,how='inner',left_on="barcodeInd",right_on="barcode")



    splice1_pos = pd.merge(bottom, df_splice1, how="right",left_on='HDMI', right_on='BARCODE')
    splice2_pos = pd.merge(bottom, df_splice2, how="right",left_on='HDMI', right_on='BARCODE')
    splice3_pos = pd.merge(bottom, df_splice3, how="right",left_on='HDMI', right_on='BARCODE')

    unsplice1_pos = pd.merge(bottom, df_unsplice1, how="right",left_on='HDMI', right_on='BARCODE')
    unsplice2_pos = pd.merge(bottom, df_unsplice2, how="right",left_on='HDMI', right_on='BARCODE')
    unsplice3_pos = pd.merge(bottom, df_unsplice3, how="right",left_on='HDMI', right_on='BARCODE')

    os.chdir(workingdir)

    pd.to_csv(splice1_pos,'splice1_pos.csv')
    pd.to_csv(splice2_pos,'splice2_pos.csv')
    pd.to_csv(splice3_pos,'splice3_pos.csv')
    pd.to_csv(unsplice1_pos,'unsplice1_pos.csv')
    pd.to_csv(unsplice2_pos,'unsplice2_pos.csv')
    pd.to_csv(unsplice3_pos,'unsplice3_pos.csv')

    subplot(splice1_pos,tiles,alpha,vmin,vmax,'splice_subset_1')
    subplot(splice1_pos,tiles,alpha,vmin,vmax,'splice_subset_2')
    subplot(splice1_pos,tiles,alpha,vmin,vmax,'splice_subset_3')
    subplot(unsplice1_pos,tiles,alpha,vmin,vmax,'unsplice_subset_1')
    subplot(unsplice2s_pos,tiles,alpha,vmin,vmax,'unsplice_subset_2')
    subplot(unsplice3_pos,tiles,alpha,vmin,vmax,'unsplice_subset_3')



DGEdir=sys.argv[1] 
workingdir=sys.argv[2] 
spatial=sys.argv[3] 
seqscope1st=sys.argv[4] 
tiles=sys.argv[5]   #how to pass a vector of integer
alpha=sys.argv[6]
vmin=sys.argv[7]
vmax=sys.argv[8]
# arr = sys.argv[1].split(',')
# print(arr[2])
#tiles

subCellularAna(DGEdir,workingdir,spatial,seqscope1st,tiles,alpha,vmin,vmax):
