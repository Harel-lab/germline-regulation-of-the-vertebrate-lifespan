import pandas as pd
import numpy as np
import os
import seaborn as sns
import matplotlib.pylab as plt
import cmasher as cmr

path = '/ZFcomparision/'
method = '' 

convertKF2H = pd.read_csv('/data/genes_names7.csv')
convertZF2H = pd.read_csv(os.path.join(path, 'mart_export.txt'))

print(convertZF2H)


def convertGenesNamesKF2H(geneList, source='NCBI', target='Final symbol'):
    if type(geneList) is not list:
        geneList = geneList.to_list()
    return convertKF2H.set_index(source).loc[geneList].reset_index()[target].to_list()


def convertGenesNamesZF2H(geneList, source='Zebrafish gene name', target='Gene name'):
    if type(geneList) is not list:
        geneList = geneList.to_list()
    return convertZF2H.set_index(source).loc[geneList].reset_index()[target].to_list()


maleZF = pd.read_csv(os.path.join(path, 'maleZf.csv'))
print(len(set(maleZF['gene'])))
maleZF = maleZF.merge(convertZF2H[['Zebrafish gene name', 'Gene name']], left_on='gene', right_on='Zebrafish gene name', how='inner')
print(len(set(maleZF['gene'])))
maleClusters = ['late round spermatid', 'spermatocyte', 'elongated spermatid', 'SPG', 'early round spermatid', 'SPG2',
                'middle round spermatid', 'Sertoli', 'erythrocyte', 'Leydig']
femaleZF = pd.read_csv(os.path.join(path, 'FemaleZF.csv'))
print(len(set(femaleZF['gene'])))
femaleZF = femaleZF.merge(convertZF2H[['Zebrafish gene name', 'Gene name']], left_on='gene', right_on='Zebrafish gene name', how='inner')
print(len(set(femaleZF['gene'])))
femaleZF.rename(columns={"avg_logFC": "avg_log2FC"}, inplace=True)
femaleClusters = set(femaleZF['cluster'])
print(femaleZF[femaleZF['p_val'] > 0.01])     
print(femaleZF['p_val'].max())    

print(maleZF)
print(femaleZF)

somaticClusters = ['Sertoli', 'Fibroblast', 'Neutophils', 'Granulosa', 'SMC', 'Mono I', 'Mono II', 'leydig', 'mixed',
                  'OE', 'Erythroid', 'Theca']
GermClusters = ['c','g','a','e','d','b','f']

# ,names,scores,pvals,pvals_adj,logFC,NCBI,Human,NCBI Definition

def calculateCorrTable(ZFtable, KFgroup, KFclusters, indexZF, ZFclusters, nKF):
    corrTable = pd.DataFrame(index=indexZF, columns=range(0, nKF))

    for i in range(0, nKF):
        kfCluster = pd.read_csv(os.path.join(path, 'killifish', 'WT%s%dleiden.csv' % (KFgroup, i)))
        # kfCluster = kfCluster[(kfCluster['logFC'] > 0)]
        # kfCluster = kfCluster[(kfCluster['logFC'] > 0.25) & (kfCluster['pvals'] < 0.01)]
        for j in indexZF:
            zfCluster = ZFtable[ZFtable['cluster'] == j]
            a = kfCluster.merge(zfCluster, left_on='Human', right_on='Gene name')
            corrTable.loc[j, i] = np.corrcoef(a['logFC'], a['avg_log2FC'])[0, 1]

    print()
    corrTable = corrTable[corrTable.columns].astype(float)
    corrTable.index = ZFclusters
    corrTable.columns = KFclusters
    corrTable.sort_index(axis=1, inplace=True)
    corrTable.sort_index(axis=0, inplace=True)
    return corrTable


cmap = cmr.get_sub_cmap('Blues', 0.0, 1.0)

corrTable = calculateCorrTable(maleZF, 'germ', GermClusters, range(0, len(maleClusters)), maleClusters, 7)
corrTable = corrTable.loc[['SPG', 'SPG2', 'spermatocyte', 'early round spermatid', 'middle round spermatid', 'late round spermatid', 'elongated spermatid']]
print(corrTable)
plt.figure()
hmap = sns.heatmap(corrTable, cmap=cmap)
plt.savefig(os.path.join(path, method + 'GermVsMaleZF2.pdf'), bbox_inches='tight')

# ----- compare the marker genes in ZF papers
germMarkers = pd.read_csv(os.path.join(path, '../markersGerm.csv'))
print(germMarkers['Final Symbol'])

humanMarkers = convertGenesNamesKF2H(germMarkers['Final Symbol'], 'Final symbol', 'Human')
print(maleZF[maleZF['Gene name'].isin(humanMarkers)])

germMarkers = pd.read_csv(os.path.join(path, '../markersWTcellTypes.csv'))
humanMarkers = convertGenesNamesKF2H(germMarkers['Final Symbol'], 'Final symbol', 'Human')
print(femaleZF[femaleZF['Gene name'].isin(humanMarkers)])

