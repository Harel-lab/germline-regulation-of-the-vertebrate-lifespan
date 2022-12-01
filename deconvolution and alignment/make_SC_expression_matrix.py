###############
# this code makes a gene expression matrix
###############
import os, re
import pandas as pd


orgs = ['kf'] #'mm9']#'hg38'] #,
# basePath = "/cs/icore/tehila_atlan/train/Unaligned_200521_NS500183_0650_AHY33KBGXC/"
# batches = ["Chen_MP1_", "Chen_MP2_", "Chen_MP3_"]
'''
# henrik
basePath = "/cs/icore/tehila_atlan/data/Henrik/"
batches = ["Henrik_SC_E1_", "Henrik_SC_E2_", "Henrik_SC_E3_", "Henrik_SC_C1_", "Henrik_SC_C2_", "Henrik_SC_C3_"]
'''
# Eitan germ
basePath = '/sci/labs/itamarh/tehila_atlan/icore-home/data/scRNAseq/Eitan/'
batches = ['WT_M3_', 'WT_F2_', 'KO_M3_', 'KO_F1_',
           'KilliWT_M1_', 'KilliKO_F1_', 'KilliKO_M3_', 'KilliWT_F2_',
           'KOF3_', 'KOF4_', 'KOM1_', 'KOM2_', 'WTF1_', 'WTF3_', 'WTM2_',
           'KfKOF2_', 'KfKOF3_', 'KfKOM2_', 'KfKOM4_', 'KfWTF1_', 'KfWTF3_', 'KfWTM1_', 'KfWTM2_', 'KfWTM1_'] #]#, ]

# symbols2names=pd.read_csv('/cs/icore/tehila_atlan/shirliyaScripts/UCSCtoGeneSymbol_mm9_knownGene.bed',header=None, index_col=0,delimiter='\t',names=['gene_name'])
# print(symbols2names)
savein = "raw_data6.22/"
os.popen("mkdir -p " + basePath + savein)
# basePath = "/cs/icore/shirliyadadon/data/ido_amit/"
# batch = ["SRR507"]
# files = os.popen('ls ' + basePath + batch[0] + "*.fastq").read().split('\n')[:-1]  # take names of all fastq files of this experiment
# files = [f.split('/')[-1].split('.')[0] for f in files]  # create a list of all the names without the last-empty string
# batches = files[0:len(files):2]
# batches = ['Rand_'+f for f in batches]
for org in orgs:
    for batch in batches:
        path = basePath + batch[:-1] + "/rsem_output/"
        cells = os.popen("ls " + path + "*." + org + ".genes.results").read().split('\n')[:-1]
        pathExpCounts = basePath + savein + "named_RSEM_" + org + "_matrix_expected_counts_" + batch[:-1] + ".csv"
        matExpectedCounts = open(pathExpCounts, 'w')
        for singleCell in cells:
            search = re.search(path + "(.+)\." + org + "\.genes\.results", singleCell)
            matExpectedCounts.write(search.group(1))
            with open(singleCell, 'r') as sc:
                line = sc.readline().strip()
                while line:
                    if "TPM" in line:
                        line = sc.readline().strip()
                        continue
                    lineInParts = line.split('\t')
                    matExpectedCounts.write('\t' + lineInParts[4])
                    line = sc.readline().strip()
                matExpectedCounts.write('\n')
        matExpectedCounts.close()
        names = ""
        with open(cells[0], 'r') as sc:
            line = sc.readline().strip()
            while line:
                if "TPM" in line:
                    line = sc.readline().strip()
                    continue
                lineInParts = line.split('\t')
                names += '\t' + lineInParts[0]
                line = sc.readline().strip()
        names=names.split('\t')[1:]
        # if org == 'mm9':
        #     names = symbols2names.loc[names,'gene_name']
        # print(names)
        names= '\t'+'\t'.join(list(names))
        with open(pathExpCounts, 'r') as matExpectedCounts:
            content = matExpectedCounts.read()
        with open(pathExpCounts, 'w') as matExpectedCounts:
            matExpectedCounts.write(names + '\n')
            matExpectedCounts.write(content)
print ("End Program\n")


###
"""
a = ['nkx1-1', 'nkx1-2', 'nkx2-1', 'nkx2-3', 'nkx2-5', 'nkx2-8', 'nkx3-2', 'nkx6-1', 'nkx6-2']
b = ['nkx1_1', 'nkx1_2', 'nkx2_1', 'nkx2_3', 'nkx2_5', 'nkx2_8', 'nkx3_2', 'nkx6_1', 'nkx6_2']

for f in ['C1', 'C2', 'C3', 'E1', 'E2', 'E3']:
    aaa = pd.read_table('icore-home/data/scRNAseq/Henrik/raw_data1.21/named_RSEM_kf_matrix_expected_counts_Henrik_SC_%s.csv' % f, index_col=0)
    for i in range(0,10):
        aaa.columns = [s.replace('-%d' % i, '_%d' % i) for s in aaa.columns]
        for j in range(0, len(a)):
            aaa.columns = [s.replace(b[j], a[j]) for s in aaa.columns]
    aaa.to_csv('icore-home/data/scRNAseq/Henrik/raw_data1.21/named_RSEM_kf_matrix_expected_counts_Henrik_SC_%s.csv' % f, sep='\t')
"""