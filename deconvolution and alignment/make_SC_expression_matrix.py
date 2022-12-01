###############
# this code makes a gene expression matrix
###############
import os, re
import pandas as pd

# scRNAseq germ
basePath = 'icore-home/data/scRNAseq/'
batches = ['WT_M3_', 'WT_F2_', 'KO_M3_', 'KO_F1_',
           'KilliWT_M1_', 'KilliKO_F1_', 'KilliKO_M3_', 'KilliWT_F2_',
           'KOF3_', 'KOF4_', 'KOM1_', 'KOM2_', 'WTF1_', 'WTF3_', 'WTM2_',
           'KfKOF2_', 'KfKOF3_', 'KfKOM2_', 'KfKOM4_', 'KfWTF1_', 'KfWTF3_', 'KfWTM1_', 'KfWTM2_', 'KfWTM1_'] 


savein = "raw_data/"
os.popen("mkdir -p " + basePath + savein)

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
