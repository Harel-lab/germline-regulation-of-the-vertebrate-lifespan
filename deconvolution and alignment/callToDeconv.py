# this script calls to a deconvolution of raw files script
# comment out non relevant type of experiment
import os, sys

exp = "indrop"
path = '/icore-home/data/scRNAseq/fastq/'
'''
# scRNAseq germ1
batch = ['KilliWT_M1_']#'KilliKO_F1_', 'KilliKO_M3_', 'KilliWT_F2_',

# scRNAseq germ2
batch = ['KO_M3_'] #  , 'KO_F1_', 'WT_F2_', 'WT_M3_']

# scRNAseq germ3
batch = ['KOF3_', 'KOF4_', 'KOM1_', 'KOM2_', 'WTF1_', 'WTF3_', 'WTM2_']
'''
# scRNAseq germ4
batch = ['KfKOF2_', 'KfKOF3_', 'KfKOM2_', 'KfKOM4_', 'KfWTF1_', 'KfWTF3_', 'KfWTM1_', 'KfWTM2_']


if exp == "indrop":
    for b in batch:
        os.popen(
            "sbatch --mem=50g --time=3-0 --job-name=" + b + "_DeconvInDrop" + " -o icore-home/logs/DeconvInDropLog.log --wrap=\"python "
                                                                              "Deconvoluting_inDropMP.py " + path + " " + b + "\"")
