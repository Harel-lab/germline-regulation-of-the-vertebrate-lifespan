# this script calls to a deconvolution of raw files script
# comment out non relevant type of experiment
import os, sys

exp = "indrop"
'''
# Henrik senescence
path = "/cs/icore/oren.ram/SequencingData/Unaligned_200708_NS500183_0670_AHY35JBGXF/Henrik_SC/"
batch = ["Henrik_SC_C1_", "Henrik_SC_C2_", "Henrik_SC_C3_", "Henrik_SC_E1_", "Henrik_SC_E2_", "Henrik_SC_E3_"]

# Eitan germ1
path = '/cs/icore/oren.ram/SequencingData/Unaligned_210202_NS500183_0747_AHHKJMBGXH/'
batch = ['KilliWT_M1_']#'KilliKO_F1_', 'KilliKO_M3_', 'KilliWT_F2_',

# Eitan germ2
path = '/sci/labs/itamarh/tehila_atlan/icore-home/data/scRNAseq/Eitan/fastq/2/'
batch = ['KO_M3_'] #  , 'KO_F1_', 'WT_F2_', 'WT_M3_']

# Eitan germ3
path = '/sci/labs/itamarh/tehila_atlan/icore-home/data/scRNAseq/Eitan/fastq/3/'
batch = ['KOF3_', 'KOF4_', 'KOM1_', 'KOM2_', 'WTF1_', 'WTF3_', 'WTM2_']
'''
# Eitan germ4
path = '/sci/labs/itamarh/tehila_atlan/icore-home/data/scRNAseq/Eitan/fastq/4/'
batch = ['KfKOF2_', 'KfKOF3_', 'KfKOM2_', 'KfKOM4_', 'KfWTF1_', 'KfWTF3_', 'KfWTM1_', 'KfWTM2_']


if exp == "indrop":
    for b in batch:
        os.popen(
            "sbatch --mem=50g --time=3-0 --job-name=" + b + "_DeconvInDrop" + " -o /sci/labs/itamarh/tehila_atlan/icore-home/logs/DeconvInDropLog.log --wrap=\"python "
                                                                              "/sci/labs/itamarh/tehila_atlan/icore-home/shirliyaScripts/Deconvoluting_inDropMP.py " + path + " " + b + "\"")
