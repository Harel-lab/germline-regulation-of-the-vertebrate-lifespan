###############
# this code calls to run of rsem alignment
###############
import os, sys
import re, time


path = 'icore-home/data/scRNAseq/'
""" 
# scRNAseq germ
batch = ['KilliWT_M1_', 'KilliKO_F1_', 'KilliKO_M3_', 'KilliWT_F2_']

# scRNAseq germ2
batch = ['WT_F2_', 'KO_M3_', 'WT_M3_', 'KO_F1_']

# scRNAseq germ3
batch = ['KOF3_', 'KOF4_', 'KOM1_', 'KOM2_', 'WTF1_', 'WTF3_', 'WTM2_']
"""

# scRNAseq germ4
batch = ['KfWTM1_','KfWTF1_', 'KfWTF3_', 'KfKOF2_', 'KfKOF3_', 'KfWTM2_', 'KfKOM2_', 'KfKOM4_', 'KfWTM1_']

for b in batch:
    # newPath = path + "SC_" + b[:-1]
    newPath = path + b[:-1]
    orgMus = "mm9"
    orgH = "GRCh38_RSEM"
    orgK = "Nfu_20140520"
    Fastq = os.popen("ls " + newPath + "/*.fastq").read().split("\n")[:-1]
    dirRsem = newPath + "/rsem_output"
    os.popen("mkdir -p " + dirRsem)
    dirBowtie = newPath + "/bowtie_output"
    # os.popen("mkdir -p " + dirBowtie)
    count =0
#[1:3]
    for fastq in Fastq:
        filePath = re.search(newPath + "\/(.+)\.fastq$", fastq)
        if filePath is not None:
            outRsemMus = dirRsem + "/" + filePath.group(1) + ".mm9"
            outRsemH = dirRsem + "/" + filePath.group(1) + ".hg38"
            outRsemK = dirRsem + "/" + filePath.group(1) + ".kf"
            os.popen(
                "sbatch --mem=150g --time=1-0 --job-name=RSEM" + filePath.group(1)[:2] +
                " -o icore-home/logs/RSEMLogIdo.log --wrap=\"python3 rsem.py "
                + fastq + " " + outRsemMus + " " + outRsemH + " " + outRsemK + " " + orgMus + " " + orgH + " " + orgK + "\"")
            count += 1
            if count == 505:
                time.sleep(60)
                count = 0
print("End program\n")
