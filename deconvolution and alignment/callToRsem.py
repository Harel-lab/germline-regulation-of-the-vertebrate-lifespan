###############
# this code calls to run of rsem alignment
###############
import os, sys
import re, time

""" Henrik- senescence
path = "/cs/icore/tehila_atlan/data/Henrik/"
batch = ["Henrik_SC_E1_", "Henrik_SC_E2_", "Henrik_SC_E3_", "Henrik_SC_C1_", "Henrik_SC_C2_", "Henrik_SC_C3_"]#,  # ["pertRANDOMminus2i_"]  # ["Chen_ES_"]

# Eitan germ
path = '/sci/labs/itamarh/tehila_atlan/icore-home/data/scRNAseq/Eitan/'
batch = ['KilliWT_M1_', 'KilliKO_F1_', 'KilliKO_M3_', 'KilliWT_F2_']

# Eitan germ2
path = '/sci/labs/itamarh/tehila_atlan/icore-home/data/scRNAseq/Eitan/'
batch = ['WT_F2_', 'KO_M3_', 'WT_M3_', 'KO_F1_']

# Eitan germ3
path = '/sci/labs/itamarh/tehila_atlan/icore-home/data/scRNAseq/Eitan/'
batch = ['KOF3_', 'KOF4_', 'KOM1_', 'KOM2_', 'WTF1_', 'WTF3_', 'WTM2_']
"""

# Eitan germ4
path = '/sci/labs/itamarh/tehila_atlan/icore-home/data/scRNAseq/Eitan/'
batch = ['KfWTM1_']
#'KfWTF1_', 'KfWTF3_', 'KfKOF2_', 'KfKOF3_', 'KfWTM2_', 'KfKOM2_', 'KfKOM4_', 'KfWTM1_'






# path = "/cs/icore/shirliyadadon/data/ido_amit/"
# batch = ["SRR507"]
# files = os.popen('ls ' + path + batch[0] + "*.fastq").read().split('\n')[:-1]  # take names of all fastq files of this experiment
# files = [f.split('/')[-1].split('.')[0] for f in files]  # create a list of all the names without the last-empty string
# randFiles = files[0:len(files):2]asxa
# randFiles = ['Rand_'+f for f in randFiles]
# path = sys.argv[1]
# batch = sys.argv[2:]

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
            # print("starting: " + filePath.group(1)+"\n")
            # print("sbatch --mem=150g --time=1-0 --job-name=RSEM" + filePath.group(1)[:2] +
            #     " -o /sci/labs/itamarh/tehila_atlan/icore-home/logs/RSEMLogIdo.log --wrap=\"python3 /sci/labs/itamarh/tehila_atlan/icore-home/shirliyaScripts/rsem.py "
            #     + fastq + " " + outRsemMus + " " + outRsemH + " " + outRsemK + " " + orgMus + " " + orgH + " " + orgK + "\"")
            print(fastq)
            os.popen(
                "sbatch --mem=150g --time=1-0 --job-name=RSEM" + filePath.group(1)[:2] +
                " -o /sci/labs/itamarh/tehila_atlan/icore-home/logs/RSEMLogIdo.log --wrap=\"python3 /sci/labs/itamarh/tehila_atlan/icore-home/shirliyaScripts/rsem.py "
                + fastq + " " + outRsemMus + " " + outRsemH + " " + outRsemK + " " + orgMus + " " + orgH + " " + orgK + "\"")
            count += 1
            if count == 505:
                time.sleep(60)
                count = 0
print("End program\n")
