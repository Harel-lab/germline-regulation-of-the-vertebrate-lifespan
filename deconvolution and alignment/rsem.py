#############
# this code runs bowtie2 through rsem
#############
import os, sys

fastq = sys.argv[1]
outRsemMus = sys.argv[2]
outRsemH = sys.argv[3]
outRsemK = sys.argv[4]
orgMus = sys.argv[5]
orgH = sys.argv[6]
orgK = sys.argv[7]

org = "killifish" #"mouse" #"human"  #

if org == "mouse":
    os.popen(
        "/cs/cbio/tommy/Software/RSEM/rsem-calculate-expression --num-threads 1 --bowtie2"
        " --single-cell-prior --estimate-rspd " + fastq + " /cs/icore/tehila_atlan/RefGenomes/"
        + orgMus + "/RSEM/RefSeq_RSEM " + outRsemMus)
elif org == "human":
    os.popen(
        "/cs/cbio/tommy/Software/RSEM/rsem-calculate-expression --num-threads 1 --bowtie2"
        " --single-cell-prior --estimate-rspd " + fastq + " /cs/icore/tehila_atlan/RefGenomes/"
        + orgH + "/RefSeq_RSEM " + outRsemH)
elif org == "killifish":
    rsem_line = "rsem-calculate-expression --num-threads 1 --bowtie2 --single-cell-prior --estimate-rspd " \
                + fastq + " /sci/labs/itamarh/tehila_atlan/icore-home/RefGenomes/" + orgK + "/RefSeq_RSEM " + outRsemK


    #   /sci/labs/itamarh/tehila_atlan/icore-home/Programs/RSEM/RSEM-1.3.3/
    # rsem_line = 'rsem-calculate-expression --num-threads 1 --bowtie2 --single-cell-prior --estimate-rspd /sci/labs/itamarh/tehila_atlan/icore-home/data/scRNAseq/Eitan/KO_F1/ACGCTTTAGATGGGT.fastq /sci/labs/itamarh/tehila_atlan/icore-home/RefGenomes/Nfu_20140520/RefSeq_RSEM /sci/labs/itamarh/tehila_atlan/icore-home/data/scRNAseq/Eitan/KO_F1/rsem_output/ACGCTTTAGATGGGT.kf'
    print(rsem_line)
    log = os.popen(rsem_line).read()
    os.popen(rsem_line)
    # print(log)
print("End program")


#/cs/cbio/tommy/Software/RSEM/