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

org = "killifish" 


if org == "killifish":
    rsem_line = "rsem-calculate-expression --num-threads 1 --bowtie2 --single-cell-prior --estimate-rspd " \
                + fastq + " icore-home/RefGenomes/" + orgK + "/RefSeq_RSEM " + outRsemK


    print(rsem_line)
    log = os.popen(rsem_line).read()
    os.popen(rsem_line)
print("End program")
