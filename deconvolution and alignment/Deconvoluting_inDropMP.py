####################
# works on fastq files
# this script take multiple files with multiple cells and creates each cell a specific file
# need arguments where the fastqs are- 1.path 2.batch(the start string of all files from the same experiment) 3.UGI barcodes file
####################
import os, sys, datetime
import re

path = sys.argv[1]
batch = sys.argv[2]
# newPath = path + "SC_" + batch[:-1]
newPath = "icore-home/data/scRNAseq/" + batch[:-1] # + "2"
os.popen("mkdir -p " + newPath)  # create a folder for this experiment
print("start of " + newPath)
print(datetime.datetime.now())
files = os.popen('ls ' + path + batch + "*.fastq").read()  # take names of all fastq files of this experiment
files = files.split('\n')[:-1]  # create a list of all the names without the last-empty string
log = open(newPath + "/Deconv_log.log", 'w')

minNumOfTransc = 3000  # threshold for minimum number of transcripts per cell
numUmiNoFilter = 0
numUmiWithFilt = 0
numAllReads = 0
#################
# create dictionaries
eachCellAllRead2 = {}
cellBCDict = {}
cellMultiPert = {}
noMultiPert = {}
uniqueReads = set()
#################
# start iterating through files
for i in range(0, len(files), 2):
    #################
    # OPEN FILES
    allReadOne = open(files[i], 'r')
    print("opened " + files[i])
    allReadTwo = open(files[i + 1], 'r')
    print("opened " + files[i + 1])
    # readsToFilter = files[i + 1].split('.fas')[0] + '_readsToBothOrgs.txt'
    #################
    line = "start"
    while line:  # while file not ended
        #################
        # extract one read (4 lines) from read1 file-no need after finding cell BC
        line1Read1 = allReadOne.readline()[:-1]  # go a line ahead
        line2Read1 = allReadOne.readline()[:-1]  # wanted line of seq from read1 file
        line1Read1 = allReadOne.readline()  # go more lines ahead. don't save, garbage
        line1Read1 = allReadOne.readline()  # go more lines ahead. don't save, garbage
        # extract one read (4 lines) from read2 file-in order to construct a cell fastq file
        line1Read2 = allReadTwo.readline()[:-1]  # wanted line of id
        line2Read2 = allReadTwo.readline()[:-1]  # wanted line of seq from read1 file
        line3Read2 = allReadTwo.readline()[:-1]
        line4Read2 = allReadTwo.readline()[:-1]
        line = line4Read2
        #################
        # search pattern of cell BC and UMI
        findCBC = re.search("^(\w+)GA.T.A.+C.CC.T(\w{8})(\w{5})", line2Read1)
        if findCBC is not None:
            cellbc = findCBC.group(1) + findCBC.group(2)  # connect cell barcodes, separated with 22bp
            umi = findCBC.group(3)  # extract UMI
            if cellbc in cellBCDict:  # add count of reads for this cell, if exists or not
                cellBCDict[cellbc] += 1
            else:
                cellBCDict[cellbc] = 1
            if cellbc not in eachCellAllRead2:  # if the key doesnt exist, initialize it
                eachCellAllRead2[cellbc] = ""
            # recreate read2 and take unique from each

            line2Read2 = line2Read2[6:]
            line4Read2 = line4Read2[6:]
            if re.search("N", line2Read2[:5]) is not None:
                continue
            # make an unique seq for each read to eliminate PCR duplicates
            # include cellBC to avoid match to other cells (UMIs are overlapping between cells)
            unique = cellbc + umi + line2Read2[:5]
            numUmiNoFilter += 1
            if unique not in uniqueReads:
                read = ""
                numUmiWithFilt += 1
                read += line1Read2 + ":" + umi + "\n" + line2Read2 + "\n" + line3Read2 + "\n" + line4Read2 + "\n"
                eachCellAllRead2[cellbc] += read
                uniqueReads.add(unique)
        numAllReads += 1  # just a check
    allReadOne.close()
    allReadTwo.close()

# check quality
over100 = 0
over1000 = 0
over3000 = 0
over5000 = 0
over10000 = 0

for k, v in list(cellBCDict.items()):
    if v >= 100:
        over100 += 1
    if v >= 1000:
        over1000 += 1
    if v >= 5000:
        over5000 += 1
    if v >= 10000:
        over10000 += 1
    if v >= minNumOfTransc:
        over3000 += 1
    else:
        del cellBCDict[k]

#################
# write logs of run
log.write("Number of reads: " + str(numAllReads))
log.write("\nNumber of significant unique cell barcode: " + str(len(cellBCDict)))

# table = open(newPath + "/tablePert.csv", 'w')
# table.write("CELL,PERTURBATION,NUMBER,IN-LIST\n")
# count = 0
for k in cellBCDict.keys():
    if k in eachCellAllRead2:
        with open(newPath + "/" + k + ".fastq", 'w') as singleCell:
            singleCell.write(eachCellAllRead2[k])

log.write("\nNumber of UMIs without filtering: " + str(numUmiNoFilter))
log.write("\nNumber of UMIs with filtering: " + str(numUmiWithFilt))
log.write("\nNumber of cells with over 100 transcripts: " + str(over100))
log.write("\nNumber of cells with over 1000 transcripts: " + str(over1000))
log.write("\nNumber of cells with over %d (threshold) transcripts: %d" % (minNumOfTransc, over3000))
log.write("\nNumber of cells with over 5000 transcripts: " + str(over5000))
log.write("\nNumber of cells with over 10000 transcripts: " + str(over10000))
# table.write("number of hits: " + str(count))
########
# CLOSE FILES
log.close()
# table.close()
print("end of " + newPath)
print(datetime.datetime.now())
