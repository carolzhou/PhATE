################################################################
#
# phate_genecallPhage.py
#
# Programmer: Jeff Kimbrel
#
# Description: Single command to run PhATE, Prodigal, Glimmer and GeneMarkS on a fasta file
#
# Usage:  /usr/local/bin/python3.4 annotatePhage.py fastaFile.fa outFolder
# Run this code on mpath-dev, due to using Python3.4.
#
# Notes:
#    1) Jeff's code was modified to run under python2.7 so that it can be run on mpath,
#       without my having to install a bunch of stuff on mpath-dev. Only print statements  
#       needed to be reformatted.  (cez 23 march 2017)
#    2) The call to Prodigal was modified so that Prodigal would generate the correct
#       .sco file. Now with the file in the correct format, the call to CGC_parser.py 
#       is working, and CGC produces the correct output. 
#
# Further development (cez):
#    23 March 2017: Modifications for incorporation into phate pipeline
#    13 April 2017: Modified paths as environment vars controlled by pipeline driver
#    18 April 2017: Slight modification to calling of CGC
#    25 May   2017: Inserting a "band-aid" to accommodate an error in PHATE
#    15 Aug   2017: Changing PHATE to THEA
#    16 Aug   2017: THEA v.1.0 now outputs contig name in last column; capturing this info
#
################################################################

# This code was developed by Jeff Kimbrel at Lawrence Livermore National Laboratory.
# THIS CODE IS COVERED BY THE BSD LICENSE. SEE INCLUDED FILE BSD.PDF FOR DETAILS.

import os
import sys
import re

#DEBUG = False 
DEBUG = True

# booleans control which gene finder(s) to call
GENEMARKS_CALLS = False   
PRODIGAL_CALLS  = False 
GLIMMER_CALLS   = False 
THEA_CALLS      = False 

# patterns
p_genemarks = re.compile('[gG][eE][nN][eE][mM][aA][rR][kK]')
p_glimmer   = re.compile('[gG][lL][iI][mM][mM][eE][rR]')
p_prodigal  = re.compile('[pP][rR][oO][dD][iI][gG][aA][lL]')
p_thea      = re.compile('[tT][hH][eE][aA]')

########## HOUSEKEEPING ##########

#paths
prodigalPath  = os.environ["PRODIGAL_PATH"]
glimmerPath   = os.environ["GLIMMER_PATH"] 
geneMarkSPath = os.environ["GENEMARKS_PATH"]
theaPath      = os.environ["THEA_PATH"]
cgcPath       = os.environ["CGC_PATH"]

# Data Structures

files = {
    'Raw Prodigal GFF'    : '',
    'Raw Glimmer Output'  : '',
    'Raw GeneMarkS GFF'   : '',
    }

# Input Parameters

print "There are", len(sys.argv), "input parameters:", sys.argv

if len(sys.argv) == 1:
    print ("Usage: /usr/local/bin/python3.4 annotatePhage.py fastaFile.fa outFolder")
    exit(0)

fastaFileName = sys.argv[1]
outputFolder = sys.argv[2] + "/"

# booleans to control gene finding
if len(sys.argv) == 4:
    # get instructions for which gene finders to run; typically, running in the automated pipeline here
    # argument is a string containing the names of gene callers to be run
    genecallParams = sys.argv[3]  # a string listing gene callers to use 
    match_genemarks = re.search(p_genemarks,genecallParams)
    match_prodigal  = re.search(p_prodigal, genecallParams)
    match_glimmer   = re.search(p_glimmer,  genecallParams)
    match_thea      = re.search(p_thea,     genecallParams)
    if match_genemarks:
        GENEMARKS_CALLS = True 
    if match_prodigal:
        PRODIGAL_CALLS = True
    if match_glimmer:
        GLIMMER_CALLS = True 
    if match_thea:
        THEA_CALLS = True 

#logfile = open("./phate_genecallPhage.log","w")
logfilefullpath = outputFolder + "phate_genecallPhage.log"
logfile = open(logfilefullpath,"w")
logfile.write("%s%s\n" % ("Input parameters are:",sys.argv))
workingFolder = os.getcwd()
print "workingFolder is", workingFolder
if not os.path.exists(outputFolder):
    os.mkdir(outputFolder)
logfile.write("%s%s\n" % ("output folder is ", outputFolder))
logfile.write("%s%s\n" % ("working folder is ", workingFolder))
resultsFile = outputFolder + "results.txt"
results = open(resultsFile,"w")
files = {'Results File' : resultsFile}  #*** CHECK THIS
logfile.write("%s%s\n" % ("results file is ",resultsFile))
logfile.write("%s%s\n" % ("GENEMARKS_CALLS is ",GENEMARKS_CALLS))
logfile.write("%s%s\n" % ("PRODIGAL_CALLS is ",PRODIGAL_CALLS))
logfile.write("%s%s\n" % ("GLIMMER_CALLS is ",GLIMMER_CALLS))
logfile.write("%s%s\n" % ("THEA_CALLS is ",THEA_CALLS))
logfile.write("%s%s\n" % ("DEBUG is ",DEBUG))

callCounts = {'prodigal' : 0, 'glimmer' : 0, 'genemarks' : 0, 'thea' : 0}

iterateGlimmer = False
#iterateGlimmer = True 
runCGC = False  # Turn on (below) if there are at least 2 gene callers being used 

########## CLASSES ##########
class geneCall:
    geneCallList = []

    def __init__(self, ID, method, contig, start, stop, strand, score):
        geneCall.geneCallList.append(self)
        self.ID = ID
        self.method = method
        self.contig = contig
        self.start = start
        self.stop = stop
        self.strand = strand
        self.score = score

    def printall(self):
        print self.ID, self.method, self.contig, self.start, self.stop, self.strand, self.score
        for gene in geneCallList:
            print gene

########## FUNCTIONS ##########

def systemCall(command):
    #print("\nSYSTEM CALL: "+command)
    print "\nSYSTEM CALL: ", command
    logfile.write("%s%s\n" % ("command is ",command))
    os.system(command)

def getProdigalId(x):
    colonSplit = x.split(";")
    equalSplit = colonSplit[0].split("=")
    return(equalSplit[1])

def getGeneMarkSId(x):
    colonSplit = x.split(",")
    equalSplit = colonSplit[0].split("=")
    return(equalSplit[1])

def processProdigal(line):
    if not line.startswith("#"):
        lineSplit = line.split("\t")

        contig = lineSplit[0]
        method = lineSplit[1]
        start = lineSplit[3]
        stop = lineSplit[4]
        score = lineSplit[5]
        strand = lineSplit[6]

        ID = getProdigalId(lineSplit[8])

	print "   ID=",ID,"method=",method,"contig=",contig,"start/stop=",start,'/',stop,"strand=",strand,"score=",score
        geneCall(ID, method, contig, start, stop, strand, score)
        callCounts['prodigal'] += 1

def processGlimmer(line,contig):
    contigSplit = contig.split(" ")
    contig = contigSplit[0]

    lineSplit = line.split()
    method = "glimmer3"
    ID = lineSplit[0]

    start = lineSplit[1]
    stop = lineSplit[2]
    strand = lineSplit[3][0]
    if strand == "-":
        start = lineSplit[2]
        stop = lineSplit[1]

    score = lineSplit[4]
    geneCall(ID, method, contig, start, stop, strand, score)

    callCounts['glimmer'] += 1

def processGeneMarkS(line):
    if not line.startswith("#"):
        lineSplit = line.split("\t")

        if len(lineSplit) > 1:

            contigSplit = lineSplit[0].split(" ")
            contig = contigSplit[0]

            method = lineSplit[1]
            start  = lineSplit[3]
            stop   = lineSplit[4]
            score  = lineSplit[5]
            strand = lineSplit[6]

            ID = getGeneMarkSId(lineSplit[8])

            geneCall(ID, method, contig, start, stop, strand, score)
            callCounts['genemarks'] += 1

def processThea(line):
    if not line.startswith('#'):

        split = line.split("\t")
        geneCall("NA", "THEA", split[3], split[0], split[1], split[2], "NA")
        #geneCall("NA", "THEA", "NA", split[0], split[1], split[2], "NA")
        callCounts['thea'] += 1

########## PRODIGAL ##########

# system call to prodigal
## the '-p' option will currently give an error that there is not enough sequence info to train.
## This will be addressed in version 3 (https://github.com/hyattpd/Prodigal/issues/11)

if PRODIGAL_CALLS:
    print("\n########## Prodigal ##########")
    logfile.write("%s\n" % ("Processing Prodigal"))

    command = prodigalPath + 'prodigal -q -i ' + fastaFileName + ' -o ' + outputFolder + 'prodigal.gff -f gff -p meta'
    systemCall(command)

    prodigalPeptideFile = outputFolder + "prodigal.proteins.faa"
    prodigalPotentialFile = outputFolder + "prodigal.genes.potential"

    command = prodigalPath + 'prodigal -i ' + fastaFileName + ' -o ' + outputFolder + 'prodigal.genes.sco -f sco -p meta -d ' + prodigalPeptideFile + ' -s ' + prodigalPotentialFile
    systemCall(command)
    files['Raw Prodigal GFF'] = outputFolder + 'prodigal.gff'

    prodigalFile = open(files['Raw Prodigal GFF'], 'r')
    lines = prodigalFile.read().splitlines()
    for line in lines:
        processProdigal(line)
    logfile.write("%s\n" % ("Prodigal processing complete."))

else:
    logfile.write("%s\n" % ("Not running Prodigal gene calling"))

########## GLIMMER ##########

## This runs the standard version of glimmer3.02b ("from scratch")

if DEBUG:
    logfile.write("%s\n" % ("Preparing to process Glimmer calls"))
    logfile.write("%s%s\n" % ("GLIMMER_CALLS is ",GLIMMER_CALLS))
if GLIMMER_CALLS:
    print("\n########## Glimmer ##########")
    logfile.write("%s\n" % ("Processing Glimmer"))
    logfile.write("%s%s\n" % ("glimmerPath is ",glimmerPath))
    logfile.write("%s%s\n" % ("fastaFileName is ",fastaFileName))
    logfile.write("%s%s\n" % ("outputFolder is ",outputFolder))

    systemCall(glimmerPath + 'long-orfs -n -t 1.15 ' + fastaFileName + ' ' + outputFolder + 'glimmer.longorfs')
    systemCall(glimmerPath + 'extract -t '           + fastaFileName + ' ' + outputFolder + 'glimmer.longorfs > ' + outputFolder + 'glimmer.train')
    systemCall(glimmerPath + 'build-icm -r ' + outputFolder + 'glimmer.icm < ' + outputFolder + 'glimmer.train')
    systemCall(glimmerPath + 'glimmer3 -o50 -g110 -t30 ' + fastaFileName + ' ' + outputFolder + 'glimmer.icm ' + outputFolder + 'glimmer')
    systemCall('tail -n +2 ' + outputFolder + 'glimmer.predict > ' + outputFolder + 'glimmer.coords') 

    glimmerOutputHandle = "glimmer"
    
    if iterateGlimmer == True:
        systemCall('tail -n +2 ' + outputFolder + 'glimmer.predict > ' + outputFolder + 'long-orfs.coords')    
        systemCall(glimmerPath + '/scripts/upstream-coords.awk 25 0 ' + outputFolder + 'long-orfs.coords | ' + glimmerPath + '/bin/extract ' + fastaFileName + ' - > ' + outputFolder + 'glimmer.upstream')
        systemCall('/data/data1/softwares/ELPH/bin/Linux-i386/elph ' + outputFolder + 'glimmer.upstream LEN=6 | ' + glimmerPath + '/scripts/get-motif-counts.awk > ' + outputFolder + 'glimmer.motif')
        systemCall('set startuse = \'' + glimmerPath + '/bin/start-codon-distrib -3 ' + fastaFileName + ' ' + outputFolder + '/long-orfs.coords\'')
        systemCall(glimmerPath + '/bin/glimmer3 -o50 -g110 -t30 -b ' + outputFolder + 'glimmer.motif -P $startuse ' + fastaFileName + ' ' + outputFolder + 'glimmer.icm ' + outputFolder + 'glimmerIterative')
        glimmerOutputHandle = "glimmerIterative"


    files['Raw Glimmer Output'] = outputFolder + '' + glimmerOutputHandle + '.predict'
    logfile.write("%s%s\n" % ("Raw Glimmer Output is ", files['Raw Glimmer Output']))

    currentGlimmerContig = ""
    for line in open(outputFolder + '' + glimmerOutputHandle + '.predict', 'rt'):
        line = line.rstrip()
        if line.startswith(">"):
            currentGlimmerContig = line[1:]
        else:
            processGlimmer(line,currentGlimmerContig)
    logfile.write("%s\n" % ("Processing Glimmer complete."))
else:
    logfile.write("%s\n" % ("Not running Glimmer gene calling"))

########## GENEMARKS ##########

if DEBUG:
    logfile.write("%s\n" % ("Preparing to process Genemarks calls"))
    logfile.write("%s%s\n" % ("GENEMARKS_CALLS is ",GENEMARKS_CALLS))
if GENEMARKS_CALLS:
    print("\n########## GeneMarkS ##########")
    logfile.write("%s\n" % ("Processing GeneMarkS"))

    systemCall(geneMarkSPath + 'gmhmmp -m ' + geneMarkSPath + '/heu_11.mod ' + fastaFileName + ' -o ' + outputFolder + 'geneMarkS.lst -r')
    systemCall(geneMarkSPath + 'gmhmmp -m ' + geneMarkSPath + '/heu_11.mod ' + fastaFileName + ' -o ' + outputFolder + 'geneMarkS.gff -r -f G')
    files['Raw GeneMarkS GFF'] = outputFolder + 'geneMarkS.gff'

    for line in open(outputFolder + 'geneMarkS.gff', 'rt'):
        line = line.rstrip()
        processGeneMarkS(line)
    logfile.write("%s\n" % ("Processing Genemarks complete."))
else:
    logfile.write("%s\n" % ("Not running GeneMarkS gene calling"))

########## THEA ##########

if DEBUG:
    logfile.write("%s\n" % ("Preparing to process THEA calls"))
    logfile.write("%s%s\n" % ("THEA_CALLS is ",THEA_CALLS))
if THEA_CALLS:
    print("\n########## THEA ##########")
    logfile.write("%s\n" % ("Processing THEA"))

    os.chdir(theaPath)
    systemCall('python thea.py ' + fastaFileName + ' > ' + outputFolder + 'theaOutput.txt' )
    os.chdir(workingFolder)

    for line in open(outputFolder + 'theaOutput.txt', 'rt'):
        line = line.rstrip()
        processThea(line)
    logfile.write("%s\n" % ("Processing THEA complete."))
else:
    logfile.write("%s\n" % ("Not running THEA gene calling"))

########## RESULTS ##########

print("\n########## RESULTS ##########")
logfile.write("%s\n" % ("Preparing results"))
results.write("ID\tSTART\tSTOP\tSTRAND\tSCORE\tMETHOD\tCONTIG\n")

if DEBUG:
    logfile.write("%s\n" % ("GENE_CALL_TESTING"))
    print "Printing all gene calls:"
    #geneCall.printall()

logfile.write("%s\n" % ("Writing genecall data to results file..."))
for gene in geneCall.geneCallList:
    results.write(gene.ID + "\t" + gene.start + "\t" + gene.stop + "\t" + gene.strand + "\t" + gene.score + "\t" + gene.method + "\t" + gene.contig + "\n")
results.close

logfile.write("%s\n" % ("Printing tallied genecall method call counts..."))
for method in callCounts:
    print(method + ": " + str(callCounts[method]))

print()

logfile.write("%s\n" % ("Printing file info..."))
for fileType in files:
    #print(fileType,files[fileType],sep=": ")
    print fileType, ": ", files[fileType],
print

logfile.write("%s\n" % ("Parsing genecall files into CGC format..."))
callerCount = 0
if GENEMARKS_CALLS:
    callerCount += 1
    systemCall('python ' + cgcPath + '/CGC_parser.py Genemarks ' + outputFolder + 'geneMarkS.gff ' + outputFolder + 'genemark.cgc')
if PRODIGAL_CALLS:
    callerCount += 1
    systemCall('python ' + cgcPath + '/CGC_parser.py Prodigal ' + outputFolder + 'prodigal.genes.sco ' + outputFolder + 'prodigal.cgc')
if GLIMMER_CALLS:
    callerCount += 1
    systemCall('python ' + cgcPath + '/CGC_parser.py Glimmer ' + outputFolder + 'glimmer.coords ' + outputFolder + 'glimmer.cgc')
if THEA_CALLS:
    callerCount += 1
    systemCall('python ' + cgcPath + '/CGC_parser.py THEA ' + outputFolder + 'theaOutput.txt ' + outputFolder + 'thea.cgc')

logfile.write("%s%s\n" % ("callerCount is ",callerCount))
if callerCount >= 2:
    runCGC = True

if runCGC:
    systemCall('python ' + cgcPath + '/CGC_main.py ' + outputFolder + '*.cgc > ' + outputFolder + 'CGC_results.txt')
else:
    logfile.write("%s\n" % ("Not running CGC code: too few gene callers to meaningfully compare"))


logfile.close()
