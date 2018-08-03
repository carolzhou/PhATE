#!/usr/bin/env python

################################################################
#
# pullSeq.py
#
# Description: Pulls specified sequence from a fasta file. 
#
# Usage:  python pullSeq.py fastaFile start end 
#
# Programmers: 
#    Carol E. Zhou - 
#
################################################################

# This code was developed by Carol L. Ecale Zhou at Lawrence Livermore National Laboratory.

import sys, os, re, string, copy, time, datetime
from subprocess import call

# Get input parameters - expect exactly 3 arguments: name of genome file, start coordinate, end coordinate

genomeFile = ""
startCoord = 0
endCoord = 0
outCharList = []
outString = ""

if len(sys.argv) != 4:
    print "Input the name of your genome fasta file followed by start and end coordinates" 
else:
    genomeFile = sys.argv[1]
    startCoord = sys.argv[2]
    endCoord   = sys.argv[3]

# Open genome fasta file and output designated bases

GENOME = open(genomeFile,"r")
gLines = GENOME.read().splitlines()

cCount = 0
for gLine in gLines:
    charList = list(gLine)
    for char in charList:
        cCount += 1
        if cCount >= int(startCoord) and cCount <= int(endCoord):
            outCharList.append(char)
outString = ''.join(outCharList)
print outString
GENOME.close()
