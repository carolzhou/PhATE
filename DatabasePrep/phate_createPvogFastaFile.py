#!/usr/bin/env python

#####################################################################
#
# phate_createPvogFastaFile.py
#
# Description: This code reads input files, AllFamilyProteinList.tab and NCBI_Phage.faa,
#    and creates a file containing pVOG-tagged protein fastas.
#
# Usage:  python createPvogFsataFile.py
#
# Programmers Notes:
#
# Updates:
#    23 June 2017: begin
#
# Programmer:  CEZhou
#
####################################################################

# This code was developed by Carol L. Ecale Zhou at Lawrence Livermore National Laboratory.
# THIS CODE IS COVERED BY THE BSD LICENSE. SEE INCLUDED FILE BSD.PDF FOR DETAILS.


import sys
import os
import re
import string
import copy
import time
import datetime
from subprocess import call 

# Control
CHATTY = True

CODE_BASE = "createPvogFastaFile"

# Paths to files; when running on mpath

PVOG_BASE_DIR         = "/home/zhou4/PhATE/Databases/pVOGs/"

# Files

PROTEIN_LIST          = PVOG_BASE_DIR + "AllFamilyProteinList.tab"
NCBI_PHAGE_FAA        = PVOG_BASE_DIR + "NCBI_phage.faa"
PVOG_TAGGED_PROTEINS  = PVOG_BASE_DIR + "pVOGs.faa"
PVOG_MISSING_PROTEINS = PVOG_BASE_DIR + "pVOGs_missing.lst"
LOG                   = PVOG_BASE_DIR + CODE_BASE + ".log"
PVOG_LIB_DIR          = PVOG_BASE_DIR + "Allvogtables/"

# Dat# Environment variables, which are global to any instance of this code's execution
MPATH = True
SANDBOX = False

if MPATH:
    BASE_DIR = "/home/zhou4/PhATE/Code/PhATE_pipeline/"
    os.environ["BLAST_HOME"]               = "/usr/bin/"
    os.environ["KEGG_VIRUS_BASE_DIR"]      = "/home/zhou4/BioRemediation/Phage/Kegg_virus/"
    os.environ["KEGG_VIRUS_BLAST_HOME"]    = os.environ["KEGG_VIRUS_BASE_DIR"] + "BlastDBs/T40000.pep"
    os.environ["NCBI_VIRUS_BASE_DIR"]      = "/home/zhou4/BioRemediation/Phage/NCBI_virus_genomes/"
    os.environ["NCBI_VIRUS_BLAST_HOME"]    = os.environ["NCBI_VIRUS_BASE_DIR"] + "viral.1.1.genomic.fna"
    os.environ["PHANTOME_BASE_DIR"]        = "/home/zhou4/BioRemediation/Phage/Phantome/"
    os.environ["PHANTOME_BLAST_HOME"]      = os.environ["PHANTOME_BASE_DIR"] + "Phantome_Phage_genes.faa"
    os.environ["UNIPARC_BASE_DIR"]         = "/data/data1/data/UniParc/"
    os.environ["UNIPARC_VIRUS_BLAST_HOME"] = os.environ["UNIPARC_BASE_DIR"] + "insertSubDir/insertName"  #***uniparc_active.fasta ???
    os.environ["NR_BLAST_HOME"]            = "/data/data1/sandbox/BLAST/nr"
    os.environ["EMBOSS_HOME"]              = "/data/data1/softwares/EMBOSS-6.6.0/emboss/"
    os.environ["NCBI_TAXON_DIR"]           = os.environ["NCBI_VIRUS_BASE_DIR"]
    os.environ["PIPELINE_DIR"]             = BASE_DIR
    os.environ["PSAT_OUT_DIR"]             = BASE_DIR
    os.environ["PRODIGAL_PATH"]            = "/data/data1/softwares/prodigal.v2_50/"
    os.environ["GLIMMER_PATH"]             = "/data/data1/softwares/glimmer3.02/"
    os.environ["GENEMARKS_PATH"]           = "/data/data1/softwares/GeneMarkS/genemark_suite_linux_64/gmsuite/"
    os.environ["PHATE_PATH"]               = "/data/data1/softwares/PHATE/PHATE-0.5b/"
    os.environ["CGC_PATH"]                 = "/home/zhou4/PhATE/Code/PhATE_pipeline/CompareCalls/"
elif SANDBOX:
    BASE_DIR = "/data/data1/sandbox/PhATE_pipeline/"
    os.environ["BLAST_HOME"]               = "/usr/bin/"
    os.environ["KEGG_VIRUS_BASE_DIR"]      = "/data/data1/sandbox/Phage_pipeline/BlastDBs/Kegg_virus/"
    os.environ["KEGG_VIRUS_BLAST_HOME"]    = os.environ["KEGG_VIRUS_BASE_DIR"] + "BlastDBs/T40000.pep"
    os.environ["NCBI_VIRUS_BASE_DIR"]      = "/data/data1/sandbox/Phage_pipeline/BlastDBs/NCBI_virus_genomes/"
    os.environ["NCBI_VIRUS_BLAST_HOME"]    = os.environ["NCBI_VIRUS_BASE_DIR"] + "viral.1.1.genomic.fna"
    os.environ["PHANTOME_BASE_DIR"]        = "/data/data1/sandbox/Phage_pipeline/BlastDBs/Phantome/"
    os.environ["PHANTOME_BLAST_HOME"]      = os.environ["PHANTOME_BASE_DIR"] + "Phantome_Phage_genes.faa"
    os.environ["UNIPARC_BASE_DIR"]         = "/data/data1/data/UniParc/"
    os.environ["NR_BLAST_HOME"]            = "/data/data1/sandbox/BLAST/"
    os.environ["EMBOSS_HOME"]              = "/data/data1/softwares/Emboss/EMBOSS-6.6.0/emboss/"
    os.environ["NCBI_TAXON_DIR"]           = os.environ["NCBI_VIRUS_BASE_DIR"]
    os.environ["PIPELINE_DIR"]             = BASE_DIR
    os.environ["PSAT_OUT_DIR"]             = BASE_DIR
    os.environ["PRODIGAL_PATH"]            = "/data/data1/softwares/prodigal.v2_50/"
    os.environ["GLIMMER_PATH"]             = "/data/data1/softwares/glimmer3.02/"
    os.environ["GENEMARKS_PATH"]           = "/data/data1/softwares/GeneMarkS/genemark_suite_linux_64/gmsuite/"
    os.environ["PHATE_PATH"]               = "/data/data1/softwares/PHATE/PHATE-0.5b/"
    os.environ["CGC_PATH"]                 = "/home/zhou4/PhATE/Code/PhATE_pipeline/CompareCalls/"

# Import PhATE Modules

import phate_annotation
import phate_fastaSequence
import phate_pVOG

# Open Files

PROTEIN_LIST_H          = open(PROTEIN_LIST, 'r')
NCBI_PHAGE_FAA_H        = open(NCBI_PHAGE_FAA, 'r')
PVOG_TAGGED_PROTEINS_H  = open(PVOG_TAGGED_PROTEINS, 'w')
PVOG_MISSING_PROTEINS_H = open(PVOG_MISSING_PROTEINS, 'w') 
LOG_H                   = open(LOG, 'w')
LOG_H.write("%s%s\n" % ("Processing began at ",datetime.datetime.now()))

# Capture fastas from ncbi phage fasta file

if CHATTY:
    print "Reading ncbi phage fasta file"
ncbiPhageSeqs = phate_fastaSequence.multiFasta()
fastaLines = NCBI_PHAGE_FAA_H.read().splitlines()
ncbiPhageSeqs.addFastas(fastaLines,'aa')
if CHATTY:
    print "done!"

# Read in the pVOG identifiers and their associated accession numbers

if CHATTY:
    print "Reading pVOGs from pVOG file"
pVOGdb = phate_pVOG.pVOGs()
pVOGlines = PROTEIN_LIST_H.read().splitlines()
#pVOGdb.addPvogs_old(pVOGlines,LOG_H)
pVOGdb.addPvogs(PVOG_LIB_DIR,LOG_H)
if CHATTY:
    print "pVOGs have been recorded"
LOG_H.write("%s%s%s\n" % ("There are ",len(pVOGdb.pVOGlist)," pVOGs"))
accessionCount = pVOGdb.getAccessionCount()
LOG_H.write("%s%s\n" % ("The total number of accessions is ",accessionCount))

# Visit each pVOG and find the fasta sequence that corresponds to each member accession
# Modify the header of each identified fasta to reflect its membership in the pVOG cluster

if CHATTY:
    print "Searching sequences for each pVOG-associated accession"

# Create a fasta object (to be replicated as needed)
nextFasta = phate_fastaSequence.fasta()

# For each pVOG, get its associated peptide accessions, then
# Find that fasta in the ncbi database subset, and 
# Tag the fasta header with the pVOG information
for pVOG in pVOGdb.pVOGlist:
    foundCount = 0; missingCount = 0; missingList = []
    for accession in pVOG.accessionList:                    #*** Here need to modify to produce more complete tags for headers
        if CHATTY:
            print "Processing pVOG", pVOG.pVOGid, "and accession", accession
        LOG_H.write("%s%s%s%s\n" % ("Processing pVOG ",pVOG.pVOGid," and accession ",accession))
        nextFasta = ncbiPhageSeqs.findStringInHeader(accession)
        if nextFasta:
            nextFasta.pVOGassociationList.append(pVOG.pVOGid)  
            pVOGstring = ""
            for pVOGidentifier in nextFasta.pVOGassociationList:
                pVOGstring += pVOGidentifier + '|'
            nextFasta.customHeader = pVOGstring + nextFasta.header 
            if len(nextFasta.pVOGassociationList) > 1:
                LOG_H.write("%s%s\n" % ("WARNING: Fasta associated with multiple pVOGs: ", len(nextFasta.pVOGassociationList)))
                LOG_H.write("%s%s\n" % ("   header is: ", nextFasta.customHeader))
            foundCount += 1
        else:
            if accession not in missingList:
                missingList.append(accession)
                missingCount += 1
                PVOG_MISSING_PROTEINS_H.write("%s%s\n" % ("No pVOG found for accession ", accession))
if CHATTY:
    print "pVOG accession sequences have been tagged"

# Write the pVOG-tagged fastas to file
if CHATTY:
    print "Writing pVOG-tagged fasta sequences to file"
ncbiPhageSeqs.printMultiFasta2file_custom(PVOG_TAGGED_PROTEINS_H)
if CHATTY:
    print "done!"

LOG_H.write("%s%s\n" % ("Number of fasta sequences found = ", foundCount))
LOG_H.write("%s%s\n" % ("Number missing = ", missingCount))

# Clean up

PROTEIN_LIST_H.close()
NCBI_PHAGE_FAA_H.close()
PVOG_TAGGED_PROTEINS_H.close()
PVOG_MISSING_PROTEINS_H.close()
print "Processing complete!"
LOG_H.write("%s%s\n" % ("Processing completed at ",datetime.datetime.now()))
LOG_H.close()
