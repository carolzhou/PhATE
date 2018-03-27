################################################################################################
#
# Module:  CGC_geneCall.py
#
# Programmer:  Carol Zhou
#
# Description:  Module containing classes and methods for handling data output from a gene caller 
#
# Updates:
#    Begin 2 June 2016
#    25 Jan 2018:  adding protein names
#
# Programmer's Notes:
#
# Classes and Methods:
#    GeneCall()
#        AssignGeneCall(<input parameters>)
#        PrintAll()
#        PrintAll_brief()
#    GeneCallSet()
#        AddGeneCall(newGeneCall)
#        AddGeneCalls(GENE_FILE_HANDLE)
#        IsLesser(gene1,gene2)
#        UpdateGeneCount()
#        GetGeneCalls()
#        SortGeneCalls()
#        PrintAll()
#        PrintAll_brief()
#
#################################################################################################

# This code was developed by Carol L. Ecale Zhou at Lawrence Livermore National Laboratory.
# THIS CODE IS COVERED BY THE BSD LICENSE. SEE INCLUDED FILE BSD.pdf FOR DETAILS.

import re
import copy

p_comment    = re.compile('^#')
p_caller     = re.compile('([\w\d]+)\sgene\scalls')
p_callerName = re.compile('[Gg][Ee][Nn][Ee][Mm][Aa][Rr][Kk]|[Gg][Ll][Ii][Mm][Mm][Ee][Rr]|[Pp][Rr][Oo][Dd][Ii][Gg][Aa][Ll]|[Rr][Aa][Ss][Tt]|[Tt][Hh][Ee][Aa]|[Gg][Ff][Ff][3]')
#p_callerName = re.compile('[Gg][Ee][Nn][Ee][Mm][Aa][Rr][Kk]|[Gg][Ll][Ii][Mm][Mm][Ee][Rr]|[Pp][Rr][Oo][Dd][Ii][Gg][Aa][Ll]|[Rr][Aa][Ss][Tt]|[Pp][Hh][Aa][Tt][Ee]')
p_dataLine   = re.compile('^(\d+)\t([+-])\t(\d+)\t(\d+)\t(\d+)\t([\d\w\.\-\_]+)')

class GeneCall(object):
    
    def __init__(self):
        self.geneName    = "unknwon"
        self.geneCaller  = "unknown"
        self.geneNumber  = 0
        self.strand      = 'x'  # Typically '+' or '-'; could be '?'; 'x' indicates NULL
        self.leftEnd     = 0    # may be start or stop, depending on strand (orientation)
        self.rightEnd    = 0
        self.geneLength  = 0
        self.contig      = "unknown"
        self.protein     = "unknown"

    def AssignGeneCall(self,geneName,geneCaller,geneNumber,strand,leftEnd,rightEnd,geneLength,contig="unknown",protein="unknown"):
        self.geneName    = geneName
        self.geneCaller  = geneCaller
        self.geneNumber  = geneNumber
        self.strand      = strand 
        self.leftEnd     = leftEnd 
        self.rightEnd    = rightEnd 
        self.geneLength  = geneLength
        self.contig      = contig
        self.protein     = protein
        return

    def PrintAll(self):
        print "\ngeneName =", self.geneName
        print "geneCaller =", self.geneCaller
        print "geneNumber =", self.geneNumber
        print "leftEnd =",    self.leftEnd
        print "rightEnd =",   self.rightEnd
        print "strand =",     self.strand
        print "length =",     self.geneLength
        print "contig=",      self.contig
        print "protein=",     self.proteinName
        return

    def PrintAll_brief(self):
        print "Gene No.", self.geneNumber, "gene caller: ", self.geneCaller, ", leftEnd:", self.leftEnd, ", rightEnd:", self.rightEnd, ", strand:", self.strand, ", length:", self.geneLength, ", contig:", self.contig, ", protein:", self.protein
        return

class GeneCallSet(object):

    def __init__(self):
        self.geneCaller     = ""  # Typically, 'GeneMark', 'Glimmer', 'Prodigal', 'RAST', 'THEA'
        self.geneCount      = 0
        self.geneCallList   = []  # list of GeneCall objects 
        self.geneCall_obj   = GeneCall()

    def UpdateGeneCount(self):

        self.geneCount = len(self.geneCallList)
        return self.geneCount

    def GetGeneCalls(self,fLines,geneCaller):

        for line in fLines:
            match_data = re.search(p_dataLine,line)
            if match_data:
                geneNumber = match_data.group(1)
                strand     = match_data.group(2)
                leftEnd    = match_data.group(3)
                rightEnd   = match_data.group(4)
                geneLength = match_data.group(5) 
                contig     = match_data.group(6)
                geneName   = self.geneCaller + '_' + geneNumber
                newGeneCall = copy.deepcopy(self.geneCall_obj)
                newGeneCall.AssignGeneCall(geneName,geneCaller,geneNumber,strand,leftEnd,rightEnd,geneLength,contig)
                self.AddGeneCall(newGeneCall)
        return

    def AddGeneCall(self,newGeneCall):

        self.geneCallList.append(newGeneCall)
        self.UpdateGeneCount()
        return

    def AddGeneCalls(self,GENE_FILE_HANDLE):

        fLines = GENE_FILE_HANDLE.read().splitlines()
        for line in fLines:
            match_caller = re.search(p_caller,line)
            if match_caller:
                caller = match_caller.group(1).lower()
                match_callerName = re.search(p_callerName,caller)
                if match_callerName:
                    self.geneCaller = caller 
                    self.GetGeneCalls(fLines,caller)
                else:
                    print "ERROR: gene caller not recognized in geneCall.GeneCallSet,", caller, line
        return

    # Determine which of 2 gene calls occurs first along the sequence (left to right, regardless of orientation) 
    def IsLesser(self,gene1,gene2):  # Input is 2 geneCall objects

        # Sort on left end position
        if (int(gene1.leftEnd) < int(gene2.leftEnd)):
            return True

        # Sort on end position if start positions are equal
        if (int(gene1.leftEnd) == int(gene2.leftEnd)) and (int(gene1.rightEnd) < int(gene2.rightEnd)):
            return True 

        # If you're still here, then gene1 > gene2
        return False 
 
    def SortGeneCalls(self):

        # Use insertion sort because it is fast when the list is nearly sorted to begin with 
        # Sort from lowest start position to highest; sort on end position if starts are equal
        for index in xrange(0,len(self.geneCallList)-1):
            currentValue = self.geneCallList[index]
            position = index
            while (position > 0) and self.IsLesser(currentValue,self.geneCallList[position-1]):
                self.geneCallList[position] = self.geneCallList[position-1]
                position = position - 1
            self.geneCallList[position] = currentValue
        return

    def PrintAll(self):

        print "Gene Caller: ",self.geneCaller
        for gene in self.geneCallList:
            gene.PrintAll()
        return

    def PrintAll_brief(self):

        print "Gene Caller: ",self.geneCaller
        for gene in self.geneCallList:
            gene.PrintAll_brief()
        return


