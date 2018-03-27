###################################################################################################
#
# Module:  CGC_compare.py
#
# Programmer:  Carol Zhou
#
# Description:  Module containing classes and methods for comparing results from different gene
#    callers.
#
# Updates:
#    Begin 3 June 2016
#
# Programmer's Notes:
#
# Classes and Methods:
#    Comparison
#        IdentifyCallers()
#        IdentifyCommonCore()
#        IsLesser(gene1,gene2)
#        Merge(nextGeneSet)
#        Compare()
#        PrintMergeList()
#        PrintUniqueList()
#        PrintCommonCore()
#        PrintCallerList()
#        PrintGenecallGrid()
#        PrintReport()
#        PrintStats()
#        PrintAll()
#        PrintAll_verbose()
#
###################################################################################################

# This code was developed by Carol L. Ecale Zhou at Lawrence Livermore National Laboratory.
# THIS CODE IS COVERED BY THE BSD LICENSE. SEE INCLUDED FILE BSD.pdf FOR DETAILS.

import re
import copy
import CGC_geneCall

p_comment   = re.compile('^#')

class Comparison(object):
    
    def __init__(self):
        self.commonCore = []  # list of lists of identical CGC_geneCall object calls (ie, different callers, same call) 
        self.mergeList  = []  # combined list of CGC_geneCall objects, merged by self.Merge()
        self.uniqueList = []  # list of lists of unique gene calls over all callers; each item in list is a common gene call (>=1 gene caller)
        self.callerList = []  # non-redundant list of callers 
        self.geneCall   = CGC_geneCall.GeneCall()  # a geneCall object

    # Create a non-redundant list of gene callers
    def IdentifyCallers(self):
        if self.mergeList:
            for gene in self.mergeList:
                if gene.geneCaller not in self.callerList:
                    self.callerList.append(gene.geneCaller)
            self.callerList.sort()
            return len(self.callerList)
        else:
            print "IdentifyCallers(): No callers to extract: call method Merge() to establish mergeList before calling this method" 
            return 0

    # Run Merge() and Compare() before running this method   
    def IdentifyCommonCore(self):
        # First, determine the callers used  
        if self.uniqueList:  # Must have previously called self.Compare() to fill this list 
            if self.mergeList:
                callerCount = self.IdentifyCallers() 
                if callerCount > 0:
                    for commonCalls in self.uniqueList:
                        count = 1
                        if len(commonCalls) == callerCount:
                            newCommonCoreCall = copy.deepcopy(self.geneCall)
                            geneName   = "CommonCoreGene_" + str(count)
                            strand     = commonCalls[0].strand
                            leftEnd    = commonCalls[0].leftEnd
                            rightEnd   = commonCalls[0].rightEnd
                            geneLength = commonCalls[0].geneLength
                            newCommonCoreCall.AssignGeneCall(geneName,"All_callers",count,strand,leftEnd,rightEnd,geneLength)
                            self.commonCore.append(newCommonCoreCall)
                            count += 1
                else:
                    print "IdentifyCommonCore(): callerCount is zero! cannot process"
            else:
                print "IdentifyCommonCore(): MergeList is empty:  need to run self.Merge()"
        else:
            print "IdentifyCommonCore(): No data available to identify common core"
        return 

    # Determine which gene call occurs first along the sequence
    def IsLesser(self,gene1,gene2):  # input is 2 geneCall objects
        if (int(gene1.leftEnd) < int(gene2.leftEnd)):
            return True
        if (int(gene1.leftEnd) == int(gene2.leftEnd)) and (int(gene1.rightEnd) < int(gene2.rightEnd)):
            return True
        return False

    # Merge a list of gene call objects with self.mergeList
    # Call this method once for each caller's output (i.e., loop over the set of gene caller outputs) 
    def Merge(self,nextGeneSet):  # Merge a list of gene call objects with self.mergeList
        if self.mergeList == []:
            # Trivial case: copy geneCall objects into self.mergeList
            for geneCall in nextGeneSet:
                newGeneCall = copy.deepcopy(geneCall)
                self.mergeList.append(newGeneCall)

        else: # merge next gene call set with existing self.mergeList
            temp = []   # holds merged geneCall objects
            mergeIndex = 0; mergeEnd = len(self.mergeList)
            nextIndex  = 0; nextEnd  = len(nextGeneSet)

            # Capture the gene call with smallest leftEnd location, or rightEnd if leftEnds are equal
            while mergeIndex < mergeEnd and nextIndex < nextEnd:
                if self.IsLesser(self.mergeList[mergeIndex],nextGeneSet[nextIndex]):
                    nextMergeGene = self.mergeList[mergeIndex]
                    mergeIndex += 1
                else:
                    nextMergeGene = copy.deepcopy(nextGeneSet[nextIndex])
                    nextIndex += 1
                temp.append(nextMergeGene) 

            # If this code executes, then all remaining geneCall objects from self.mergeList go next
            while mergeIndex < mergeEnd:
                temp.append(self.mergeList[mergeIndex])
                mergeIndex += 1

            # If this code executes, then all remaining geneCall objects from nextGeneSet go next
            while nextIndex < nextEnd:
                temp.append(nextGeneSet[nextIndex])
                nextIndex += 1

            # update self.mergeList
            self.mergeList = temp
        return

    # Compares the genes in self.mergeList; creates a list of ordered, unique gene calls 
    # Run this method after having merged all of your gene call sets into self.mergeList
    def Compare(self):  
        identityList = []
        if self.mergeList:
            callCount = len(self.mergeList)               # number of total gene calls, all callers
            nextCall  = copy.deepcopy(self.mergeList[0])  # capture 1st gene call
            identityList.append(nextCall)                 # each addition to list is identical to existing in list
            for i in xrange(1,callCount):                 # start with 2nd gene call
                nextCall = copy.deepcopy(self.mergeList[i])
                if self.mergeList[i].strand   != self.mergeList[i-1].strand  or \
                   self.mergeList[i].leftEnd  != self.mergeList[i-1].leftEnd or \
                   self.mergeList[i].rightEnd != self.mergeList[i-1].rightEnd:
                    self.uniqueList.append(identityList)  # all identicals for this gene call are identified
                    identityList = []                     # reset
                identityList.append(nextCall)
            if identityList:
                self.uniqueList.append(identityList)
        else:
            print "Compare(): Nothing to Compare"
        return

    def PrintMergeList(self):
        print "\n***************Merge List"
        count = 1 
        for gene in self.mergeList:
            print count,
            gene.PrintAll_brief()
            count += 1
        return
  
    def PrintUniqueList(self):
        print "\n***************Unique List"
        count = 1 
        for list in self.uniqueList:
            print "Unique gene call", count
            for gene in list:
                print "   ", 
                gene.PrintAll_brief()
            count += 1
        return

    def PrintCommonCore(self):
        print "\n***************Common Core"
        count = 1 
        for gene in self.commonCore:
            print count,
            gene.PrintAll_brief()
            count += 1
        return

    def PrintCallerList(self):
        print "\n***************List of Gene Callers"
        for caller in self.callerList:
            print caller 

    # Formats the unique calls list and prints to standard out; this is the final comparison data set
    def PrintGenecallGrid(self): # Prints an ordered, complete list of gene calls, each caller's in a column, identical calls in same row
        if self.callerList:
            if self.uniqueList: # Recall, uniqueList is list of unique gene calls, many of which were called by >1 caller
                count = 1

                # Print column headers, for as many gene callers as we have
                print "count\t",
                for i in xrange(0,len(self.callerList)):
                    print "caller\tstrand\tleftEnd\trightEnd\tlength\tcontig\t",
                print 

                # Format each gene call as a single line of output, arranging gene callers in order left to right
                for geneList in self.uniqueList:

                    # Create an empty array for printing identical gene calls in order going across by gene caller
                    printArray = [] # initialize
                    for i in xrange(0,len(self.callerList)): # Make the array as big as it needs to be for current gene call 
                        printArray.append('')

                    # Fill the print array in order by gene caller # will be at least 1 gene caller's call 
                    for i in xrange(0,len(geneList)):
                        currentCaller = geneList[i].geneCaller
                        printColumn = self.callerList.index(currentCaller) # capture index of this gene caller in self.callerList
                        printArray[printColumn] = geneList[i].geneCaller + '\t' + geneList[i].strand   + '\t' \
                                                + geneList[i].leftEnd    + '\t' + geneList[i].rightEnd + '\t' \
                                                + geneList[i].geneLength + '\t' + geneList[i].contig + '\t'

                    # Print the current row: horizontal list of identical gene calls 
                    print count, '\t',
                    for geneCallString in printArray:
                        if geneCallString == '':
                            print "\t\t\t\t\t\t",
                        else:
                            print geneCallString, 
                    print 
                    count += 1
            else:
                print "PrintGenecallGrid(): uniqueList is empty"
        else:
            print "PrintGenecallGrid(): callerList is empty"
        return

    def PrintReport(self):  # Final output
        self.PrintStats()
        self.PrintGenecallGrid()
        return

    def PrintStats(self):

        # Print a list of the callers 
        print "The following gene callers were considered:",
        for caller in self.callerList:
            print ',', caller,
        print
        print "The number of distinct gene calls over all gene callers is", len(self.uniqueList)
        print "The number of gene calls in common among all callers is", len(self.commonCore) 

        # Calculate number of gene calls that are not shared between any 2 gene callers
        loneCallCount = 0
        for callerList in self.uniqueList:
            if len(callerList) == 1:
                loneCallCount += 1
        print "The number of unique (non-matching) gene calls is", loneCallCount 

        # For each gene caller, calculate the number of calls it made 
        for caller in self.callerList:
            callCount        = 0
            cumulativeLength = 0 
            maxLength        = 0
            minLength        = 1000000 
            aveLength        = 0
            for call in self.mergeList:
                if (call.geneCaller == caller):
                    callCount += 1
                    intLength = int(call.geneLength)
                    cumulativeLength += intLength
                    if maxLength < intLength:
                        maxLength = intLength
                    if minLength > intLength:
                        minLength = intLength
            aveLength = cumulativeLength / callCount
            print "Caller", caller, "produced", callCount, "gene calls."
            print "Caller", caller, "gene-call length stats:  min:", minLength, ", max:", maxLength, ", ave:", aveLength

    def PrintAll(self):  # Print a dump of everything (debug/diagnostic) 
        self.PrintCallerList()
        self.PrintMergeList()
        self.PrintUniqueList()
        self.PrintCommonCore()
        self.PrintGenecallGrid()
        return

    def PrintAll_verbose(self):
        self.PrintCallerList()
        print "\n***************Merge List"
        for gene in self.mergeList:
            gene.PrintAll()
        print "\n***************Unique List"
        for list in self.uniqueList:
            for gene in list:
                gene.PrintAll()
        print "\n***************Common Core"
        for gene in self.commonCore:
            gene.PrintAll()
        print "\n***************GeneCall Grid"
        self.PrintGenecallGrid()
        return


