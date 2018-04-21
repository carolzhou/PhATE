#!/usr/bin/env python

################################################################
#
# Next:
# *) Consolidate log files.
# *) Implement user-configurable blast parameters.
# *) Incorporate blast against refseqgene, swissprot
#
# phate_runPipeline.py
#
# Description: Runs the phate annotation pipeline.  This code runs under Python 2.7, and requires
#    dependent packages.
#
# Usage:  python phate_runPipeline.py phate.config
#    (see phate_sample.config)
#
# Setup:
#     CompareCalls/         - code for comparing gene caller results
#     DatabasePrep/         - code for preparing custom databases
#     GeneCalling/          - mini-pipeline runs gene-call programs
#     SequenceAnnotation/   - PhATE sequence annotation codes
#     phate_runPipeline.py  - pipeline driver
#     PipelineInput/        - contains myGenome.fasta (and optionally, myGenome.psat)
#     PipelineOutput/       - output files are written here to a subdirectory specified in config file
#     phate.config          - configuration file 
#
# Programmer's Notes:
#    This code uses a running log; need to occasionally clean it out
#
# Programmers: 
#    Carol E. Zhou - pipeline programmer: CompareCalls/, DatabasePrep/, SequenceAnnotation, phate_runPipeline.py
#    Jeff Kimbrel  - GeneCalling/
#
################################################################

# This code was developed by Carol L. Ecale Zhou at Lawrence Livermore National Laboratory.
# THIS CODE IS COVERED BY THE BSD LICENSE. SEE INCLUDED FILE BSD.PDF FOR DETAILS.

import sys, os, re, string, copy, time, datetime
from subprocess import call


# Constants and Configurables

#various parameters
CONSENSUS_CALLS_FILE     = 'thea.cgc' #*** For now this is THEA calls, though may be consensus calls in future
GENE_FILE                = 'gene.fnt'                              #
PROTEIN_FILE             = 'protein.faa'                           # 
GENETIC_CODE             = '11'      # default is bacterial (11)
GENE_CALLER              = 'thea'    # default is annotation of phage, so THEA is preferred gene caller; if bac, could be 'consensus', 'genemark', 'glimmer', or 'prodigal'
GENOME_TYPE              = 'phage'   # default is phage; could be 'bacterium'
NAME                     = 'unknown' # user provided
CONTIG_NAME              = 'unknown' # user provided: temporary, finished genomes/single contig only for now
SPECIES                  = 'unknown' # user provided

#blast parameters
MAX_BLAST_HIT_COUNT      = 100       # maximum number of hits to capture (user should specify far fewer than max)
MIN_BLASTP_IDENTITY      = '20'      # default; sets a lower limit based on value at which a structure model can provide information
MAX_BLASTP_HIT_COUNT     = '100'     # default; sets an upper limit; user's value should typically be well below this 
MAX_BLASTN_HIT_COUNT     = '10'      # default; sets an upper limit
BLASTP_IDENTITY_DEFAULT  = '60'
BLASTP_HIT_COUNT_DEFAULT = '3'
BLASTN_HIT_COUNT_DEFAULT = '3'

#blast databases
NCBI_VIRUS_BLAST_DEFAULT         = True 
NCBI_VIRUS_PROTEIN_BLAST_DEFAULT = True
KEGG_VIRUS_BLAST_DEFAULT         = False     # Requires license
NR_BLAST_DEFAULT                 = False     # Large data set; blast run takes time
REFSEQ_PROTEIN_BLAST_DEFAULT     = True      # Large data set; blast run takes time
PHANTOME_BLAST_DEFAULT           = True
PVOGS_BLAST_DEFAULT              = True
UNIPARC_BLAST_DEFAULT            = False     # Turned 'off' for now; not yet in service
REFSEQ_GENE_BLAST_DEFAULT        = True 
SWISSPROT_BLAST_DEFAULT          = True

#gene callers
GENEMARKS_CALLS_DEFAULT          = False     # Requires license
PRODIGAL_CALLS_DEFAULT           = True
GLIMMER_CALLS_DEFAULT            = True
THEA_CALLS_DEFAULT               = True

#other
PSAT_ANNOTATION_DEFAULT          = False     # Requires LLNL processing

# Which computer/workspace am I running this on?
MYCLUSTER = False    # Your working directory
HOME_LAPTOP = True

# Environment variables, which are global to any instance of this code's execution
# You need to set these according to where the data sets and codes reside
if MYCLUSTER:  # machine running code 
    BASE_DIR = "/myBaseDir/myDir/PhATE/Code/PhATE_pipeline/" 
    os.environ["PIPELINE_DIR"]                  = BASE_DIR
    os.environ["PSAT_OUT_DIR"]                  = BASE_DIR # only on LLNL system

    # Data sets
    os.environ["KEGG_VIRUS_BASE_DIR"]           = "/Users/yourname/DEV/PhATE/Databases/KEGG/Kegg_virus_22Aug2017/"
    os.environ["KEGG_VIRUS_BLAST_HOME"]         = os.environ["KEGG_VIRUS_BASE_DIR"] + "T40000.pep"
    os.environ["NCBI_VIRUS_BASE_DIR"]           = "/Users/yourname/DEV/PhATE/Databases/NCBI/"
    os.environ["NCBI_VIRUS_BLAST_HOME"]         = os.environ["NCBI_VIRUS_BASE_DIR"] + "Virus_Genome/"  + "viral.1.1.genomic.fna"
    os.environ["NCBI_VIRUS_PROTEIN_BLAST_HOME"] = os.environ["NCBI_VIRUS_BASE_DIR"] + "Virus_Protein/" + "viral.1.protein.faa"
    os.environ["NCBI_TAXON_DIR"]                = os.environ["NCBI_VIRUS_BASE_DIR"] + "Virus_Genome/" 
    os.environ["PHANTOME_BASE_DIR"]             = "/Users/yourname/DEV/PhATE/Databases/Phantome/"
    os.environ["PHANTOME_BLAST_HOME"]           = os.environ["PHANTOME_BASE_DIR"] + "Phantome_Phage_genes.faa"
    os.environ["PVOGS_BASE_DIR"]                = "/Users/yourname/DEV/PhATE/Databases/pVOGs/"
    os.environ["PVOGS_BLAST_HOME"]              = os.environ["PVOGS_BASE_DIR"] + "pVOGs.faa"
    os.environ["UNIPARC_BASE_DIR"]              = "/Users/yourname/DEV/PhATE/Databases/UniParc/" # Uniparc not yet in service
    os.environ["UNIPARC_VIRUS_BLAST_HOME"]      = os.environ["UNIPARC_BASE_DIR"] + "insertSubDir/insertName"  #***uniparc_active.fasta ???
    os.environ["NR_BLAST_BASE_DIR"]             = "/Users/yourname/DEV/PhATE/Databases/NR/"
    os.environ["NR_BLAST_HOME"]                 = os.environ["NR_BLAST_BASE_DIR"] + "nr"
    os.environ["REFSEQ_PROTEIN_BASE_DIR"]       = "/Users/yourname/DEV/PhATE/Databases/Refseq/"
    os.environ["REFSEQ_PROTEIN_BLAST_HOME"]     = os.environ["REFSEQ_PROTEIN_BASE_DIR"] + "refseq_protein"

    # Gene calling
    os.environ["PRODIGAL_PATH"]                 = "/usr/local/bin/"
    os.environ["GLIMMER_PATH"]                  = "/Users/yourname/DEV/PhATE/OtherCodes/Glimmer/glimmer3.02/bin/"
    os.environ["GENEMARKS_PATH"]                = "/Users/yourname/DEV/PhATE/OtherCodes/GeneMarkS/genemark_suite_linux_64/gmsuite/"  
    os.environ["THEA_PATH"]                     = "/Users/yourname/DEV/PhATE/OtherCodes/THEA/THEA-master/"
    os.environ["CGC_PATH"]                      = "/Users/yourname/DEV/PhATE/Code/CompareCalls/"

    # Blast
    os.environ["BLAST_HOME"]                    = "/Users/yourname/DEV/PhATE/OtherCodes/Blast/ncbi-blast-2.7.1+/bin/"
    os.environ["MIN_BLASTP_IDENTITY"]           = MIN_BLASTP_IDENTITY
    os.environ["MAX_BLASTP_HIT_COUNT"]          = MAX_BLASTP_HIT_COUNT 
    os.environ["MAX_BLASTN_HIT_COUNT"]          = MAX_BLASTN_HIT_COUNT 
    os.environ["BLASTP_IDENTITY_DEFAULT"]       = BLASTP_IDENTITY_DEFAULT
    os.environ["BLASTP_HIT_COUNT_DEFAULT"]      = BLASTP_HIT_COUNT_DEFAULT
    os.environ["BLASTN_HIT_COUNT_DEFAULT"]      = BLASTN_HIT_COUNT_DEFAULT

elif HOME_LAPTOP:  # machine running code 
    BASE_DIR = "/Users/carolzhou/DEV/PhATE/Code/" 
    os.environ["EMBOSS_HOME"]                   = "/Users/carolzhou/DEV/PhATE/OtherCodes/EMBOSS/EMBOSS-6.6.0/emboss/" #*** name of code is transeq.c
    os.environ["PIPELINE_DIR"]                  = BASE_DIR
    os.environ["PSAT_OUT_DIR"]                  = BASE_DIR # only on LLNL system

    # Data sets
    os.environ["KEGG_VIRUS_BASE_DIR"]           = "/Users/carolzhou/DEV/PhATE/Databases/KEGG/Kegg_virus_22Aug2017/"
    os.environ["KEGG_VIRUS_BLAST_HOME"]         = os.environ["KEGG_VIRUS_BASE_DIR"] + "T40000.pep"
    os.environ["NCBI_VIRUS_BASE_DIR"]           = "/Users/carolzhou/DEV/PhATE/Databases/NCBI/"
    os.environ["NCBI_VIRUS_BLAST_HOME"]         = os.environ["NCBI_VIRUS_BASE_DIR"] + "Virus_Genome/"  + "viral.1.1.genomic.fna"
    os.environ["NCBI_VIRUS_PROTEIN_BLAST_HOME"] = os.environ["NCBI_VIRUS_BASE_DIR"] + "Virus_Protein/" + "viral.1.protein.faa"
    os.environ["NCBI_TAXON_DIR"]                = os.environ["NCBI_VIRUS_BASE_DIR"] + "Virus_Genome/" 
    os.environ["PHANTOME_BASE_DIR"]             = "/Users/carolzhou/DEV/PhATE/Databases/Phantome/"
    os.environ["PHANTOME_BLAST_HOME"]           = os.environ["PHANTOME_BASE_DIR"] + "Phantome_Phage_genes.faa"
    os.environ["PVOGS_BASE_DIR"]                = "/Users/carolzhou/DEV/PhATE/Databases/pVOGs/"
    os.environ["PVOGS_BLAST_HOME"]              = os.environ["PVOGS_BASE_DIR"] + "pVOGs.faa"
    os.environ["UNIPARC_BASE_DIR"]              = "/Users/carolzhou/DEV/PhATE/Databases/UniParc/" # Uniparc not yet in service
    os.environ["UNIPARC_VIRUS_BLAST_HOME"]      = os.environ["UNIPARC_BASE_DIR"] + "insertSubDir/insertName"  #***uniparc_active.fasta ???
    os.environ["NR_BLAST_BASE_DIR"]             = "/Users/carolzhou/DEV/PhATE/Databases/NR/"
    os.environ["NR_BLAST_HOME"]                 = os.environ["NR_BLAST_BASE_DIR"] + "nr"
    os.environ["REFSEQ_PROTEIN_BASE_DIR"]       = "/Users/carolzhou/DEV/PhATE/Databases/Refseq/Protein/"
    os.environ["REFSEQ_PROTEIN_BLAST_HOME"]     = os.environ["REFSEQ_PROTEIN_BASE_DIR"] + "refseq_protein"
    os.environ["REFSEQ_GENE_BASE_DIR"]          = "/Users/carolzhou/DEV/PhATE/Databases/Refseq/Gene/"
    os.environ["REFSEQ_GENE_BLAST_HOME"]        = os.environ["REFSEQ_GENE_BASE_DIR"] + "refseqgene"
    os.environ["SWISSPROT_BASE_DIR"]            = "/Users/carolzhou/DEV/PhATE/Databases/Swissprot/"
    os.environ["SWISSPROT_BLAST_HOME"]          = os.environ["SWISSPROT_BASE_DIR"] + "swissprot"

    # Gene calling
    os.environ["PRODIGAL_PATH"]                 = "/usr/local/bin/"
    os.environ["GLIMMER_PATH"]                  = "/Users/carolzhou/DEV/PhATE/OtherCodes/Glimmer/glimmer3.02/bin/"
    os.environ["GENEMARKS_PATH"]                = "/Users/carolzhou/DEV/PhATE/OtherCodes/GeneMarkS/genemark_suite_linux_64/gmsuite/"  
    os.environ["THEA_PATH"]                     = "/Users/carolzhou/DEV/PhATE/OtherCodes/THEA/THEA-master/"
    os.environ["CGC_PATH"]                      = "/Users/carolzhou/DEV/PhATE/Code/CompareCalls/"

    # Blast
    os.environ["BLAST_HOME"]                    = "/Users/carolzhou/DEV/PhATE/OtherCodes/Blast/ncbi-blast-2.7.1+/bin/"
    os.environ["MIN_BLASTP_IDENTITY"]           = MIN_BLASTP_IDENTITY
    os.environ["MAX_BLASTP_HIT_COUNT"]          = MAX_BLASTP_HIT_COUNT 
    os.environ["MAX_BLASTN_HIT_COUNT"]          = MAX_BLASTN_HIT_COUNT 
    os.environ["BLASTP_IDENTITY_DEFAULT"]       = BLASTP_IDENTITY_DEFAULT
    os.environ["BLASTP_HIT_COUNT_DEFAULT"]      = BLASTP_HIT_COUNT_DEFAULT
    os.environ["BLASTN_HIT_COUNT_DEFAULT"]      = BLASTN_HIT_COUNT_DEFAULT

else:
    print """You need to set the environment variables in """ + CODE + """\n"""

# Set environment variables

CODE_BASE   = "phate_runPipeline"
CODE        = CODE_BASE + ".py"
CONFIG_FILE = "phate.config"  # by default, but user should name their own, ending in ".config"

# Subordinate codes
GENECALL_CODE_DIR       = BASE_DIR + "GeneCalling/"         # Here resides GENECALL_CODE plus CGC (Compare Gene Calls) codes
SEQANNOT_CODE_DIR       = BASE_DIR + "SequenceAnnotation/"  # Performs sequence annotation via blast and incorporates PSAT output
GENECALL_CODE           = GENECALL_CODE_DIR + "phate_genecallPhage.py"
SEQANNOTATION_CODE_BASE = "phate_sequenceAnnotation_main"
SEQANNOTATION_CODE      = SEQANNOT_CODE_DIR + SEQANNOTATION_CODE_BASE + ".py"

# Configurables set by values in phate.config input file
DEFAULT_PIPELINE_INPUT_DIR  = 'PipelineInput/'                        # Default
DEFAULT_PIPELINE_OUTPUT_DIR = 'PipelineOutput/'                       # Default
PIPELINE_INPUT_DIR          = BASE_DIR + DEFAULT_PIPELINE_INPUT_DIR   # Default
PIPELINE_OUTPUT_DIR         = BASE_DIR + DEFAULT_PIPELINE_OUTPUT_DIR  # Default
PIPELINE_OUTPUT_SUBDIR      = 'unknown'                               # Will be read in from phate.config file 
GENOME_FILE                 = 'unknown'                               # Will be read in from phate.config file
GENOME_DIR                  = PIPELINE_INPUT_DIR                      # Default; can be overridden via phate.config
PSAT_FILE                   = 'unknown'                               # Will be read in from phate.config file, if specified
PSAT_DIR                    = PIPELINE_INPUT_DIR                      # Default; can be overridden via phate.config 

# In/out files 

logfile = PIPELINE_OUTPUT_DIR + CODE_BASE + ".log" #*** Should be converted to a generic pipeline log
outfile = PIPELINE_OUTPUT_DIR + CODE_BASE + ".out"
infile  = CONFIG_FILE   # by default, but this is specified by input parameter
LOGFILE = open(logfile,"a")  # running log; nead to clean it out occasionally
LOGFILE.write("%s%s%s\n" % ("Begin log file ",datetime.datetime.now(), " **************************** "))
OUTFILE = open(outfile,"a")  # eventually this file will contain instructions for where to find the various outputs from each module
OUTFILE.write("%s%s%s\n" % ("Begin out file ",datetime.datetime.now(), " **************************** "))

##### PATTERNS and CONTROL

# Patterns
p_comment               = re.compile("^#")
p_blank                 = re.compile("^$")
p_help                  = re.compile("help")
p_input                 = re.compile("input")
p_usage                 = re.compile("usage")
p_detail                = re.compile("detail")
p_config                = re.compile("config")
p_outputSubdir          = re.compile("output_subdir='(.*)'")
p_genomeFile            = re.compile("genome_file='(.*)'")
p_geneCaller            = re.compile("gene_caller='(.*)'")
p_psatFile              = re.compile("psat_file='(.*)'")
p_geneticCode           = re.compile("genetic_code='(\d+)'")
p_translateOnly         = re.compile("translate_only='(.*)'")
p_genomeType            = re.compile("genome_type='(.*)'")
p_name                  = re.compile("name='(.*)'")  
p_contig                = re.compile("contig='(.*)'")  #*** For now, finished genome, single contig only
p_species               = re.compile("species='(.*)'")
p_blastpIdentity        = re.compile("blast_identity='(\d+)'")   #*** For now; but should distinguish between blastn/blastp
p_blastpHitCount        = re.compile("blastp_hit_count='(\d+)'")
p_blastnHitCount        = re.compile("blastn_hit_count='(\d+)'")
p_ncbiVirusBlast        = re.compile("ncbi_virus_blast='(.*)'")
p_ncbiVirusProteinBlast = re.compile("ncbi_virus_protein_blast='(.*)'")
p_keggVirusBlast        = re.compile("kegg_virus_blast='(.*)'")
p_nrBlast               = re.compile("nr_blast='(.*)'")
p_refseqProteinBlast    = re.compile("refseq_protein_blast='(.*)'")
p_refseqGeneBlast       = re.compile("refseq_gene_blast='(.*)'")
p_swissprotBlast        = re.compile("swissprot_blast='(.*)'")
p_phantomeBlast         = re.compile("phantome_blast='(.*)'")
p_pvogsBlast            = re.compile("pvogs_blast='(.*)'")
p_uniparcBlast          = re.compile("uniparc_blast='(.*)'")
p_refseqGeneBlast       = re.compile("refseq_gene_blast='(.*)'")
p_genemarksCalls        = re.compile("genemarks_calls='(.*)'")
p_glimmerCalls          = re.compile("glimmer_calls='(.*)'")
p_prodigalCalls         = re.compile("prodigal_calls='(.*)'")
p_theaCalls             = re.compile("thea_calls='(.*)'")
p_psatAnnotation        = re.compile("psat_annotation='(.*)'")

# Booleans 

# Verbosity
CHATTY = True
#CHATTY = False

DEBUG = True
#DEBUG = False

# Other boolean control
PSAT = False            # This turns True if a psat output file is specified in the config file AND TRANSLATE_ONLY is False
TRANSLATE_ONLY = False  # User will specify 'True' in config file if only generating gene & protein files

##### HELP STRINGS

HELP_STRING = """This code, """ + CODE + """, runs a phage annotation pipeline, comprising 1) gene calling by 4 gene callers (THEA, GeneMarkS, Glimmer3, and Prodigal), followed by identification of closest phage genome by means of blast against an NCBI-phage database, and sequence-based functional annotation by means of blastp against several peptide databases (NR, KEGG-phage, and Phantome). If a PSAT output file is provided, then those annotations are merged with the blast results.\nType: python """ + CODE + """ usage - for more information about constructing the command line.\nType: python """ + CODE + """ detail - for more information about how this code can be run.\n"""

INPUT_STRING = """The input files and other parameters for running this code are specified in a configuration file, which is provided as the only input parameter. See sample configuration file (phate.config.sample) for details on how to set up the configuration file.\n"""

USAGE_STRING = """Usage: python """ + CODE + """ phate.config\n"""

DETAIL_STRING = """Currently the PSAT module is run separately as a web service. In order to incorporate PSAT output into your annotations, you should first run this pipeline specifying "translation_only" in the configuration file. Then, use the generated peptide/protein fasta file as input for PSAT processing. Once you have the PSAT output, save it to the pipeline input directory, and re-run this pipeline, specifying that translation_only is false.\n"""

##### GET INPUT PARAMETERS #####

if len(sys.argv) != 2:
    print HELP_STRING
    dateTime = os.popen('date')
    LOGFILE.write("%s%s%s%s\n" % ("Incorrect number of input parameters: ", len(sys.argv), ". End log ",dateTime))
    LOGFILE.close(); exit(0)
else:
    match_config = re.search(p_config,sys.argv[1])
    if match_config:
        configFile = sys.argv[1]
        LOGFILE.write("%s%s\n" % ("Config file is ",configFile))
    else: 
        match_input  = re.search(p_input,  sys.argv[1].lower())
        match_usage  = re.search(p_usage,  sys.argv[1].lower())
        match_detail = re.search(p_detail, sys.argv[1].lower())
        if match_input:
            print INPUT_STRING
        elif match_usage:
            print USAGE_STRING
        elif match_detail:
            print DETAIL_STRING
        else:
            print HELP_STRING
        LOGFILE.write("%s%s\n" % ("A help string was provided to user; End log ",datetime.datetime.now()))
        LOGFILE.close(); exit(0)

# Open and check input file

fileError = False
try:
    CONFIG = open(configFile,"r")
except IOError as e:
    fileError = True
    print e

if fileError:
    print "Check your config file."
    print HELP_STRING
    LOGFILE.write("%s%s\n" % ("A help string was provided to user; End log ",datetime.datetime.now()))
    LOGFILE.close(); exit(0)

##### Read input parameters from configuration file

# First, set as defaults; note: setting these values in config file is optional

geneticCode           = GENETIC_CODE
geneCaller            = GENE_CALLER
genomeType            = GENOME_TYPE
name                  = NAME
contigName            = CONTIG_NAME
species               = SPECIES
blastpIdentity        = BLASTP_IDENTITY_DEFAULT
blastpHitCount        = BLASTP_HIT_COUNT_DEFAULT
blastnHitCount        = BLASTN_HIT_COUNT_DEFAULT
ncbiVirusBlast        = NCBI_VIRUS_BLAST_DEFAULT
ncbiVirusProteinBlast = NCBI_VIRUS_PROTEIN_BLAST_DEFAULT
keggVirusBlast        = KEGG_VIRUS_BLAST_DEFAULT
nrBlast               = NR_BLAST_DEFAULT
refseqProteinBlast    = REFSEQ_PROTEIN_BLAST_DEFAULT
refseqGeneBlast       = REFSEQ_GENE_BLAST_DEFAULT
phantomeBlast         = PHANTOME_BLAST_DEFAULT
pvogsBlast            = PVOGS_BLAST_DEFAULT
uniparcBlast          = UNIPARC_BLAST_DEFAULT
swissprotBlast        = SWISSPROT_BLAST_DEFAULT
genemarksCalls        = GENEMARKS_CALLS_DEFAULT
prodigalCalls         = PRODIGAL_CALLS_DEFAULT
glimmerCalls          = GLIMMER_CALLS_DEFAULT
theaCalls             = THEA_CALLS_DEFAULT
psatAnnotation        = PSAT_ANNOTATION_DEFAULT

# Capture user's configured values

cLines = CONFIG.read().splitlines()
for cLine in cLines:
    match_comment               = re.search(p_comment,cLine)
    match_blank                 = re.search(p_blank,cLine)
    match_outputSubdir          = re.search(p_outputSubdir,cLine)
    match_genomeFile            = re.search(p_genomeFile,cLine)
    match_psatFile              = re.search(p_psatFile,cLine)
    match_geneticCode           = re.search(p_geneticCode,cLine)
    match_translateOnly         = re.search(p_translateOnly,cLine)
    match_geneCaller            = re.search(p_geneCaller,cLine)
    match_genomeType            = re.search(p_genomeType,cLine)
    match_name                  = re.search(p_name,cLine)
    match_contig                = re.search(p_contig,cLine)
    match_species               = re.search(p_species,cLine)
    match_blastpIdentity        = re.search(p_blastpIdentity,cLine)
    match_blastpHitCount        = re.search(p_blastpHitCount,cLine)
    match_blastnHitCount        = re.search(p_blastnHitCount,cLine)
    match_ncbiVirusBlast        = re.search(p_ncbiVirusBlast,cLine)
    match_ncbiVirusProteinBlast = re.search(p_ncbiVirusProteinBlast,cLine)
    match_keggVirusBlast        = re.search(p_keggVirusBlast,cLine)
    match_nrBlast               = re.search(p_nrBlast,cLine)
    match_refseqProteinBlast    = re.search(p_refseqProteinBlast,cLine)
    match_refseqGeneBlast       = re.search(p_refseqGeneBlast,cLine)
    match_phantomeBlast         = re.search(p_phantomeBlast,cLine)
    match_pvogsBlast            = re.search(p_pvogsBlast,cLine)
    match_uniparcBlast          = re.search(p_uniparcBlast,cLine)
    match_swissprotBlast        = re.search(p_swissprotBlast,cLine)
    match_refseqGeneBlast       = re.search(p_refseqGeneBlast,cLine)
    match_genemarksCalls        = re.search(p_genemarksCalls,cLine)
    match_prodigalCalls         = re.search(p_prodigalCalls,cLine)
    match_glimmerCalls          = re.search(p_glimmerCalls,cLine)
    match_theaCalls             = re.search(p_theaCalls,cLine)
    match_psatAnnotation        = re.search(p_psatAnnotation,cLine)
 
    if (match_comment or match_blank):
        continue 

    elif match_outputSubdir: #*** Note that if the output dir is not read before subdir; depends on user not changing order in config - Clean this up!
        value = match_outputSubdir.group(1)
        if value != '':
            PIPELINE_OUTPUT_SUBDIR = PIPELINE_OUTPUT_DIR + value 
        LOGFILE.write("%s%s\n" % ("PIPELINE_OUTPUT_SUBDIR is ",PIPELINE_OUTPUT_SUBDIR))

    elif match_genomeFile:
        value = match_genomeFile.group(1)
        if value != '':
            GENOME_FILE = value 
        LOGFILE.write("%s%s\n" % ("GENOME_FILE is ",GENOME_FILE))

    elif match_psatFile:
        value = match_psatFile.group(1)
        if value != '':
            PSAT_FILE = value 
            PSAT = True   # Yes, a psat file will be passed to subordinate code
        LOGFILE.write("%s%s\n" % ("PSAT_FILE is ",PSAT_FILE))

    elif match_geneticCode:
        value = match_geneticCode.group(1)
        if value != '':
            geneticCode = value

    elif match_translateOnly:
        value = match_translateOnly.group(1)
        if value.lower() == 'yes' or value.lower() == 'true' or value.lower() == 'on':
            TRANSLATE_ONLY = True
        elif value.lower() == 'no' or value.lower() == 'false' or value.lower() == 'off' or value == '':
            TRANSLATE_ONLY = False
        else:
            print "Invalid string following translate_only parameter in config file:", value
            LOGFILE.write("%s%s\n" % ("Invalid string following translate_only parameter in config file: ", value))

    elif match_geneCaller:
        value = match_geneCaller.group(1)
        if value.lower() == 'thea':
            geneCaller = 'thea'
            CONSENSUS_CALLS_FILE = 'thea.cgc'
        elif value.lower() == 'consensus':
            geneCaller = 'consensus'
            CONSENSUS_CALLS_FILE = 'consensus.cgc'
        elif value.lower() == 'genemarks' or value.lower() == 'genemark':
            geneCaller = 'genemarks'
            CONSENSUS_CALLS_FILE = 'genemark.cgc'
        elif value.lower() == 'glimmer2':
            geneCaller = 'glimmer2'
            CONSENSUS_CALLS_FILE = 'glimmer.cgc'
        elif value.lower() == 'glimmer3' or value.lower() == 'glimmer':
            geneCaller = 'glimmer3'
            CONSENSUS_CALLS_FILE = 'glimmer.cgc'
        elif value.lower() == 'prodigal':
            geneCaller = 'prodigal'
            CONSENSUS_CALLS_FILE = 'prodigal.cgc'
        elif value.lower() == 'rast':
            geneCaller = 'rast'
            CONSENSUS_CALLS_FILE = 'rast.cgc'
        print "Gene caller has been set to:", geneCaller

    elif match_genomeType:
        value = match_genomeType.group(1)
        if value.lower() == 'phage' or value.lower() == 'bacteriophage':
            genomeType = 'phage' 
        elif value.lower() == 'virus' or value.lower() == 'viral' or value.lower() == 'viridae':
            genomeType = 'virus'
        elif value.lower() == 'bacteria' or value.lower() == 'bacterium' or value.lower() == 'bacterial':
            genomeType = 'bacterium' 

    elif match_name:
        value = match_name.group(1)
        name = value

    elif match_contig:
        value = match_contig.group(1)
        contigName = value

    elif match_species:
        value = match_species.group(1)
        species = value

    elif match_blastpIdentity:
        value = match_blastpIdentity.group(1)
        if int(value) > int(MIN_BLASTP_IDENTITY) and int(value) <= 100:
            blastpIdentity = value

    elif match_blastpHitCount:
        value = match_blastpHitCount.group(1)
        if int(value) > 0 and int(value) <= MAX_BLASTP_HIT_COUNT:
            blastpHitCount = value

    elif match_blastnHitCount:
        value = match_blastnHitCount.group(1)
        if int(value) > 0 and int(value) <= MAX_BLASTN_HIT_COUNT:
            blastnHitCount = value

    elif match_ncbiVirusBlast:
        value = match_ncbiVirusBlast.group(1)
        if value.lower() == 'true' or value.lower() == 'yes' or value.lower() == 'on':
            ncbiVirusBlast = True
        else:
            ncbiVirusBlast = False

    elif match_ncbiVirusProteinBlast:
        value = match_ncbiVirusProteinBlast.group(1)
        if value.lower() == 'true' or value.lower() == 'yes' or value.lower() == 'on':
            ncbiVirusProteinBlast = True
        else:
            ncbiVirusProteinBlast = False

    elif match_keggVirusBlast:
        value = match_keggVirusBlast.group(1)
        if value.lower() == 'true' or value.lower() == 'yes' or value.lower() == 'on':
             keggVirusBlast = True
        else:
             keggVirusBlast = False

    elif match_nrBlast:
        value = match_nrBlast.group(1)
        if value.lower() == 'true' or value.lower() == 'yes' or value.lower() == 'on':
             nrBlast = True
        else:
             nrBlast = False 

    elif match_refseqProteinBlast:
        value = match_refseqProteinBlast.group(1)
        if value.lower() == 'true' or value.lower() == 'yes' or value.lower() == 'on':
            refseqProteinBlast = True
        else:
            refseqProteinBlast = False

    elif match_refseqGeneBlast:
        value = match_refseqGeneBlast.group(1)
        if value.lower() == 'true' or value.lower() == 'yes' or value.lower() == 'on':
            refseqGeneBlast = True
        else:
            refseqGeneBlast = False

    elif match_phantomeBlast:
        value = match_phantomeBlast.group(1)
        if value.lower() == 'true' or value.lower() == 'yes' or value.lower() == 'on':
            phantomeBlast = True
        else:
            phantomeBlast = False 

    elif match_pvogsBlast:
        value = match_pvogsBlast.group(1)
        if value.lower() == 'true' or value.lower() == 'yes' or value.lower() == 'on':
            pvogsBlast = True
        else:
            pvogsBlast = False

    elif match_uniparcBlast:
        value = match_uniparcBlast.group(1).lower()
        if value == 'true' or value == 'yes' or value == 'on':
            uniparcBlast = True
        else:
            uniparcBlast = False

    elif match_swissprotBlast:
        value = match_swissprotBlast.group(1).lower()
        if value == 'true' or value == 'yes' or value == 'on':
            swissprotBlast = True
        else:
            swissprotBlast = False

    elif match_genemarksCalls:
        value = match_genemarksCalls.group(1)
        if value.lower() == 'true' or value.lower() == 'yes' or value.lower() == 'on':
            genemarksCalls = True
        else:
            genemarksCalls = False

    elif match_prodigalCalls:
        value = match_prodigalCalls.group(1)
        if value.lower() == 'true' or value.lower() == 'yes' or value.lower() == 'on':
            prodigalCalls = True
        else:
            prodigalCalls = False

    elif match_glimmerCalls:
        value = match_glimmerCalls.group(1)
        if value.lower() == 'true' or value.lower() == 'yes' or value.lower() == 'on':
            glimmerCalls = True
        else:
            glimmerCalls = False

    elif match_theaCalls:
        value = match_theaCalls.group(1)
        if value.lower() == 'true' or value.lower() == 'yes' or value.lower() == 'on':
            theaCalls = True
        else:
            theaCalls = False

    elif match_psatAnnotation:
        value = match_psatAnnotation.group(1)
        if value.lower() == 'true' or value.lower() == 'yes' or value.lower() == 'on':
            psatAnnotation = True
        else:
            psatAnnotation = False

    else:
        LOGFILE.write("%s%s\n" % ("ERROR: Unrecognized line in config file: ", cLine))
        print "ERROR: unrecognized line in config file:", cLine

if DEBUG:
    print "After reading config file, contigName is", contigName
    print "After reading config file, genemarksCalls is", genemarksCalls
    print "After reading config file, ncbiVirusProteinBlast is", ncbiVirusProteinBlast

# Create objects for passing blast and genecall parameters to subordinate codes 

blastParameters = {
    "ncbiVirusBlast"        : ncbiVirusBlast,
    "ncbiVirusProteinBlast" : ncbiVirusProteinBlast,
    "keggVirusBlast"        : keggVirusBlast,
    "nrBlast"               : nrBlast,
    "refseqProteinBlast"    : refseqProteinBlast,
    "refseqGeneBlast"       : refseqGeneBlast,
    "phantomeBlast"         : phantomeBlast,
    "pvogsBlast"            : pvogsBlast,
    "uniparcBlast"          : uniparcBlast,
    "swissprotBlast"        : swissprotBlast,
    }

genecallParameters = {
    "genemarksCalls"     : genemarksCalls,
    "prodigalCalls"      : prodigalCalls,
    "glimmerCalls"       : glimmerCalls,
    "theaCalls"          : theaCalls,
    }

# Double check: issue warning if necessary, but continue processing assuming this is what the user intends.
if GENOME_TYPE == 'PHAGE' and CONSENSUS_CALLS_FILE != 'thea.cgc':
    print "WARNING: If genome type is phage, the consensus gene-call file should be thea.cgc! Yours is", CONSENSUS_CALLS_FILE
    LOGFILE.write("%s%s\n" % ("WARNING:  User has selected genome type as phage, but consensus gene-call file as ", CONSENSUS_CALLS_FILE))

CONFIG.close()

# Turn PSAT back to False if the user only needs translations (even if they provided a .psat file)
if PSAT and TRANSLATE_ONLY:  
    PSAT = False

if CHATTY:
    print "PIPELINE_INPUT_DIR is", PIPELINE_INPUT_DIR
    print "PIPELINE_OUTPUT_DIR is", PIPELINE_OUTPUT_DIR
    print "PIPELINE_OUTPUT_SUBDIR is", PIPELINE_OUTPUT_SUBDIR
    print "GENOME_FILE is", GENOME_FILE
    print "GENE_FILE is", GENE_FILE
    print "PROTEIN_FILE is", PROTEIN_FILE
    print "geneticCode is", geneticCode 
    print "TRANSLATE_ONLY is", TRANSLATE_ONLY
    if PSAT:
        print "PSAT_FILE is", PSAT_FILE
    else:
        print "PSAT_FILE was not provided." 
    print "geneCaller is", geneCaller 
    print "CONSENSUS_CALLS_FILE is", CONSENSUS_CALLS_FILE
    print "genomeType is", genomeType 
    print "name is", name 
    print "contigName is", contigName
    print "species is", species 
    print "blastpIdentity is", blastpIdentity 
    print "blastpHitCount is", blastpHitCount 
    print "blastnHitCount is", blastnHitCount
    print "ncbiVirusBlast is", ncbiVirusBlast
    print "ncbiVirusProteinBlast is", ncbiVirusProteinBlast
    print "keggVirusBlast is", keggVirusBlast
    print "nrBlast is", nrBlast
    print "refseqProteinBlast is", refseqProteinBlast
    print "refseqGeneBlast is", refseqGeneBlast
    print "phantomeBlast is", phantomeBlast
    print "pvogsBlast is", pvogsBlast
    print "uniparcBlast is", uniparcBlast
    print "swissprotBlast is", swissprotBlast
    print "genemarksCalls is", genemarksCalls
    print "prodigalCalls is", prodigalCalls
    print "glimmerCalls is", glimmerCalls
    print "theaCalls is", theaCalls
    print "psatAnnotation is", psatAnnotation

LOGFILE.write("%s\n" % ("Input parameters:"))
LOGFILE.write("%s%s\n" % ("   PIPELINE_INPUT_DIR: ", PIPELINE_INPUT_DIR))
LOGFILE.write("%s%s\n" % ("   PIPELINE_OUTPUT_DIR: ", PIPELINE_OUTPUT_DIR))
LOGFILE.write("%s%s\n" % ("   PIPELINE_OUTPUT_SUBDIR: ", PIPELINE_OUTPUT_SUBDIR))
LOGFILE.write("%s%s\n" % ("   GENOME_FILE: ", GENOME_FILE))
LOGFILE.write("%s%s\n" % ("   GENE_FILE: ", GENE_FILE))
LOGFILE.write("%s%s\n" % ("   PROTEIN_FILE: ", PROTEIN_FILE))
LOGFILE.write("%s%s\n" % ("   PSAT_FILE: ", PSAT_FILE))
LOGFILE.write("%s%s\n" % ("   geneticCode: ", geneticCode))
LOGFILE.write("%s%s\n" % ("   Status of boolean TRANSLATE_ONLY is ",TRANSLATE_ONLY))
LOGFILE.write("%s%s\n" % ("   Status of boolean PSAT is ",PSAT))
LOGFILE.write("%s%s\n" % ("   geneCaller is ",geneCaller))
LOGFILE.write("%s%s\n" % ("   CONSENSUS_CALLS_FILE is ",CONSENSUS_CALLS_FILE))
LOGFILE.write("%s%s\n" % ("   genomeType is ",genomeType))
LOGFILE.write("%s%s\n" % ("   name is ",name))
LOGFILE.write("%s%s\n" % ("   contigName is ",contigName))
LOGFILE.write("%s%s\n" % ("   species is ",species))
LOGFILE.write("%s%s\n" % ("   blastpIdentity is ",blastpIdentity))
LOGFILE.write("%s%s\n" % ("   blastpHitCount is ",blastpHitCount))
LOGFILE.write("%s%s\n" % ("   blastnHitCount is ",blastnHitCount))
LOGFILE.write("%s%s\n" % ("   ncbiVirusBlast is ",ncbiVirusBlast))
LOGFILE.write("%s%s\n" % ("   ncbiVirusProteinBlast is ",ncbiVirusProteinBlast))
LOGFILE.write("%s%s\n" % ("   keggVirusBlast is ",keggVirusBlast))
LOGFILE.write("%s%s\n" % ("   nrBlast is ",nrBlast))
LOGFILE.write("%s%s\n" % ("   refseqProteinBlast is ",refseqProteinBlast))
LOGFILE.write("%s%s\n" % ("   refseqGeneBlast is ",refseqGeneBlast))
LOGFILE.write("%s%s\n" % ("   phantomeBlast is ",phantomeBlast))
LOGFILE.write("%s%s\n" % ("   pvogsBlast is ",pvogsBlast))
LOGFILE.write("%s%s\n" % ("   uniparcBlast is ",uniparcBlast))
LOGFILE.write("%s%s\n" % ("   swissprotBlast is ",swissprotBlast))
LOGFILE.write("%s%s\n" % ("   genemarksCalls is ",genemarksCalls))
LOGFILE.write("%s%s\n" % ("   prodigalCalls is ",prodigalCalls))
LOGFILE.write("%s%s\n" % ("   glimmerCalls is ",glimmerCalls))
LOGFILE.write("%s%s\n" % ("   theaCalls is ",theaCalls))
LOGFILE.write("%s%s\n" % ("   psatAnnotation is ",psatAnnotation))

# Create user's output subdirectory, if doesn't already exist

try:
    os.stat(PIPELINE_OUTPUT_SUBDIR)
except:
    os.mkdir(PIPELINE_OUTPUT_SUBDIR)

# Open and check input file(s)

inputDir     = PIPELINE_INPUT_DIR
genomeFile   = PIPELINE_INPUT_DIR + GENOME_FILE
genecallFile = PIPELINE_OUTPUT_SUBDIR + CONSENSUS_CALLS_FILE
geneFile     = PIPELINE_OUTPUT_SUBDIR + GENE_FILE
proteinFile  = PIPELINE_OUTPUT_SUBDIR + PROTEIN_FILE
if PSAT:
    psatFile = PIPELINE_INPUT_DIR + PSAT_FILE
else:
    psatFile = ''
outputDir    = PIPELINE_OUTPUT_SUBDIR

LOGFILE.write("%s%s\n" % ("inputDir is ",inputDir))
LOGFILE.write("%s%s\n" % ("genomeFile is ",genomeFile))
LOGFILE.write("%s%s\n" % ("genecallFile is ",genecallFile))
LOGFILE.write("%s%s\n" % ("geneFile is ",geneFile))
LOGFILE.write("%s%s\n" % ("proteinFile is ",proteinFile))
LOGFILE.write("%s%s\n" % ("psatFile is ",psatFile))

# Check PSAT file

LOGFILE.write("%s\n" % ("Checking files..."))
fileError = False
if PSAT:
    try:
        PSAT_H = open(psatFile,"r")
    except IOError as e:
        fileError = True
        print e

    if fileError:
        print "Check your PSAT file,", psatFile
        print USAGE_STRING
        LOGFILE.write("%s%s%s%s\n" % ("ERROR:  PSAT file could not be opened: ", psatFile, "; End log ", datetime.datetime.now()))
        LOGFILE.close(); exit(0)
    PSAT_H.close()

# Check genome file

fileError = False
try:
    GENOME_H = open(genomeFile,"r")
except IOError as e:
    fileError = True
    print e 

if fileError:
    print USAGE_STRING
    print "Check your genome file,", genomeFile
    LOGFILE.write("%s%s%s%s\n" % ("ERROR:  Genome file could not be opened: ", genomeFile, "; End log ", datetime.datetime.now()))
    LOGFILE.close(); exit(0)
GENOME_H.close()

# Copy config file to the user's results directory
configSave = outputDir + configFile
command = "cp " + configFile + ' ' + configSave
os.system(command)

##### BEGIN MAIN ########################################################################################

##### Run Gene-calling Module

LOGFILE.write("%s\n" % ("Preparing to run genecall module..."))

param2 = outputDir[:-1]  # remove terminal '/' because subordinate code adds it explicitly

param3 = ''  # can't pass a dict; future re-write genecalling module as class
if genecallParameters["genemarksCalls"]:
    param3 += "genemarks_"
if genecallParameters["prodigalCalls"]:
    param3 += "prodigal_"
if genecallParameters["glimmerCalls"]:
    param3 += "glimmer_"
if genecallParameters["theaCalls"]:
    param3 += "thea_"

command = "python " + GENECALL_CODE + ' ' + genomeFile + ' ' + param2 + ' ' + param3 

print "Calling the gene-call module...command is,", command
LOGFILE.write("%s%s\n" % ("Calling the gene-call module. Command is ", command))
result = os.system(command)
if CHATTY:
    print "Done!"
LOGFILE.write("%s%s\n" % ("Gene-call processing complete at ", datetime.datetime.now()))

##### Run Sequence Annotation Module

LOGFILE.write("%s\n" % ("Preparing to call sequence annotation module..."))
if DEBUG:
    print "Before constructing command line to invoke sequence annotation code, contigName is", contigName

# Construct command line parameter string

# First, construct string listing the names of databases to be blasted
blastParameterString = '_'
if blastParameters['ncbiVirusBlast']:
    blastParameterString += '_ncbiVirusGenome'
if blastParameters['ncbiVirusProteinBlast']:
    blastParameterString += '_ncbiVirusProtein'
if blastParameters['nrBlast']:
    blastParameterString += '_nr'
if blastParameters['keggVirusBlast']:
    blastParameterString += '_kegg'
if blastParameters['refseqProteinBlast']:
    blastParameterString += '_refseqProtein'
if blastParameters['refseqGeneBlast']:
    blastParameterString += '_refseqGene'
if blastParameters['pvogsBlast']:
    blastParameterString += '_pvogs'
if blastParameters['phantomeBlast']:
    blastParameterString += '_phantome'
if blastParameters['uniparcBlast']:
    blastParameterString += '_uniparc'
if blastParameters['swissprotBlast']:
    blastParameterString += '_swissprot'

commandRoot1 = "python " + SEQANNOTATION_CODE + " -o " + outputDir  # code and output direction
commandRoot2 = " -G " + genomeFile        + " -g " + geneFile       + " -p " + proteinFile    # genome files
commandRoot3 = " -c " + geneCaller        + " -f " + genecallFile   + " -C " + contigName     # gene-call information
commandRoot4 = " -t " + genomeType        + " -n " + name           + " -s " + species        # genome meta-data
commandRoot5 = " -i " + blastpIdentity    + " -h " + blastpHitCount + " -H " + blastnHitCount # blast parameters
commandRoot6 = " -d " + blastParameterString                                                  # databases to blast against
commandRoot = commandRoot1 + commandRoot2 + commandRoot3 + commandRoot4 + commandRoot5 + commandRoot6

# As appropriate, append additional parameters
if TRANSLATE_ONLY:
    command = commandRoot + " -x true "        # setting TRANSLATE_ONLY to True
elif PSAT:
    command = commandRoot + " -P " + psatFile  # including a PSAT results file
else:
    command = commandRoot 

# Communicate and execute
if CHATTY:
    print "Calling the sequence annotation module...command is,", command
LOGFILE.write("%s%s\n" % ("Calling the sequence annotation module. Command is ", command))
result = os.system(command)
if CHATTY:
    print "Done!"
LOGFILE.write("%s%s\n" % ("Sequence annotation processing complete at ", datetime.datetime.now()))

##### CLEAN UP

if CHATTY:
    print "Code completed at", datetime.datetime.now()
OUTFILE.write("%s%s\n" %("Pipeline output is in output file created by code ",SEQANNOTATION_CODE))
OUTFILE.close()
LOGFILE.write("%s%s\n" % ("Code completed at ", datetime.datetime.now()))
LOGFILE.close()
