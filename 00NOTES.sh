# EPN, Thu Mar 26 15:50:58 2015
# version 0.01
# 
# Example of running esl-fetch-cds.pl to fetch
# CDS sequences. 
#
# This is a Perl that uses the Bio-Easel library which 
# uses Inline to interact and call functions from 
# the Easel sequence analysis library from Sean Eddy's
# group. 
#
#
#######################
# CODE-RELATED TODOs
# - clean up code
# - add option to specify path to idfetch (currently hardcoded)
# - see additional TODOs in the code (esl-fetch-cds.pl)
# - write a test script
# 
#######################
# More information
#######################
#
# See /home/nawrocke/notebook/15_0324_dnaorg_esl_fetch_cds/00LOG.txt
# for notes on development and testing of this script.
# 
#######################
# Prerequisites
#######################
# 
# Directories that include the BioEasel perl modules must be part of your
#  $PERL5LIB environment variable.
# 
# To modify your PERL5LIB environment variable appropriately:
# For bash shell users
source /panfs/pan1/dnaorg/programs/setup-bio-easel.bash.sh
# For C shell or C shell compatible users
source /panfs/pan1/dnaorg/programs/setup-bio-easel.csh.sh

#######################
# Usage and options
#######################
#esl-fetch-cds.pl [OPTIONS] <input coordinate file>
#	OPTIONS:
#		-odir <s>: write temporary files to directory <s>, instead of cwd
#		-onlyaccn: name sequences only as the accession (no coord info)
#		-nocodon : do not include codonstart in sequence name
#
############################
# Example command and output
############################
# 
esl-fetch-cds.pl sample.in > sample.cds.fa
#
# This will create a fasta file 'sample.cds.fa' with the CDS sequences
# specified in smn1.in, in format described above. These sequences will be named
# as follows:
# ><token1>:codon_start<codon_start>:<accn>:<start>:<stop>:<strand>
# 
# With values taken from the input file. <token1> is 
# the first token from the input coordinate file. 
# <accn>:<start><stop> are derived from similar
# tokens in the input file, and <strand> will
# be '-' if 'complement' exists in the corresponding
# input file line, and otherwise '+'.
# 
# For example for the input file sample.in:
# 
# AAH45158.1	BC045158.1:4..870
# CAA91173.1	complement(join(CU329670.1:822332..822429,CU329670.1:822489..822849))
# AFJ92633.1	JQ690867.1:<1..400	2
# BAJ21116.1	AB590961.1:10..>891
# 
# The four sequences are named:
# 
# >AAH45158.1:codon_start1:BC045158.1:4:870:+:
# >CAA91173.1:codon_start1:CU329670.1:822332:822429:-:CU329670.1:822489:822849:-:
# >AFJ92633.1:codon_start2:JQ690867.1:<1:400:+:
# >BAJ21116.1:codon_start1:AB590961.1:10:>891:+:
#
######################
# Input files
######################
# 
# AAH45158.1	BC045158.1:4..870
# CAA91173.1	complement(join(CU329670.1:822332..822429,CU329670.1:822489..822849))
# AFJ92633.1	JQ690867.1:<1..400	2
# BAJ21116.1	AB590961.1:10..>891
#
# The first token is (in this case) the protein accession the CDS
# derives from, it is used as part of the name of the eventual output
# CDS sequence only.
#
# Token 2 is from a 'coded_by' INSDQualifier field of a 'CDS'
# INSDFeature_key in an NCBI 'gpc' format file for a protein
# sequence. 
#
# There is an optional third token, which is present in one of the
# sample input lines. This token specifies the 'codon_start'
# INSDQualifier field of a 'CDS' INSDFeature_key the position at which
# to start translating the CDS. By default, this value is '1', so for
# lines without a third token this value is set to '1'.
#
# Two examples of fetching such information using edirect
# tools:
#
# 'esearch -db protein -query AAH45158.1 | efetch -format gpc | xtract -insd CDS coded_by codon_start' 
# AAH45158.1	BC045158.1:4..870
#
# 'esearch -db protein -query AFJ92633.1 | efetch -format gpc | xtract -insd CDS coded_by codon_start' 
# AFJ92633.1	JQ690867.1:<1..400	2
#
# The second token completely specifies how to build a CDS from >=
# 1 accession.versions. It is of the following format:
# 
# 0 or 1 of 'complement(' : if present, the CDS is on the reverse strand 
# 0 or 1 of 'join('       : multi-segment, more than one <accn>:<start>..<stop> follow
# >= 1 of '<accn>:<start>..<stop>', ' : source accession, start position, stop position 
#
# With an optional '<' character before <start> and
# an optional '>' character before the <end>. These
# characters indicate whether the CDS information
# is 'incomplete' at the start or end.
# 
# Ref:  http://www.ncbi.nlm.nih.gov/projects/Sequin/table.html
# 
# From that page:
# ~~~~~~~~~~~~~~~~~~
# Locations of partial (incomplete) features are indicated with a ">" or
# "<" next to the number. In this example, the first gene, CDS, and mRNA
# all begin upstream of the start of the nucleotide sequence. The "<"
# symbol indicates that they are 5' partial features. Furthermore, for
# the protein to translate correctly, the correct reading frame must be
# indicated with the qualifier "codon_start" on the first CDS. There is
# no need to indicate the codon_start on complete CDSs, as it is assumed
# that the translation starts at the first nucleotide of the interval if
# no codon_start is provided.
# ~~~~~~~~~~~~~~~~~~
# 
# This format is *probably* a standard NCBI way of representing
# multiple exons or subsequences to be joined together.
#
##############################################
# Last updated: EPN, Mon Apr  6 16:05:10 2015
##############################################
# Log of changes:
##############################################
# EPN, Mon Apr 13 09:45:51 2015
#
# Added -odir <s> option to output temporary files
# to pre-existing directory <s>.
#
# previous version: ./bkups/15_0413-1-before-update/
# *this* version:   ./bkups/15_0413-2-after-update/
##############################################
# EPN, Mon Apr  6 15:59:20 2015
#
# Changed input and output format. Input now
# allows an optional third character, the codon_start
# value. Sequences are now named differently too.
# Names include the codon_start value and also
# '>' or '<' characters indicating incomplete
# CDS information.
# previous version: ./bkups/15_0406-1-before-update/
# *this* version:   ./bkups/15_0406-2-after-update/
##############################################
# EPN, Tue Mar 31 15:37:42 2015
#
# Fixed major bug fix for fetching sequences
# with >1 exons on negative strand, previously
# I wasn't reversing the order of exons,
# which is necessary. 
# previous version: ./bkups/15_0331-1-before-update/
# *this* version:   ./bkups/15_0331-2-after-update/
##############################################

