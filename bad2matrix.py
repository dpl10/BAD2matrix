import sys
import os
import re
import warnings
from utils import Term_data, Partition

in_dir = ""
root_name = ""
full_fasta_names = False
infiles = []
term_names = []
name_map = {}
code_indels = True

for iar,ar in enumerate(sys.argv):

	if ar == '-d':
		if os.path.exists(sys.argv[iar+1]):
			in_dir = sys.argv[iar+1]
		else:
			raise ValueError("Input directory (-d) could not be read!")

	if ar == '-n':
		root_name = sys.argv[iar+1]

	if ar == '-f':
		full_fasta_names = True

	if ar == '-i':
		code_indels = False

print(root_name, in_dir)

if in_dir:
	for d, s, f in os.walk(in_dir):
		for file in f:
			infiles.append(os.path.join(d,file))

	if len(infiles) == 0:
		raise ValueError("Input directory (-d) does not contain any files!")

if len(infiles) > 0 and len(root_name) > 0:
	

	for file in infiles:

		########################################################################
		# File naming convention of inputs should be stated in the instructions.
		# Two classes of input matrices: molecular or morphological/gene dupli-
		# cation events. User should mention which is which through the file
		# extension.
		######################################################################## 
		if re.search(r'\.fas?t?a?$', file):

			with open(file, 'r') as fhandle:
				thname_map = {}

				for line in fhandle:

					if line.startswith(">"):
						# All cleaning procedures of terminal names should be done here.
						raw_name = line.strip()
						name = raw_name
						if not full_fasta_names:
							name = re.split(r'#+', name)[0]
						name = re.sub(r'[\s\/\-\\]+', '_', name)
						name = re.sub(r'[^\w\._]', '', name)
						thname_map[raw_name] = name

				if len(thname_map.values()) < len(thname_map):

					thnames = list(thname_map.keys())
					dup_count = list(map(lambda x: thnames.count(x), thnames))
					prob = [x for x,y in zip(thnames, dup_count) if y > 1]
					err = '\n'.join(prob)
					raise ValueError(f"Check the following sequence names in ´{file}´, there could be duplicates:\n{err}.\n")

				name_map.update(thname_map)

	term_names = sorted(list(set(name_map.values()))) # Why sort should be done in reverse order?
	spp_data = {name: Term_data(name) for name in term_names}

	for file in infiles:
		partition = Partition(file, name_map)
		print(partition.data)

		if code_indels:
			#size = len(list(partition.values())[0])
			#part_table = [size, 'molecular']
			partition.indel_coder()
			print(partition.data)

		# Parse all data to each species file
		for name in spp_data:
			spp_data[name].feed(partition)


	"""
	for each fasta file:
		for each taxa in the total count:
			parse to the temporary file
			compute accesory data (indels, aa counts) 
	"""


else:

	print("""
A Python script for merging and translating FASTA alignments into TNT extended
XREAD, FastTree FASTA, and RAxML/IQ-Tree extended PHYLIP with indel characters
(optionally, .tnt & .phy) coded using the 'simple' gap coding method of Simmons
and Ochoterena (2000; Gaps as characters in sequence-based phylogenetic
analysis. Systematic Biology 49: 369-381. DOI 10.1080/10635159950173889), and
with coded gene content (absence/presence) characters (optionally, .tnt & .phy).
This script is slower than 2matrix.pl due to more disk use, but will not run out
of RAM (hopefully). The RAxML .part file is either in original RAxML (-r) or
RAxML-NG format.

USAGE:	BAD2matrix.pl [ -a 2|3|4|5|6|6dso|6kgb|6sr|8|10|11|12|15|18|20 ]
	-d <directory> [ -f ] [ -g ] [ -i ] [ -m int ] -n <root-name>
	[ -o speciesA,speciesB... ] [ -r ]

OPTIONS:
-a	Number of amino acid states (default = 20). Reduction with option '6dso'
	follows Dayhoff, Schwartz, and Orcutt (1978. A model of evolutionary
	change in proteins. In Atlas of Protein Sequence and Structure, Dayhoff
	ed. pp. 345–352); option '6kgb' follows Kosiol, Goldman, and Buttimore
	(2004. A new criterion and method for amino acid classification. Journal
	of Theoretical Biology 228:97–106. DOI 10.1016/j.jtbi.2003.12.010);
	option '6sr' follows Susko and Roger (2007. On reduced amino acid
	alphabets for phylogenetic inference. Molecular Biology and Evolution 24:
	2139–2150. DOI 10.1093/molbev/msm144); option '11' follows Buchfink, Xie,
	and Huson (2015. Fast and sensitive protein alignment using DIAMOND.
	Nature Methods 12: 59–60. DOI 10.1038/nmeth.3176); and all other options
	follow Murphy, Wallqvist, and Levy (2000. Simplified amino acid alphabets
	for protein fold recognition and implications for folding. Protein
	Engineering 13: 149–52. DOI 10.1093/protein/13.3.149).
-d	Specifies a directory of aligned FASTA files. By default, names should
	use the OrthologID naming convention ('>species#sequenceID').
	Characters other than letters, numbers, periods, and underscores will be
	deleted. Use -f for an alternate naming convention.
-f	Use full FASTA names. Characters other than letters, numbers, periods,
	and underscores will be deleted.
-g	Do NOT code gene content.
-i	Do NOT code indels.
-m	Retain the upper x percentile of genes in the distribution of missing
	sequences (default = 1; i.e. include all genes).
-n	<root-name> for output files.
-o	Outgroup(s) for rooting trees.
-r	Output original RAxML formatted .part file.

	""")


exit()