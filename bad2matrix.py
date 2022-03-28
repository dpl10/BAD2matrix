import sys
import os
import re
import warnings
from utils import Term_data, Partition, clean_name, get_name_map

in_dir = ""
root_name = ""
full_fasta_names = False
infiles = []
term_names = []
name_map = {}
code_indels = True
outgroups = []
aa_encoding = "20"
code_gene_content = True
missing_percentile = 1
part_file = 0
tsv_file = None
debbug = False

for iar,ar in enumerate(sys.argv):

	if ar == '-a':
		if re.search(r'^(2|3|4|5|6|6dso|6kgb|6sr|8|10|11|12|15|18|20)$', sys.argv[iar+1]):
			aa_encoding = 20

	elif ar == '-d':
		if os.path.exists(sys.argv[iar+1]):
			in_dir = sys.argv[iar+1]
		else:
			raise ValueError("Input directory (-d) could not be read!")

	elif ar == '-f':
		full_fasta_names = True

	elif ar == '-g':
		code_gene_content = False

	elif ar == '-i':
		code_indels = False

	elif ar == '-m':
		val = int(re.sub(r'\D', '', sys.argv[iar+1])) 
		if 0 < val <= 100:
			missing_percentile = val

	elif ar == '-n':
		root_name = sys.argv[iar+1]

	elif ar == '-o':
		outgroups = [clean_name(x) for x in sys.argv[iar+1].split(',')]

	elif ar == '-r':
		part_file = 1

	elif ar == '-t':
		if os.path.exists(sys.argv[iar+1]):
			tsv_file = sys.argv[iar+1]
		else:
			raise ValueError("Input tsv file (-t) could not be read!")

	elif ar == '-debbug':
		debbug = True


if in_dir:
	for d, s, f in os.walk(in_dir):
		for file in f:
			infiles.append(os.path.join(d,file))

	if len(infiles) == 0:
		raise ValueError("Input directory (-d) does not contain any files!")

if len(infiles) > 0 and len(root_name) > 0:
	
	name_map = get_name_map(infiles, full_fasta_names)
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
	[ -o speciesA,speciesB... ] [ -r ] [ -t ]

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

-t  Data matrix in tsv format, encapsulating ortholog duplication encoding.

	""")


exit()