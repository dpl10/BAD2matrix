import sys
import os
import re
from functools import reduce
import warnings
from utils import Term_data, Partition, Polymorphs, clean_name, get_name_map,aa_redux_dict

help_text = """

A Python script for merging and translating FASTA alignments into TNT (extended
XREAD) and RAxML-NG/IQ-Tree (extended PHYLIP) input matrices. Optionally, it 
encodes indel characters using the 'simple' gap coding method of Simmons
and Ochoterena (2000; Gaps as characters in sequence-based phylogenetic
analysis. Systematic Biology 49: 369-381. DOI 10.1080/10635159950173889), and 
gene content as binary characters (absence/presence). This script is slower than 
2matrix.pl due to more disk use, but will not run out of RAM (hopefully). 

USAGE:	python bad2matrix.py [ -a 2|3|4|5|6|6dso|6kgb|6sr|8|10|11|12|15|18|20 ]
	-d <directory> [ -f ] [ -g ] [ -i ] [ -m int ] -n <root-name>
	[ -t <directory> ]

OPTIONS:
-a	Number of amino acid states (default = 20). Reduction with option '6dso'
	follows Dayhoff, Schwartz, and Orcutt (1978. A model of evolutionary
	change in proteins. In Atlas of Protein Sequence and Structure, Dayhoff
	ed. pp. 345-352); option '6kgb' follows Kosiol, Goldman, and Buttimore
	(2004. A new criterion and method for amino acid classification. Journal
	of Theoretical Biology 228:97-106. DOI 10.1016/j.jtbi.2003.12.010);
	option '6sr' follows Susko and Roger (2007. On reduced amino acid
	alphabets for phylogenetic inference. Molecular Biology and Evolution 24:
	2139-2150. DOI 10.1093/molbev/msm144); option '11' follows Buchfink, Xie,
	and Huson (2015. Fast and sensitive protein alignment using DIAMOND.
	Nature Methods 12: 59-60. DOI 10.1038/nmeth.3176); and all other options
	follow Murphy, Wallqvist, and Levy (2000. Simplified amino acid alphabets
	for protein fold recognition and implications for folding. Protein
	Engineering 13: 149-52. DOI 10.1093/protein/13.3.149).

-d	Specifies a directory of aligned FASTA files. By default, names should
	use the OrthologID naming convention ('>species#sequenceID').
	Characters other than letters, numbers, periods, and underscores will be
	deleted. Use -f for an alternate naming convention.

-f	Use full FASTA names. Characters other than letters, numbers, periods,
	and underscores will be deleted. Default is the species root name: the 
	string before the first sharp ('#').

-g	Do NOT code gene content.

-i	Do NOT code indels.

-m	Retain the upper x percentile of genes in the distribution of missing
	sequences (default = 1; i.e. include all genes).

-n	<root-name> for output files.

-t	Folder containing data matrices in tsv format, encapsulating ortholog 
	duplication encoding. Polymorphic encodings should use pipes (`|`) to 
	separate states.

"""


in_dir = ""
in_dir_morph = ""
root_name = ""
full_fasta_names = False
infiles = []
infiles_morph = []
term_names = []
name_map = {}
code_indels = True
aa_encoding = "20"
code_gene_content = True
keep_percentile = 1
part_file = 0
tsv_file = None
raxml_bffr = ""

debbug = False

# parse command line arguments
for iar,ar in enumerate(sys.argv):

	if ar == '-a':
		if re.search(r'^(2|3|4|5|6|6dso|6kgb|6sr|8|10|11|12|15|18|20)$', sys.argv[iar+1]):
			aa_encoding = sys.argv[iar+1]

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
			keep_percentile = val / 100

	elif ar == '-n':
		root_name = sys.argv[iar+1]

	elif ar == '-r':
		part_file = 1

	elif ar == '-t':
		if os.path.exists(sys.argv[iar+1]):
			in_dir_morph = sys.argv[iar+1]
		else:
			raise ValueError("Input directory (-t) could not be read!")

	elif ar == '-debbug':
		debbug = True


# Check input directory contents
if in_dir_morph:
	for d, s, f in os.walk(in_dir_morph):
		for file in f:
			infiles_morph.append(os.path.join(d,file))

	if len(infiles_morph) == 0:
		raise ValueError("Input directory (-t) does not contain any files!")

if in_dir:
	for d, s, f in os.walk(in_dir):
		for file in f:
			infiles.append(os.path.join(d,file))

	if len(infiles) == 0:
		raise ValueError("Input directory (-d) does not contain any files!")


if len(infiles) > 0 and len(root_name) > 0:

	#TODO Check if output files already exist

	translation_dict = aa_redux_dict(aa_encoding)
	raxml_main = root_name + '.phy'
	raxml_part = root_name + '.part'
	iqtree_nexus = root_name + '.nex'
	tnt_main = root_name + '.ss'
	(name_map, act_files) = get_name_map(infiles, full_fasta_names, keep_percentile, infiles_morph)
	term_names = sorted(list(set(name_map.values()))) #? Why sort should be done in reverse order?
	longest = len(max(term_names, key = len))
	spp_data = {name: Term_data(name, code_gene_content) for name in term_names}
	non_informative_partitions = []
	final_spp_count = 0
	part_collection = {'size': [], 'type': [], 'states': []}
	polys = Polymorphs()

	for file in act_files:

		partition = Partition(file, name_map, translation_dict, polys)

		if code_indels and partition.filetype == 'fasta':
			partition.indel_coder()

		partition.informative_stats()
		tot_inf = len(reduce(lambda x, y: x + y, partition.metadata["informative_chars"]))

		if tot_inf == 0:
			warnings.warn(f"Dataset in file {file} has no informative characters, therefore it will not be further processed and its data completelly excluded from the output files.")
			non_informative_partitions.append(file)

		else:
			# Parse all data to each species file		
			for name in spp_data:
				spp_data[name].feed(partition)
			part_collection['size'] += partition.metadata['size']
			part_collection['type'] += partition.metadata['type']
			part_collection['states'] += partition.metadata['states']


	# remove uninformative files and spp 
	if not code_indels:
		act_files = [x for x in act_files if not x in non_informative_partitions]
	else:
		to_rm = []
		for file in act_files:
			if file in non_informative_partitions and f"{file}_indels" in non_informative_partitions:
				to_rm.append(file)
		act_files = [x for x in act_files if not x in to_rm]
	
	spp_data = {sp: spp_data[sp] for sp in spp_data if 
				len([x for x in spp_data[sp].metadata["presence"] if x]) > 0}
	
	#for sp in spp_data:
	#	print(sp, spp_data[sp].metadata['type'], spp_data[sp].metadata['presence'])

	if code_gene_content:
		gene_number = 0
		inf_chars = []
		state_counts = {}

		for isp, sp in enumerate(spp_data):
			counter = 0

			with open (spp_data[sp].file, 'a') as fhandle:
				for tipo, present in zip(spp_data[sp].metadata['type'], spp_data[sp].metadata['presence']):
					if tipo == 'nucleic' or tipo == 'peptidic':
						
						if isp == 0:
							state_counts[counter] = {0:0, 1:0}

						if present:
							fhandle.write('1')
							state_counts[counter][1] += 1
						
						else:
							fhandle.write('0')
							state_counts[counter][0] += 1

						if isp == 0:
							gene_number += 1
						
						counter += 1
			
			spp_data[sp].metadata['type'].append('gene_content')
			spp_data[sp].metadata['size'].append(gene_number)			
			#TODO count informative characters in gene content
			spp_data[sp].metadata["informative_chars"].append(inf_chars) # Dummy value
			spp_data[sp].metadata["presence"].append(True)

		for char in state_counts:
			#print(f'{char=} - {state_counts[char]=}')
			if min(list(state_counts[char].values())) > 1: # informative
				inf_chars.append(char)

		part_collection['size'].append(gene_number)
		part_collection['type'].append('gene_content')
		part_collection['states'].append(2)

	#print(f'{part_collection=}')

	#for spe in spp_data:
	#	print(f'{spp_data[spe].metadata=}')

	# Write IQtree phylip files
	iqtree_sets = set(part_collection['type'])

	for settype in iqtree_sets:

		thfile = f'{root_name}_{settype}.phy'
		th_sizes = [part_collection['size'][d] for d in range(len(part_collection['size'])) 
				if part_collection['type'][d] == settype]
		tot_size = sum(th_sizes)
		
		header = f" {len(spp_data)} {tot_size} \n"
		with open(thfile, "a") as oh:
			oh.write(header)
		for sp in spp_data:
			spp_data[sp].parse_phylip_block(thfile, name_space = (longest + 10), 
				partition_type=settype, polymorphs=polys)
			
	

	# Write IQtree nexus file
	with open(iqtree_nexus, 'w') as iqhandle:
		#init = 0
		init = {'nucleic':0, 'peptidic':0, 'indel':0, 'morphological':0, 'gene_content': 0} 
		partinfo = "#nexus\nbegin sets;\n"
		model_spec = "\tcharpartition mine = "

		for ix, thtype in enumerate(part_collection['type']):
			
			if thtype == 'nucleic':
				model_spec += f'GTR+I+G:part{ix+1}, '
				
			elif thtype == 'peptidic':
				model_spec += f'Blosum62:part{ix+1}, '
			
			elif thtype == 'indel':
				model_spec += f'GTR2:part{ix+1}, '

			elif thtype == 'morphological':
				model_spec += f'MK:part{ix+1}, '
							
			elif thtype == 'gene_content':
				model_spec += f'GTR2:part{ix+1}, '
							
			partinfo += f"\tcharset part{ix+1} = {root_name}_{thtype}.phy: {init[thtype]+1}-{init[thtype] + part_collection['size'][ix]};\n"
			init[thtype] += part_collection['size'][ix]

		model_spec = model_spec.rstrip(', ')
		partinfo += model_spec + ';\nend;\n'
		iqhandle.write(partinfo)



	# Write RAxML single phylip matrix
	tot_size = sum(part_collection['size'])
	raxml_header = f" {len(spp_data)} {tot_size} \n"
	with open(raxml_main, "a") as oh:
		oh.write(raxml_header)
	for sp in spp_data:
		spp_data[sp].parse_phylip_block(raxml_main, name_space = (longest + 10), polymorphs=polys)


	# Write RAxML partition file
	with open(raxml_part, 'w') as ph:
		init = 0
		partinfo = ""

		for ix, thtype in enumerate(part_collection['type']):
			
			if thtype == 'nucleic':
				partinfo += 'GTR+I+G, '
				
			elif thtype == 'peptidic':
				partinfo += 'Blosum62, '
			
			elif thtype == 'indel':
				partinfo += 'BIN, '

			elif thtype == 'morphological' or thtype == 'gene_content':
				states = part_collection['states'][ix]
				
				if states == 2:
					partinfo += 'BIN, '
				elif states > 2:
					partinfo += f"MULTI{states}_GTR, "
				else:
					raise ValueError(f"{thtype.capitalize()} partition is uninformative.")
			
			partinfo += f"p{ix+1} = {init+1}-{init + part_collection['size'][ix]}\n"
			init += part_collection['size'][ix]

		ph.write(partinfo)

	# Write xread file

	tnt_header = f"'xread file processed with BAD2matrix.'\n{tot_size} {len(spp_data)}\n"
	
	with open(tnt_main, 'w') as tnt_handle:
		tnt_handle.write(tnt_header)

	for sp in spp_data:
		spp_data[sp].parse_tnt_block(tnt_main, name_space = (longest + 10), polymorphs=polys)

	with open(tnt_main, 'a') as tnt_handle:
		tnt_handle.write(';\n')


	# Remove temporary files
	for name in spp_data:
		spp_data[name].clean()

else:

	print(help_text)


exit()