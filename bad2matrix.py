import sys
import os
import re
import warnings
from functools import reduce
from itertools import combinations
from typing import List, Dict
from datetime import datetime

documentation = {

	'DOI': '[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.10028408.svg)](https://doi.org/10.5281/zenodo.10028408)',

	'Preamble': 'A Python script for merging and translating FASTA alignments into TNT (extended XREAD), RAxML-NG/IQ-Tree (extended PHYLIP) and FastTree (fasta) input matrices. Optionally, it encodes indel characters using the `simple` gap coding method of Simmons and Ochoterena (2000; Gaps as characters in sequence-based phylogenetic analysis. Systematic Biology 49: 369-381. DOI 10.1080/10635159950173889), and gene content as binary characters (absence/presence). This script is slower than 2matrix.pl due to more disk use, but will not run out of RAM (hopefully).',

	'Installation': 'Simply clone the GitHub repository or download the main script (`bad2matrix.py`). A Python 3 interpreter is required.',

	'Usage': {
		'command': 'python bad2matrix.py -d <directory> -n <root-name> [-a 2|3|4|5|6|6dso|6kgb|6sr|8|10|11|12|15|18|20] [-f] [-g] [-i] [-m int] [-t <directory>]',

		'options':{
			'-a': 'Number of amino acid states (default = 20). Reduction with option `6dso` follows [Dayhoff et al. (1978)](http://chagall.med.cornell.edu/BioinfoCourse/PDFs/Lecture2/Dayhoff1978.pdf); option `6kgb` follows [Kosiol et al. (2004)](https://doi.org/10.1016/j.jtbi.2003.12.010); option `6sr` follows [Susko and Roger (2007)](https://doi.org/10.1093/molbev/msm144); option `11` follows [Buchfink et al. (2015)](https://doi.org/10.1038/nmeth.3176); and all other options follow [Murphy et al. (2000)](https://doi.org/10.1093/protein/13.3.149).',
			
			'-d': 'The input directory of aligned FASTA files. The default behavior aggregates sequences of the same species across partitions, in which case names should use the following convention: `>species#sequenceID`. This implies that a species can be represented by only one read within each FASTA file. Characters other than letters, numbers, periods, and underscores will be deleted. Use `-f` for an alternate naming convention.',
			
			'-f': 'Use full FASTA names rather than default settings (see `-d` description for default). Sequences from the same species but different reads will not be aggregated and will be considered distinct OTUs. Characters other than letters, numbers, periods, and underscores will be deleted.',
			
			'-g': 'Do not code gene content (absence/presence). If this flag is not set, gene content is coded.',
			
			'-i': 'Do not code indels. If this flag is not set, indels are coded using simple indel coding following [Simmons and Ochoterena (2000)](https://doi.org/10.1080/10635159950173889).',
			
			'-m': 'Retain the upper `x` percentile of genes in the distribution of missing sequences. By default `x` = 1 (i.e. include all genes with four or more sequences).',
			
			'-n': 'Specify the root-name for output files.',
			
			'-r': 'Folder containing a morphological matrix or a set of ortholog duplication matrices. Datasets should be saved a .tsv tables. Multiple states (polymorphic characters) should be separated by pipes (`\|`).'

			}
	},

	'Input specification': 'The program only accepts alignments in fasta format, and their file extensions should be `fa`, `fas`, or `fasta`. Alignments should also be located in a single folder, which path is parsed with the option `-d`.\n\nOther kinds of data (e.g. a morphological matrix or multiple ortholog encoding tables) can also be concatenated. The first row of each table should contain the character names, and the first column the names of the terminals. These matrices should be saved as tab-separated values: simple text tables with tabs as column separators and `tsv` extension. These tables should be located in a single directory, which is parsed with the option `-r`. Polymorphisms should be separated by pipes (`|`): `leaves_pinnate|leaves_bipinnate`. Be aware that polymorphisms are only fully supported for the TNT output, they will be encoded as missing data in the rest of output datasets.',

	'Sample input/output': 'python bad2matrix.py -d test-data/fastas -f -g -n test',

	'Citation':  'Little, D. P. & N. R. Salinas. 2023. BAD2matrix: better phylogenomic matrix concatenation, indel coding, gene content coding, reduced amino acid alphabets, and occupancy filtering. Software distributed by the authors. DOI: 10.5281/zenodo.10028408.',

	'Additional files': 'Python scripts used for development.\n\n`make_readme.py`: Updates README.md file from source documentation.\n\n`sim_data.py`: Generates fasta alignments, executes BAD2matrix contatenation, and verifies gene content encoding.\n\n`test.py`: Preliminary test functions.',

	'License': '[GPL2](https://github.com/dpl10/BAD2matrix/blob/master/LICENSE)',

	'Related repositories': '[2matrix](https://github.com/nrsalinas/2matrix)',

	'Contact':{
		'Nelson R. Salinas': [ 
        	'nrsalinas@gmail.com', 
    		'New York Botanical Garden'
		],

		'Damon P. Little' : [
			'dlittle@nybg.org', 
			'New York Botanical Garden'
		]
	}	
}

def justify(text, tab_number, text_width):
	out = ''
	tabbies = '\t' * tab_number
	out += tabbies
	thwidth = 0
	breakme = False

	for ichar, char in enumerate(text):
		#print(ichar, char)
		thwidth += 1

		if thwidth % text_width == 0:
			breakme = True

		if breakme and char == ' ':
			out += '\n' + tabbies
			breakme = False
			thwidth = 0

		else:
			out += char

	#for i in range(0, len(text), text_width):
	#	tabbies = '\t' * tab_number
	#	out += f'{tabbies}{text[i:(i+text_width)]}\n'

	return out

def get_cli_help():
	out = '\nbad2matrix.py\n\n'
	preamble_formatted = justify(documentation['Preamble'], 0, 80)
	out += preamble_formatted + '\n\n'
	out += 'USAGE\n\n'
	cl_formatted = justify(documentation['Usage']['command'], 1, 76)
	out += f"{cl_formatted}\n\n"
	out += 'WHERE\n\n'
	for opt, exp in documentation['Usage']['options'].items():
		exp_formatted = justify(exp, 1, 76)
		out += f'{opt}{exp_formatted}\n\n'
	
	citation_formatted = justify(documentation['Citation'], 0, 80)
	out += f"Citation: {citation_formatted}\n\n"
	
	return out


# Handy data dictionaries/lists

valid_fasta_ext = ['fa', 'fasta', 'fan']

nucl_amb_codes = {
	'R': ['A' , 'G'],
	'Y': ['C' , 'T'],
	'S': ['G' , 'C'],
	'W': ['A' , 'T'],
	'K': ['G' , 'T'],
	'M': ['A' , 'C'],
	'B': ['C' , 'G' , 'T'],
	'D': ['A' , 'G' , 'T'],
	'H': ['A' , 'C' , 'T'],
	'V': ['A' , 'C' , 'G'],
	'N': ['A' , 'C' , 'G' , 'T']
}

prot_amb_codes = {
	'B': ['D', 'N'], 
	'J': ['I', 'L'], 
	'X': ['A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L', 'M', 'N', 'P', 
		  'Q', 'R', 'S', 'T', 'V', 'W', 'Y'],
	'Z': ['E', 'Q']
}

nucl2numb = {
	'A': '0',
	'C': '1',
	'G': '2',
	'T': '3',
	'R': '[02]', # ['A' , 'G'],
	'Y': '[13]', # ['C' , 'T'],
	'S': '[12]', # ['G' , 'C'],
	'W': '[03]', # ['A' , 'T'],
	'K': '[23]', # ['G' , 'T'],
	'M': '[01]', # ['A' , 'C'],
	'B': '[123]', # ['C' , 'G' , 'T'],
	'D': '[023]', # ['A' , 'G' , 'T'],
	'H': '[013]', # ['A' , 'C' , 'T'],
	'V': '[012]', # ['A' , 'C' , 'G'],
	'N': '?', # ['A' , 'C' , 'G' , 'T']
}

pep2numb = {
	'A': '0', 
	'C': '1', 
	'D': '2', 
	'E': '3', 
	'F': '4', 
	'G': '5', 
	'H': '6', 
	'I': '7', 
	'K': '8', 
	'L': '9', 
	'M': 'A', 
	'N': 'B', 
	'P': 'C', 
	'Q': 'D', 
	'R': 'E', 
	'S': 'F', 
	'T': 'G', 
	'V': 'H', 
	'W': 'I', 
	'Y': 'J',
	'B': '[2B]', # ['D', 'N'], 
	'J': '[79]', # ['I', 'L'], 
	'Z': '[3D]', # ['E', 'Q']
	'X': '?'
}

def get_name_map(infiles: List[str], full_fasta_names: bool, keep: float = 1.0, 
		 infiles_morph: List[str] = []) -> dict:

	name_map = {}
	file2terms = {}
	
	for file in sorted(infiles + infiles_morph):

		###########################    TODO    #################################
		# File naming convention of inputs should be stated in the instructions.
		# Two classes of input matrices: molecular or morphological/gene dupli-
		# cation events. User should mention which is which through the file
		# extension.
		######################################################################## 
		file_type = None
		thname_map = {}
		file2terms[file] = []
		pattern = '|'.join(valid_fasta_ext)

		if re.search(f'\.({pattern})$', file, re.I):
			file_type = 'fasta'
		elif re.search(r'\.tsv$', file):
			file_type = 'tsv'
		else:
			warnings.warn(f"File `{file}` skipped.")
			continue

		#print(f'{file=}, {file_type=}')
		with open(file, 'r') as fhandle:

			for line_num, line in enumerate(fhandle):
				raw_name = None

				if file_type == 'tsv' and line_num > 0:
					bits = re.split(r'\t', line)
					raw_name = bits[0]
					#print(f'{raw_name=}')

				elif file_type == 'fasta' and line.startswith(">"):
					line = line.lstrip('>')
					raw_name = line.strip()

				if raw_name: # All cleaning procedures of terminal names should be done here.
					name = raw_name
					if not full_fasta_names: name = re.split(r'#+', name)[0]
					name = clean_name(name)
					#print(f'{raw_name=}, {name=}')
					thname_map[raw_name] = name
					file2terms[file].append(name)				

			#print(f'{thname_map=}')
			if len(set(thname_map.values())) < len(thname_map):

				dup_count = {v:0 for v in thname_map.values()}
				for k in thname_map:
					dup_count[thname_map[k]] += 1
				dup_count = {k:v for k,v in dup_count.items() if v > 1}
				err = '\n'.join([k for k in thname_map if thname_map[k] in dup_count])
				raise ValueError(f"Check the following sequence names in ´{file}´, there could be duplicates:\n{err}\n")

			name_map.update(thname_map)

	if keep < 1:
		name_set = set()

		# Filter data files
		file2terms = {z: file2terms[z] for z in	sorted(file2terms, reverse=True, 
			key=lambda x: len(file2terms[x]))}
		new_size = int(keep * len(file2terms))
		file2terms = {x: file2terms[x] for x in list(file2terms.keys())[:new_size]}
		
		# Update taxa names
		for terms in file2terms.values():
			seqs = [x for x,y in name_map.items() if y in terms]
			name_set.update(seqs)

		to_rm = []
		for raw_name in name_map:
			if not raw_name in name_set:
				to_rm.append(raw_name)
		name_map = {x:name_map[x] for x in name_map if not x in to_rm}
	
	#print(f'\n{file2terms=}\n{name_map=}\n')

	return (name_map, list(file2terms.keys()))


def clean_name(name:str) -> str:
	name = re.sub(r'[\s\/\-\\#]+', '_', name) # Verify sharp replacement
	name = re.sub(r'[^\w\._]', '', name)
	return name

#dna_codes = 'ABCD  GH  K MN   RST  VW Y'  #-> T
#rna_codes = 'ABCD  GH  K MN   RS  UVW Y'  #-> U # U can be Selenocisteine!
#aa_codes =  'ABCDEFGHIJKLMNOPQRST UVWXYZ' #-> EFILOPQ JZX


def aa_redux_dict(redux_code: str) -> Dict[str, str]:
	trans_dict = {}
	orig = "OU" # remove pyrrolycine and selenocisteine
	dest = "KC"

	if redux_code == '2':
		orig += "LVIMCAGSTPFYW EDNQKRH"
		dest += "FFFFFFFFFFFFF EEEEEEE"
		trans_dict = str.maketrans(orig, dest)

	elif redux_code == '3':
		orig += "LASGVTIPMC EKRDNQH FYW"
		dest += "IIIIIIIIII EEEEEEE FFF"
		trans_dict = str.maketrans(orig, dest)

	elif redux_code == '4':
		orig += "LVIMC AGSTP FYW EDNQKRH"
		dest += "IIIII PPPPP FFF EEEEEEE"
		trans_dict = str.maketrans(orig, dest)

	elif redux_code == '5':
		orig += "LVIMC ASGTP FYW EDNQ KRH"
		dest += "IIIII PPPPP FFF EEEE KKK"
		trans_dict = str.maketrans(orig, dest)

	elif redux_code == '6':
		orig += "LVIM ASGT PHC FYW EDNQ KR"
		dest += "IIII AAAA PPP FFF EEEE KK"
		trans_dict = str.maketrans(orig, dest)

	elif redux_code == '6dso':
		orig += "AGPST DENQ HKR MIVL WFY C"
		dest += "PPPPP EEEE HHH IIII FFF C"
		trans_dict = str.maketrans(orig, dest)

	elif redux_code == '6kgb':
		orig += "AGPS DENQHKRT MIL W FY CV"
		dest += "PPPP EEEEEEEE III W FF CC"
		trans_dict = str.maketrans(orig, dest)

	elif redux_code == '6sr':
		orig += "APST DENG QKR MIVL WC FYH"
		dest += "PPPP EEEE QQQ IIII CC FFF"
		trans_dict = str.maketrans(orig, dest)

	elif redux_code == '8':
		orig += "LVIMC AG ST P FYW EDNQ KR H"
		dest += "IIIII AA SS P FFF EEEE KK H"
		trans_dict = str.maketrans(orig, dest)

	elif redux_code == '10':
		orig += "LVIM C A G ST P FYW EDNQ KR H"
		dest += "IIII C A G SS P FFF EEEE KK H"
		trans_dict = str.maketrans(orig, dest)
	
	elif redux_code == '11':
		orig += "KREDQN C G H ILV M F Y W P STA"
		dest += "EEEEEE C G H III M F Y W P AAA"
		trans_dict = str.maketrans(orig, dest)

	elif redux_code == '12':
		orig += "LVIM C A G ST P FY W EQ DN KR H"
		dest += "IIII C A G SS P FF W EE DD KK H"
		trans_dict = str.maketrans(orig, dest)

	elif redux_code == '15':
		orig += "LVIM C A G S T P FY W E D N Q KR H"
		dest += "IIII C A G S T P FF W E D N Q KK H"
		trans_dict = str.maketrans(orig, dest)

	elif redux_code == '18':
		orig += "LM VI C A G S T P F Y W E D N Q K R H"
		dest += "LL II C A G S T P F Y W E D N Q K R H"
		trans_dict = str.maketrans(orig, dest)

	orig = orig.replace(" ", "")
	dest = dest.replace(" ", "")
	trans_dict = str.maketrans(orig, dest)

	return trans_dict


def min_steps_char(count_dict, char_type):
	
	min_steps = None
	
	if char_type in ['nucleic', 'peptidic']:
		stand_states = None
		projector = None

		if char_type == 'nucleic':
			stand_states = ['A', 'C', 'G', 'T']
			projector = nucl_amb_codes
		elif char_type == 'peptidic':
			stand_states = ['A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L', 
				'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y']
			projector = prot_amb_codes

		proj_states = {x for x in count_dict if x in stand_states}
		ambs = {x for x in count_dict if not x in stand_states}
		
		# Check if ambiguities are already represented in projected set
		to_rm = []
		for amb in ambs:
			thset = set(projector[amb])
			inter = proj_states & thset
			if len(inter) > 0:
				to_rm.append(amb)
		ambs = {x for x in ambs if not x in to_rm}

		# Check if two ambiguities have common aa/nucleotide
		#? What if there is only one amb left?
		l0 = 0
		l1 = 1
		while len(ambs) > 1 and l0 != l1:
			l0 = len(ambs)
			to_rm = []
			for i,d in combinations(ambs, 2):
				iset = set(projector[i])
				dset = set(projector[d])
				inter = iset & dset
				if len(inter) > 0:
					proj_states.update(next(iter(inter)))
					to_rm += [i, d]
					break
			ambs = {x for x in ambs if not x in to_rm}
			l1 = len(ambs)

		# Just add the non-empathetic ambiguities
		for amb in ambs:
			proj_states.update(next(iter(amb)))

		min_steps = len(proj_states) - 1

	else:
		min_steps = len(count_dict) - 1

	return min_steps


def max_steps_char(count_dict, char_type):

	max_steps = None

	if char_type in ['nucleic', 'peptidic']:
		stand_states = None
		projector = None

		if char_type == 'nucleic':
			stand_states = ['A', 'C', 'G', 'T']
			projector = nucl_amb_codes
		elif char_type == 'peptidic':
			stand_states = ['A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L', 
				'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y']
			projector = prot_amb_codes

		count_dict = {x: count_dict[x] for x in sorted(count_dict, reverse=True, key=lambda y: count_dict[y])}
		new_count = {x: count_dict[x] for x in count_dict if x in stand_states} #stand_states?
		ambs = {x: count_dict[x] for x in count_dict if not x in stand_states}

		to_rm = [0]
		while len(to_rm) > 0:
			to_rm = []
			break_ext = False
			for amb in ambs:
				for sym in new_count:
					if sym in projector[amb]:
						new_count[sym] += ambs[amb]
						to_rm.append(amb)
						break_ext = True
						break
				if break_ext: break	
			ambs = {x: ambs[x] for x in ambs if not x in to_rm}

		if len(ambs) > 0:
			new_count.update(ambs)
			new_count = {x: new_count[x] for x in sorted(new_count, reverse=True, key=lambda y: new_count[y])}

		max_steps = sum(list(new_count.values())[1:])

	else:
		count_dict = {x: count_dict[x] for x in sorted(count_dict, reverse=True, key=lambda y: count_dict[y])}

		max_steps = sum(list(count_dict.values())[1:])

	return max_steps


class Polymorphs:

	def __init__(self):
		"""
		Simple map to keep record of TNT polymorphic encodings. 
		"""
		self.counter = 130
		self.mapping = {}

	def add_poly_encoding(self, poly_str: str):
		"""
		Adds a new TNT polymorphic encoding into the morphological partition.
		First and last characters should be square brackets. Returns a single 
		unicode character that will be inserted in the phylip-like, temporary 
		storing file that the script handles in the background.
		"""
		thchar = chr(self.counter)
		self.counter += 1
		self.mapping[thchar] = poly_str
		return thchar


class Partition:

#TODO###########################################################################
#
#		- Make state_translations an attribute of Partition class, but save
#		the data if it is actual morphological data, not if it is ortholog
#		duplication encodings.
#
#		- Create a method in Term_data class to retrieve state_translations
#		from a Partition and write the tnt char set names block
#
#TODO###########################################################################


	def __init__(self, filename: str, name_map: dict, translation_dict: dict=None,
		polymorphs: Polymorphs=None):

		self.data = {}
		self.filetype = None
		self.origin = filename
		
		self.metadata = { # Describes "subpartitions"
			"size": [], # Char length
			"type": [], # Char types: `nucleic`, `peptidic`, `indel`, `morphological`, `gene_content`
			"informative_chars": [], # Number of Informative positions
			"origin": [filename], 
			"character_names": [],
			"states" : []
			}
		
		#TODO######   Include name in metadata   ###########

		pattern = '|'.join(valid_fasta_ext)
		if re.search(f'\.({pattern})$', self.origin, re.I):
			self.filetype = 'fasta'

		elif re.search(r'\.tsv$', self.origin):
			self.filetype = 'tsv'

		with open(self.origin, 'r') as fhandle:
			#print(f'{self.origin=}')
			char_lens = {}
			th_term = ''
			th_seq = ''
			types = {}

			if self.filetype == 'fasta':
				
				for line in fhandle:
					line = line.strip()
					#print(f'{line=}')
					if line.startswith('>'):

						if th_term and th_seq:
							char_lens[len(th_seq)] = 0

							if len(char_lens.keys()) > 1:
								raise ValueError(f"Sequences in {filename} have different lengths, probably they are not aligned.")

							th_seq = th_seq.upper()
							thtype = self.seq_type(th_seq)
							types[thtype] = 0
							self.data[th_term] = th_seq
							th_term = ''
							th_seq = ''

						th_term = name_map[line.lstrip('>')]

					else:
						th_seq += line.upper()

				# Capture data of last fasta entry
				if th_term and th_seq:
					char_lens[len(th_seq)]

					if len(char_lens.keys()) > 1:
						raise ValueError(f"Sequences in {filename} have different lengths, probably they are not aligned.")

					th_seq = th_seq.upper()
					thtype = self.seq_type(th_seq)
					types[thtype] = 0
					self.data[th_term] = th_seq

				self.metadata['states'].append(None)

			elif self.filetype == 'tsv':

				types['morphological'] = 0
				charset = set()
				state_translations = {} # { character : { original state : new state } }
				max_state = [] # count of states per character (starts at zero)

				for line_idx, line in enumerate(fhandle):
					line = line.strip()

					if line_idx == 0:
						self.metadata["character_names"] = re.split(r'\t', line)[1:]

					else:
						bits = re.split(r'\t', line)
						th_seq = ''
						char_lens[len(bits[1:])] = 0
						if len(char_lens.keys()) > 1:
							raise ValueError(f"OTUs in {filename} have different character observations, check for missing data.")

						#######################  polymorphsm block  ########################## 

						for ichar, char in enumerate(bits[1:]):
							if not ichar in state_translations:
								state_translations[ichar] = {}
								max_state.append(0)
							tidbits = re.split(r'\|', char)
							#print(f'{tidbits=}')
							if len(tidbits) > 1: #polymorphic
								#print("It is polymorphic")
								warnings.warn(f"File `{self.origin}` contains polymorphisms. They will be encoded as missing data in the phylip matrix.")
								poly = '['
								for sta in tidbits:
									if not sta in state_translations[ichar]:
										state_translations[ichar][sta] = str(max_state[ichar])
										max_state[ichar] += 1
									poly += state_translations[ichar][sta]
								poly += ']'
								single_char = polymorphs.add_poly_encoding(poly)
								th_seq += single_char
							else:
								#print("It is not polymorphic")
								#print(f'{char=}')
								if '?' != char:
									if not char in state_translations[ichar]:
										state_translations[ichar][char] = str(max_state[ichar])
										max_state[ichar] += 1
									#print(f'{state_translations=}')
									th_seq += state_translations[ichar][char]
								else:
									th_seq += '?'

						#######################  polymorphsm block  ########################## 

						th_term = name_map[bits[0]]
						#th_seq = ''.join(bits[1:])
						charset.update(set(th_seq))
						self.data[th_term] = th_seq

				#Checking proper state conventions
				if '?' in charset:
					charset.remove('?')
								
				self.metadata['states'].append(len(charset))

		#print(f'{char_lens=}')

		types = list(set(types.keys()))
		self.metadata["size"].append(list(char_lens.keys())[0])
		self.metadata["informative_chars"].append([])

		if 'peptidic' in types:
			self.metadata["type"].append("peptidic")
		
		elif 'nucleic' in types:
			self.metadata["type"].append("nucleic")
		
		elif 'morphological' in types:
			self.metadata["type"].append("morphological")

		else:

			#TODO ########    Raise the appropriate error here    ##############

			raise ValueError('WTF')

		# Peptidic state reduction
		if translation_dict and self.metadata["type"][-1] == "peptidic":
			for term in self.data:
				self.data[term] = self.data[term].translate(translation_dict)



	def indel_coder(self):
		indel_map = {name:{} for name in self.data}

		# find indel morphology
		for taxon in self.data:
			valid_init = 0
			valid_end = 0
			in_gap = False
			thgap = [None, None]
			thindels = {}

			for ichar, char in enumerate(self.data[taxon]):
			
				if char != '-':
					valid_init = ichar
					break
			
			for ichar, char in reversed(list(enumerate(self.data[taxon]))):
			
				if char != '-':
					valid_end = ichar
					break
			
			for ichar in range(valid_init, valid_end+1):
			
				if self.data[taxon][ichar] == '-' and not in_gap:
					thgap[0] = ichar
					in_gap = True
			
				elif self.data[taxon][ichar] != '-' and in_gap:
					thgap[1] = ichar
					thindels[tuple(thgap)] = 0
					thgap = [None, None]
					in_gap = False
			
			indel_map[taxon] = thindels
		
		#print('\nindel_map', indel_map)

		# count character frequency
		indel_count = {}
		
		for taxon in indel_map:
		
			for indel in indel_map[taxon]:
		
				if indel in indel_count:
					indel_count[indel] += 1
		
				else:
					indel_count[indel] = 1
		
		#print('\nindel_count', indel_count)

		# find ambiguity dependency among indels
		indel_set = sorted([i for i in indel_count if indel_count[i] > 1])
		ambiguity_map = {x:[] for x in indel_set}
		
		for indel in indel_set:
		
			for indel_ in indel_set:
		
				if indel != indel_:
		
					if indel[0] <= indel_[0] and indel[1] >= indel_[1]:
						ambiguity_map[indel].append(indel_)
		
		#print('\nambiguity_map', ambiguity_map)

		# remove non-informative
		for indel in indel_count:
		
			if indel_count[indel] == 1:
		
				for taxon in indel_map:
		
					if indel in indel_map[taxon]:
						indel_map[taxon].pop(indel)
		
		#print('\nindel_map', indel_map)

		# code characters
		for taxon in self.data:
			thcodes = {x:None for x in indel_set}
			ambiguous = []
		
			for indel in indel_set:
		
				if indel in indel_map[taxon]:
					thcodes[indel] = '1'
		
					if len(ambiguity_map[indel]) > 0:
						ambiguous.append(indel)
		
				else:
					thcodes[indel] = '0'

			for amb in ambiguous:
				thcodes[amb] = '?'

			#print('\nthcodes', thcodes)
			self.data[taxon] += ''.join(list(thcodes.values()))

		self.metadata["size"].append(len(indel_set))
		self.metadata["type"].append("indel")	
		self.metadata["informative_chars"].append([])
		self.metadata['states'].append(None)
		self.metadata["origin"].append(self.origin)


	def seq_type(self, sequence):

		seq_type = None
		no_valid = re.compile(r'[^ABCDEFGHIJKLMNOPQRSTUVWXYZ\-]')
		prot = re.compile(r'[EFILOPQJZX]')
		robj_no_valid = no_valid.search(sequence)
		
		if robj_no_valid:
			raise ValueError(f"Sequence contains a symbol not valid for molecular partition: `{robj_no_valid.group(0)}`.")
		
		elif prot.search(sequence):
			seq_type = 'peptidic'
		
		else:
			seq_type = 'nucleic'
		
		return seq_type


	def informative_stats(self):
		acc = 0
		#print(f'{self.metadata["size"]=}')
		for sub_idx, sub_size in enumerate(self.metadata["size"]):
			#print(f'{sub_idx=}, {sub_size=}')
			for idx in range(acc, (acc + sub_size)):
				states = {}
				
				for term in self.data:
					pos = self.data[term][idx]
					if not pos in ['-', '?']:
						if pos in states:
							states[pos] += 1
						else:
							states[pos] = 1

				#print(f'{idx=}, {states=}, {self.metadata["type"][sub_idx]}')
				if len(states) > 1:

					min_steps = min_steps_char(states, self.metadata["type"][sub_idx])
					max_steps = max_steps_char(states, self.metadata["type"][sub_idx])
					#print(f"{min_steps=}, {max_steps=}")
					
					if max_steps > min_steps:
						#print('idx-acc:', (idx-acc))
						self.metadata["informative_chars"][sub_idx].append(idx-acc)
						#print(self.metadata)
			
			acc += sub_size

		return None


def get_informative_stats(part_sizes, part_types, part_info_chars, data_matrix):
	acc = 0

	for sub_idx, sub_size in enumerate(part_sizes):
		for idx in range(acc, (acc + sub_size)):
			states = {}
			
			for term in data_matrix:
				pos = data_matrix[term][idx]
				if not pos in ['-', '?']:
					if pos in states:
						states[pos] += 1
					else:
						states[pos] = 1

			if len(states) > 1:
				min_steps = min_steps_char(states, part_types[sub_idx])
				max_steps = max_steps_char(states, part_types[sub_idx])
				
				if max_steps > min_steps:
					part_info_chars[sub_idx].append(idx-acc)
		
		acc += sub_size

	return None


class Term_data:
	"""Simple class for aggregated DNA/AA data of a terminal"""
	
	def __init__(self, name: str, gene_encoding: bool):
		self.name = name
		self.file = "temporary_file_for_" + self.name + "_do_not_delete_or_you_will_die.txt"
		self.metadata = {"size": [], "type": [], "informative_chars": [], "presence": []} #, "origin": []}
		self.size = 0
		self.gene_encoding = gene_encoding


	def feed(self, part: Partition):
		
		tot_inf = len(reduce(lambda x, y: x + y, part.metadata["informative_chars"]))

		if tot_inf > 0:

			is_present = False
			
			if self.name in part.data:
				is_present = True # All subpartitions tagged as present
				
				with open(self.file, 'a') as fh:
					fh.write(part.data[self.name])

			#TODO Most of the metadata (except presence) should be held globally, it is the same across partitions
			self.metadata["size"] += part.metadata["size"]
			self.metadata["type"] += part.metadata["type"]
			self.metadata["informative_chars"] += part.metadata["informative_chars"]
			#self.metadata["origin"] += part.metadata["origin"]
			self.metadata["presence"] += [is_present for x in part.metadata["size"]]

		return None


	def clean(self):
		if os.path.exists(self.file):
			os.remove(self.file)
	

	def parse_tnt_block(self, outfile: str, name_space: int = 20, polymorphs: Polymorphs = None):
		#print(self.name)
		with open(outfile, 'a') as outhandle:
			pad = name_space - len(self.name)
			outhandle.write(self.name + " " * pad)

			with open(self.file, 'r') as inhandle:		

				for data in inhandle:
					init = 0
					end = 0
				
					for ipart, isize in enumerate(self.metadata["size"]):

						if self.metadata["presence"][ipart]: 
							#print(self.metadata["informative_chars"][ipart])
							init = end
							end = init + isize
							tmp = data[init:end]
							#print(f'{len(tmp)=}')
							tmp = [x for i,x in enumerate(list(tmp)) if i in 
			  					self.metadata["informative_chars"][ipart]]
							tmp = ''.join(tmp)
							#print(f'{len(tmp)=}')
							transdict = None

							if self.metadata['type'][ipart] == 'nucleic':
								transdict = nucl2numb

							elif self.metadata['type'][ipart] == 'peptidic':
								transdict = pep2numb
							
							if transdict: # Translate aaa or nucleic seq to TNT numbers
								for mol in transdict:
									tmp = re.sub(mol, transdict[mol], tmp)
							
							if not polymorphs is None:
								for poly_symbol in polymorphs.mapping:
									tmp = re.sub(poly_symbol, polymorphs.mapping[poly_symbol], tmp)

							tmp = re.sub(r'\-', '?', tmp) # Just to have all missing data as '?' 
							outhandle.write(tmp)

						else:
							#outhandle.write("?" * self.metadata["size"][ipart])
							outhandle.write("?" * len(self.metadata["informative_chars"][ipart]))

			outhandle.write('\n')


	def parse_phylip_block(self, outfile: str, name_space: int = 20, 
		partition_type: str = 'all', polymorphs: Polymorphs = None):

		if partition_type == 'all':
			partition_type = ['nucleic', 'peptidic','indel', 'morphological', 'gene_content']
		else:
			partition_type = [partition_type]

		with open(outfile, 'a') as ohandle:
			pad = name_space - len(self.name)
			ohandle.write(self.name + " " * pad)

			with open(self.file, 'r') as ihandle:
					
				for data in ihandle:
					init = 0
					end = 0
				
					for ipart, isize in enumerate(self.metadata["size"]):

						if self.metadata["presence"][ipart]: 
							init = end
							end = init + isize

							if self.metadata["type"][ipart] in partition_type:
								tmp = data[init:end]
								for poly_symbol in polymorphs.mapping:
									tmp = re.sub(poly_symbol, '?', tmp)
								#ohandle.write(re.sub(r'\[\w+\]', '?', data[init:end]))
								#ohandle.write(data[init:end])
								ohandle.write(tmp)

						else:
							if self.metadata["type"][ipart] in partition_type: # write missing data
								ohandle.write("-" * self.metadata["size"][ipart])
			"""
			if self.gene_encoding:

				for tipo, present in zip(self.metadata['type'], self.metadata['presence']):
				
					if tipo == 'nucleic' or tipo == 'peptidic':
				
						if present:
							ohandle.write('X')
						else:
							ohandle.write('|')
			"""
			ohandle.write('\n')



	def parse_fasta_block(self, outfile: str, polymorphs: Polymorphs = None):

		partition_type = ['nucleic', 'peptidic']

		with open(outfile, 'a') as ohandle:
			ohandle.write(f'>{self.name}\n')

			with open(self.file, 'r') as ihandle:
					
				for data in ihandle:
					init = 0
					end = 0
				
					for ipart, isize in enumerate(self.metadata["size"]):

						if self.metadata["presence"][ipart]: 
							init = end
							end = init + isize

							if self.metadata["type"][ipart] in partition_type:
								tmp = data[init:end]
								for poly_symbol in polymorphs.mapping:
									tmp = re.sub(poly_symbol, '-', tmp)
								ohandle.write(tmp)

						else:
							if self.metadata["type"][ipart] in partition_type: # write missing data
								ohandle.write("-" * self.metadata["size"][ipart])

			ohandle.write('\n')




if __name__ == '__main__':

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
	tsv_file = None
	raxml_bffr = ""
	log_bffr = "\n\nBAD2matrix execution log\n\n"
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
				if file.endswith('.tsv'):
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

		ya = datetime.now()
		ya = ya.strftime('%B %d, %Y - %H:%M:%S')
		log_bffr += f'{ya}\n\nCommand arguments: '
		log_bffr += ' '.join(sys.argv) + '\n\n'


		translation_dict = aa_redux_dict(aa_encoding)
		fasttree_main = os.path.join('fasttree_datasets', f'{root_name}.fasta')
		raxml_main = os.path.join('raxml_datasets', f'{root_name}.phy')
		raxml_part = os.path.join('raxml_datasets', f'{root_name}.part')
		iqtree_nexus = os.path.join('iqtree_datasets', f'{root_name}.nex')
		tnt_main = os.path.join('tnt_datasets', f'{root_name}.ss')
		(name_map, act_files) = get_name_map(infiles, full_fasta_names, keep_percentile, infiles_morph)
		term_names = sorted(list(set(name_map.values()))) #? Why sort should be done in reverse order?
		longest = len(max(term_names, key = len))
		spp_data = {name: Term_data(name, code_gene_content) for name in term_names}
		non_informative_partitions = []
		final_spp_count = 0
		part_collection = {'size': [], 'type': [], 'states': [], 'informative_chars': [], 'file': []}
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
				part_collection['informative_chars'] += partition.metadata['informative_chars']
				part_collection['file'] += partition.metadata['origin']


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
				spp_data[sp].metadata["informative_chars"].append(inf_chars) # Dummy value
				spp_data[sp].metadata["presence"].append(True)

			for char in state_counts:
				#print(f'{char=} - {state_counts[char]=}')
				if min(list(state_counts[char].values())) > 1: # informative
					inf_chars.append(char)

			part_collection['size'].append(gene_number)
			part_collection['type'].append('gene_content')
			part_collection['states'].append(2)
			part_collection['informative_chars'].append(inf_chars)

		#for spe in spp_data:
		#	print(f'{spp_data[spe].metadata=}')

		#print(f'{part_collection=}')

		if not os.path.exists('iqtree_datasets'):
			os.mkdir('iqtree_datasets')

		# Write IQtree phylip files
		iqtree_sets = set(part_collection['type'])

		for settype in iqtree_sets:

			thfile = f'iqtree_datasets/{root_name}_{settype}.phy'
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


		# Write FastTree fasta matrix
		if not os.path.exists('fasttree_datasets'):
			os.mkdir('fasttree_datasets')

		for sp in spp_data:
			spp_data[sp].parse_fasta_block(fasttree_main, polymorphs=polys)


		# Write RAxML single phylip matrix

		if not os.path.exists('raxml_datasets'):
			os.mkdir('raxml_datasets')

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

		if not os.path.exists('tnt_datasets'):
			os.mkdir('tnt_datasets')

		tot_size = sum([len(k) for k in part_collection['informative_chars']])

		tnt_header = f"xread\n'File processed with BAD2matrix.'\n{tot_size} {len(spp_data)}\n"
		
		with open(tnt_main, 'w') as tnt_handle:
			tnt_handle.write(tnt_header)

		for sp in spp_data:
			spp_data[sp].parse_tnt_block(tnt_main, name_space = (longest + 10), polymorphs=polys)

		with open(tnt_main, 'a') as tnt_handle:
			tnt_handle.write(';\n')


		# Remove temporary files
		for name in spp_data:
			spp_data[name].clean()


		# Body of log file
		log_bffr += 'Partitions processed:\n\n'

		for idx , filename in enumerate(part_collection['file']):
			log_bffr += f'{idx+1}: {filename}, {part_collection["type"][idx]} ({part_collection["size"][idx]} character{"s" if part_collection["size"][idx] > 1 else ""}).\n'		

		print(log_bffr)

	else:
		help_text = get_cli_help()
		print(help_text)


	exit()
