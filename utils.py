import re
from functools import reduce
import warnings

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

def get_name_map(infiles: list, full_fasta_names: bool) -> dict:

	name_map = {}
	
	for file in infiles:

		########################################################################
		# File naming convention of inputs should be stated in the instructions.
		# Two classes of input matrices: molecular or morphological/gene dupli-
		# cation events. User should mention which is which through the file
		# extension.
		######################################################################## 
		if re.search(r'\.(fas|fasta)$', file):

			with open(file, 'r') as fhandle:
				thname_map = {}

				for line in fhandle:

					if line.startswith(">"):
						# All cleaning procedures of terminal names should be done here.
						raw_name = line.strip()
						name = raw_name
						if not full_fasta_names:
							name = re.split(r'#+', name)[0]
						name = clean_name(name)
						thname_map[raw_name] = name

				if len(set(thname_map.values())) < len(thname_map):

					thnames = list(thname_map.keys())
					dup_count = list(map(lambda x: thnames.count(x), thnames))
					prob = [x for x,y in zip(thnames, dup_count) if y > 1]
					err = '\n'.join(prob)
					raise ValueError(f"Check the following sequence names in ´{file}´, there could be duplicates:\n{err}.\n")

				name_map.update(thname_map)

		else:
			warnings.warn(f"File `{file}` skipped.")

	return name_map

def clean_name(name:str) -> str:
	name = re.sub(r'[\s\/\-\\#]+', '_', name) # Verify sharp replacement
	name = re.sub(r'[^\w\._]', '', name)
	return name

#dna_codes = 'ABCD  GH  K MN  RST VW Y'  #-> T
#rna_codes = 'ABCD  GH  K MN  RS UVW Y'  #-> U
#aa_codes =  'ABCDEFGHIJKLMNPQRST VWXYZ' #-> EFILPQ JZX

class Partition:

	def __init__(self, filename: str, name_map: dict):

		self.data = {}

		# Describes "subpartitions"
		# Char length, Char type, Informative count, Informative positions
		# Char types: `nucleic`, `peptidic`, `indel`, `morphological`
		self.metadata = {"size": [], "type": [], "informative_chars": []}
		#==> Include name in metadata <==

		with open(filename, 'r') as fhandle:
			char_lens = {}
			th_term = ''
			th_seq = ''
			types = {}

			for line in fhandle:
				line = line.strip()

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

					th_term = name_map[line]

				else:
					th_seq += line.upper()

			if th_term and th_seq:
				char_lens[len(th_seq)]

				if len(char_lens.keys()) > 1:
					raise ValueError(f"Sequences in {filename} have different lengths, probably they are not aligned.")

				th_seq = th_seq.upper()
				thtype = self.seq_type(th_seq)
				types[thtype] = 0
				self.data[th_term] = th_seq

		types = list(set(types.keys()))
		self.metadata["size"].append(list(char_lens.keys())[0])
		self.metadata["informative_chars"].append([])

		if 'peptidic' in types:
			self.metadata["type"].append("peptidic")
		
		elif 'nucleic' in types:
			self.metadata["type"].append("nucleid")
		
		else:
			raise ValueError('WTF')


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


	def seq_type(self, sequence):

		seq_type = None
		no_valid = re.compile(r'[^ABCDEFGHIJKLMNPQRSTUVWXYZ\-]')
		prot = re.compile(r'[EFILPQJZX]')
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
			print(f'{sub_idx=}, {sub_size=}')
			for idx in range(acc, (acc + sub_size)):
				states = {}
				
				for term in self.data:
					seq = self.data[term]
					if not seq[idx] in ['-', '?']:
						if seq[idx] in states:
							states[seq[idx]] += 1
						else:
							states[seq[idx]] = 1
				print(f"{idx=}, {states=}")
				if len(states) > 1:
					min_steps = len(states) - 1
					max_steps = 0
					sorted_states = sorted(states, key = lambda x : states[x])
					for st in sorted_states[:-1]:
						max_steps += states[st]
					print(f"{min_steps=}, {max_steps=}, {sorted_states=}")
					
					if max_steps > min_steps:
						self.metadata["informative_chars"][sub_idx].append(idx-acc)
			#print(f"{sub_idx=}")
			#print(f"{self.metadata['informative_chars'][sub_idx]=}")
			not_inf = [x for x in range(sub_size) if not x in self.metadata['informative_chars'][sub_idx]]
			#print(f"{not_inf=}\n")
			acc += sub_size
		

class Term_data:
	"""Simple class for aggregated DNA/AA data of a terminal"""
	
	def __init__(self, name: str):
		self.name = name
		self.file = "temporary_file_for_" + self.name + "_do_not_delete_or_you_will_die.fasta"
		self.partition_table = [] # position size, datatype, informative_chars
		self.size = 0

	def feed(self, part: Partition):
		
		present = False
		
		if self.name in part.data:
			present = True
			with open(self.file, 'a') as fh:
				fh.write(part.data[self.name])

		for idx in range(len(part.metadata["size"])):
			self.partition_table.append(
				(part.metadata["size"][idx], 
				part.metadata["type"][idx],
				part.metadata["informative_chars"][idx],
				present
				)
			)


def min_steps_char(count_dict, char_type):
	""""
	Fails if the character consists of two different ambiguity states and there
	is at least one shared possible nucleotide shared in their encoding.
	"""

	min_steps = None
	
	if char_type in ['nucleic', 'peptidic']:
		stand_states = None
		projector = None

		if char_type == 'nucleic':
			stand_states = ['A', 'C', 'G', 'T']
			projector = nucl_amb_codes
		else:
			stand_states = ['A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L', 
				'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y']
			projector = prot_amb_codes

		proj_states = {x for x in count_dict if x in stand_states}
		ambs = {x for x in count_dict if not x in stand_states}
		
		if len(ambs) > 0:
			for amb in ambs:
				thset = set(projector[amb])
				inter = proj_states & thset
				if len(inter) > 0:
					break
				else:
					proj_states.update(next(iter(thset))) #==> change to for loop <==

		min_steps = len(proj_states)

	else:
		min_steps = len(count_dict)

	return min_steps



