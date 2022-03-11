import re
from functools import reduce

#dna_codes = 'ABCD  GH  K MN  RST VW Y'  #-> T
#rna_codes = 'ABCD  GH  K MN  RS UVW Y'  #-> U
#aa_codes =  'ABCDEFGHIJKLMNPQRST VWXYZ' #-> EFILPQ JZX

class Partition:

	def __init__(self, filename: str, name_map: dict):

		self.data = {}

		# Describes "subpartitions"
		# Char length, Char type
		# Char types: `nucleic`, `peptidic`, `indel`, `morphological`
		self.metadata=[]

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
		
		if 'peptidic' in types:
			self.metadata = [[list(char_lens.keys())[0], 'peptidic']]
		
		elif 'nucleic' in types:
			self.metadata = [[list(char_lens.keys())[0], 'nucleic']]
		
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

		self.metadata.append([len(thcodes), 'indel'])
	

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


class Term_data:
	"""Simple class for aggregated DNA/AA data of a terminal"""
	
	def __init__(self, name: str):
		self.name = name
		self.file = "temporary_file_for_" + self.name + "_do_not_delete_or_you_will_die.fasta"
		self.partition_table = [] # position start, position end, datatype
		self.size = 0

	def feed(self, part: Partition):

		
		if self.name in part:
			with open(self.file, 'a') as fh:
				fh.write(part[self.name])
			for subpart in part.metadata:
				self.partition_table.append(
					(self.size, self.size + subpart[0], subpart[1])
				)
				self.size += subpart[0]

		else:
			for subpart in part.metadata:
				self.partition_table.append(
					(None, None, subpart[1])
				)


			#with open(self.fasta_file, "wa") as fhandle:
			#	fhandle.write("")
		
		pass



