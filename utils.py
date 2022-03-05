import re
from functools import reduce

def get_partition(filename: str, name_map: dict) -> dict:

	#partition = {name: "" for name in name_map.values()}
	partition = {}

	with open(filename, 'r') as fhandle:
		char_lens = {}
		th_term = ''
		th_seq = ''

		for line in fhandle:
			line = line.strip()

			if line.startswith('>'):

				if th_term and th_seq:
					char_lens[len(th_seq)] = 0

					if len(char_lens.keys()) > 0:
						raise ValueError(f"Sequences in {filename} have different lengths, probably they are not aligned.")

					partition[th_term] = th_seq
					th_term = ''
					th_seq = ''

				th_term = name_map[line]

			else:
				th_seq += line

		if th_term and th_seq:
			char_lens[len(th_seq)]

			if len(char_lens.keys()) > 0:
				raise ValueError(f"Sequences in {filename} have different lengths, probably they are not aligned.")

			partition[th_term] = th_seq

	return partition


def code_indels(partition:dict) -> dict:
	indel_map = {name:{} for name in partition}

	# find indel morphology
	for taxon in partition:
		valid_init = 0
		valid_end = 0
		in_gap = False
		thgap = [None, None]
		thindels = {}
		for ichar, char in enumerate(partition[taxon]):
			if char != '-':
				valid_init = ichar
				break
		for ichar, char in reversed(list(enumerate(partition[taxon]))):
			if char != '-':
				valid_end = ichar
				break
		for ichar in range(valid_init, valid_end+1):
			if partition[taxon][ichar] == '-' and not in_gap:
				thgap[0] = ichar
			elif partition[taxon][ichar] != '-' and in_gap:
				thgap[1] = ichar-1
				thindels[tuple(thgap)] = 0
				thgap = [None, None]
		indel_map[taxon] = thindels

	# count character frequency
	indel_count = {}
	for taxon in indel_map:
		for indel in indel_map[taxon]:
			if indel in indel_count:
				indel_count[indel] += 1
			else:
				indel_count[indel] = 1

	# find ambiguity dependency among indels
	indel_set = sorted([i for i in indel_count if indel_count[i] > 1])
	ambiguity_map = {x:[] for x in indel_set}
	for indel in indel_set:
		for indel_ in indel_set:
			if indel != indel_:
				if indel[0] <= indel_[0] and indel[1] >= indel_[1]:
					ambiguity_map[indel].append(indel_)

	# remove non-informative
	for indel in indel_count:
		if indel_count[indel] == 1:
			for taxon in indel_map:
				indel_map[taxon].pop(indel)
	
	# code characters
	for taxon in partition:
		thcodes = [None for x in indel_set]
		ambiguous_idx = []
		for idx, indel in enumerate(indel_set):
			if indel in indel_map[taxon]:
				thcodes.append('1')
				if len(ambiguity_map[indel]) > 0:
					ambiguous_idx.append(idx)
			else:
				thcodes.append('0')
		thcodes = [x if not x in ambiguous_idx else '?' for i,x in enumerate(thcodes)]
		partition[taxon].append(''.join(thcodes))

	return partition


class Term_data:
	"""Simple class for aggregated DNA/AA data of a terminal"""
	
	def __init__(self, name: str):
		self.name = name
		self.fasta_file = "temporary_file_for_" + self.name + "_do_not_delete_or_you_will_die.fasta"
		self.partition_table = [] # pos init, pos end, datatype
		
		with open(self.fasta_file, "w") as fhandle:
			fhandle.write("")

	def insert(self, data, seqtype):
		if seqtype == 'molecular':
			data = re.sub(r'[ABCDEFGHIJKLMNOPQRSTUVWXYZ]', '' , data)
			
		elif seqtype == 'morpho':
			#
			# Implement morpho/indel coding
			#
			pass

		with open(self.fasta_file, "wa") as fhandle:
			fhandle.write("")
		
		pass



