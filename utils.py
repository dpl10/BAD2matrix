import re

def get_partition(filename: str, name_map: dict) -> dict:

	partition = {name: "" for name in name_map.values()}

	with open(filename, 'r') as fhandle:
		char_lens = {}
		th_term = ''
		th_seq = ''

		for line in fhandle:
			line = line.strip()

			if line.startswith('>'):

				if th_term and th_seq:
					char_lens[len(th_seq)]

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

	partition = {n: partition[n] for n in partition if len(partition[n]) > 0}

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



