import os
from bad2matrix import get_name_map, clean_name, Partition, Term_data

infiles = ['alg_test_0.fasta', 'alg_test_1.fasta', 'alg_test_2.fasta']

dummy = ['''
>sp0#sample0
TATTCCTCTATTA---GTAATTGGGCTTCTACTT-TTCCAAGAGCAACTAAAAATATTCGGCGTATC---
>sp1#sample0
-ATTCCTCTATTA---GTAATTGGGCTTCTACTT-TTCCATGAGCAACTAAAAATATACGGCGTATCT--
>sp2#sample0
--TTCCTCTATTAA-AGGAATTGGG--TCTACATTTTCCACGAGCAACTACCCATATTCGGCGTATCTG-
>sp3#sample0
---TCCTCTATTAA-AGGAATTGGG--TCTACATTTTCCAAGAGCAACTACCCATATTC-GCGTATCTGG
''',
'''
>sp0#sample0
CTTCTTTGCATTTATTACGATCGATTCTCCATGAATG------TAGTTTTAGTAAAGAAAATTTGCAGAAATCTCTGATT
>sp1#sample1
CTTCGTTGCATTTATTACGATC-ATTCTCCATGAATG------TAGTTTTAGTAAAGAATATTTGCAGAAATCTCTGACT
>sp2#sample0
CTTCGCTGCAT--ATTACGATC-ATTCTCCATGAATG------TAGTTTTAGTTAAGAATATTTGCAGAAATCTCTGACT
>sp4#sample0
CTTCGCTGCAT--ATTACGATC-ATTCTCCATGAATGATAATCTAGTTTTAGTTAAGAATATTTGCAGAAATCTCTGATT
''',
'''
>sp0#sample0
TATTCCTCTATTA---GTAATTGGGCTTCTACTT-TTCCAAGAGGGCGGAAAAATATTCGGCGTATC---
>sp5#sample0
TATTCCTCTATTA---GTAATTGGGCTTCTACTT-TTCCAAGAGCGGGGAAAAATATTCGGCGTATC---
'''
]

for inf, du in zip(infiles, dummy):
	with open(inf, 'w') as fh:
		fh.write(du)


def test_dummy_files():
	assert os.path.exists(infiles[0]) == True
	assert os.path.exists(infiles[1]) == True

def test_name_map():
	mn, tc = utils.get_name_map(infiles, True)
	assert mn == {
		'>sp0#sample0': 'sp0_sample0',
		'>sp1#sample0': 'sp1_sample0',
		'>sp2#sample0': 'sp2_sample0',
		'>sp3#sample0': 'sp3_sample0',
		'>sp1#sample1': 'sp1_sample1',
		'>sp4#sample0': 'sp4_sample0',
		'>sp5#sample0': 'sp5_sample0'
	}

	assert tc == ['alg_test_0.fasta', 'alg_test_1.fasta', 'alg_test_2.fasta']

	mn, tc = get_name_map(infiles, False)
	assert mn == {
		'>sp0#sample0': 'sp0',
		'>sp1#sample0': 'sp1',
		'>sp2#sample0': 'sp2',
		'>sp3#sample0': 'sp3',
		'>sp1#sample1': 'sp1',
		'>sp4#sample0': 'sp4',
		'>sp5#sample0': 'sp5'
	}
	assert tc == ['alg_test_0.fasta', 'alg_test_1.fasta', 'alg_test_2.fasta']

	mn, tc = get_name_map(infiles, False, 0.8)
	assert mn == {
		'>sp0#sample0': 'sp0',
		'>sp1#sample0': 'sp1',
		'>sp2#sample0': 'sp2',
		'>sp3#sample0': 'sp3',
		'>sp1#sample1': 'sp1',
		'>sp4#sample0': 'sp4'
	}

	mn, tc = get_name_map(infiles, True, 0.4)
	assert mn == {
		'>sp0#sample0': 'sp0_sample0',
		'>sp1#sample0': 'sp1_sample0',
		'>sp2#sample0': 'sp2_sample0',
		'>sp3#sample0': 'sp3_sample0'
	}

def test_clean_name():
	assert clean_name(">sp3#sample|*0-sub  subsample") == "sp3_sample0_sub_subsample"
	assert clean_name("sample.valid..name_0") == "sample.valid..name_0"


name_map, term_count = utils.get_name_map(infiles, False)

def test_partition_instant():
	part0 = Partition(infiles[0], name_map)
	part1 = Partition(infiles[1], name_map)

	assert list(part0.data.keys()) == ['sp0', 'sp1', 'sp2', 'sp3']

	assert list(part1.data.keys()) == ['sp0', 'sp1', 'sp2', 'sp4']

	assert part0.data['sp0'] == 'TATTCCTCTATTA---GTAATTGGGCTTCTACTT-TTCCAAGAGCAACTAAAAATATTCGGCGTATC---'
	assert part0.data['sp3'] == '---TCCTCTATTAA-AGGAATTGGG--TCTACATTTTCCAAGAGCAACTACCCATATTC-GCGTATCTGG'

	assert part1.data['sp0'] == 'CTTCTTTGCATTTATTACGATCGATTCTCCATGAATG------TAGTTTTAGTAAAGAAAATTTGCAGAAATCTCTGATT'

	assert part1.data['sp4'] == 'CTTCGCTGCAT--ATTACGATC-ATTCTCCATGAATGATAATCTAGTTTTAGTTAAGAATATTTGCAGAAATCTCTGATT'


part0 = Partition(infiles[0], name_map)
part0.indel_coder()
part1 = Partition(infiles[1], name_map)
part1.indel_coder()

def test_partition_indel():
	
	assert part0.data['sp0'] == 'TATTCCTCTATTA---GTAATTGGGCTTCTACTT-TTCCAAGAGCAACTAAAAATATTCGGCGTATC---?001'

	assert part0.data['sp1'] == '-ATTCCTCTATTA---GTAATTGGGCTTCTACTT-TTCCATGAGCAACTAAAAATATACGGCGTATCT--?001'

	assert part0.data['sp2'] == '--TTCCTCTATTAA-AGGAATTGGG--TCTACATTTTCCACGAGCAACTACCCATATTCGGCGTATCTG-0110'

	assert part0.data['sp3'] == '---TCCTCTATTAA-AGGAATTGGG--TCTACATTTTCCAAGAGCAACTACCCATATTC-GCGTATCTGG0110'

	assert part0.metadata['size'] == [70, 4]

	assert part0.metadata['type'] == ['nucleic', 'indel']


def test_min_steps_char():

	assert utils.min_steps_char({"A": 1, "T": 2}, "nucleic") == 1
	assert utils.min_steps_char({"A": 1, "R": 2}, "nucleic") == 0
	assert utils.min_steps_char({"G": 1, "H": 2}, "nucleic") == 1
	assert utils.min_steps_char({"A": 1, "H": 2}, "nucleic") == 0


def test_max_steps_char():
	assert utils.max_steps_char({'A': 3, 'T': 4, 'G':50}, "nucleic") == 7
	assert utils.max_steps_char({'A': 3, 'T': 44, 'G':50}, "nucleic") == 47
	assert utils.max_steps_char({'A': 3, 'T': 44, 'R':50}, "nucleic") == 53
	assert utils.max_steps_char({'A': 3, 'T': 44, 'K':50}, "nucleic") == 3
	assert utils.max_steps_char({'A': 53, 'T': 4, 'W':50}, "nucleic") == 4
	assert utils.max_steps_char({'A': 3, 'T': 44, 'W':50}, "nucleic") == 3


def test_informative_stats():
	
	part0.informative_stats()
	part1.informative_stats()
	assert part0.metadata['informative_chars'][0] == [17, 32, 50, 51, 52]
	assert part0.metadata['informative_chars'][1] == [1, 2, 3]
	assert part1.metadata['informative_chars'][0] == [5, 53, 78]
	assert part1.metadata['informative_chars'][1] == [0]



def test_term_data():

	sp0dat = utils.Term_data('sp0')
	sp0dat.feed(part0)
	sp0dat.feed(part1)
	sp3dat = utils.Term_data('sp3')
	sp3dat.feed(part0)
	sp3dat.feed(part1)
	sp4dat = utils.Term_data('sp4')
	sp4dat.feed(part0)
	sp4dat.feed(part1)

	assert sp0dat.partition_table == [
		(70, 'nucleic', [17, 32, 50, 51, 52], True), 
		(4, 'indel', [1, 2, 3], True),
		(80, 'nucleic', [5, 53, 78], True),
		(3, 'indel', [0], True)]

	assert sp3dat.partition_table == [
		(70, 'nucleic', [17, 32, 50, 51, 52], True), 
		(4, 'indel', [1, 2, 3], True),
		(80, 'nucleic', [5, 53, 78], False),
		(3, 'indel', [0], False)]

	assert sp4dat.partition_table == [
		(70, 'nucleic', [17, 32, 50, 51, 52], False), 
		(4, 'indel', [1, 2, 3], False),
		(80, 'nucleic', [5, 53, 78], True),
		(3, 'indel', [0], True)]

sp0dat = Term_data('sp0')
sp3dat = Term_data('sp3')
sp4dat = Term_data('sp4')

def test_final_cleanup():

	for fi in infiles:
		os.remove(fi)

	os.remove(sp0dat.file)
	os.remove(sp3dat.file)
	os.remove(sp4dat.file)

