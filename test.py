import os
import utils

infiles = ['alg_test_0.fasta', 'alg_test_1.fasta']

dummy_0 = '''
>sp0#sample0
TATTCCTCTATTA---GTAATTGGGCTTCTACTT-TTCCAAGAGCAACTAAAAATATTCGGCGTATC---
>sp1#sample0
-ATTCCTCTATTA---GTAATTGGGCTTCTACTT-TTCCATGAGCAACTAAAAATATACGGCGTATCT--
>sp2#sample0
--TTCCTCTATTAA-AGGAATTGGG--TCTACATTTTCCACGAGCAACTACCCATATTCGGCGTATCTG-
>sp3#sample0
---TCCTCTATTAA-AGGAATTGGG--TCTACATTTTCCAAGAGCAACTACCCATATTC-GCGTATCTGG
'''

dummy_1 = '''
>sp0#sample0
CTTCTTTGCATTTATTACGATCGATTCTCCATGAATG------TAGTTTTAGTAAAGAAAATTTGCAGAAATCTCTGATT
>sp1#sample1
CTTCGTTGCATTTATTACGATC-ATTCTCCATGAATG------TAGTTTTAGTAAAGAATATTTGCAGAAATCTCTGACT
>sp2#sample0
CTTCGCTGCAT--ATTACGATC-ATTCTCCATGAATG------TAGTTTTAGTTAAGAATATTTGCAGAAATCTCTGACT
>sp4#sample0
CTTCGCTGCAT--ATTACGATC-ATTCTCCATGAATGATAATCTAGTTTTAGTTAAGAATATTTGCAGAAATCTCTGATT
'''

with open(infiles[0], 'w') as fh:
	fh.write(dummy_0)

with open(infiles[1], 'w') as fh:
	fh.write(dummy_1)


def test_dummy_files():
	assert os.path.exists(infiles[0]) == True
	assert os.path.exists(infiles[1]) == True

def test_name_map():
	assert utils.get_name_map(infiles, True)  == {
		'>sp0#sample0': 'sp0_sample0',
		'>sp1#sample0': 'sp1_sample0',
		'>sp2#sample0': 'sp2_sample0',
		'>sp3#sample0': 'sp3_sample0',
		'>sp1#sample1': 'sp1_sample1',
		'>sp4#sample0': 'sp4_sample0'
	}

	assert utils.get_name_map(infiles, False)  == {
		'>sp0#sample0': 'sp0',
		'>sp1#sample0': 'sp1',
		'>sp2#sample0': 'sp2',
		'>sp3#sample0': 'sp3',
		'>sp1#sample1': 'sp1',
		'>sp4#sample0': 'sp4'
	}

def test_clean_name():
	assert utils.clean_name(">sp3#sample|*0-sub  subsample") == "sp3_sample0_sub_subsample"
	assert utils.clean_name("sample.valid..name_0") == "sample.valid..name_0"


name_map = utils.get_name_map(infiles, False)

def test_partition_instant():
	part0 = utils.Partition(infiles[0], name_map)
	part1 = utils.Partition(infiles[1], name_map)

	assert list(part0.data.keys()) == ['sp0', 'sp1', 'sp2', 'sp3']

	assert list(part1.data.keys()) == ['sp0', 'sp1', 'sp2', 'sp4']

	assert part0.data['sp0'] == 'TATTCCTCTATTA---GTAATTGGGCTTCTACTT-TTCCAAGAGCAACTAAAAATATTCGGCGTATC---'
	assert part0.data['sp3'] == '---TCCTCTATTAA-AGGAATTGGG--TCTACATTTTCCAAGAGCAACTACCCATATTC-GCGTATCTGG'

	assert part1.data['sp0'] == 'CTTCTTTGCATTTATTACGATCGATTCTCCATGAATG------TAGTTTTAGTAAAGAAAATTTGCAGAAATCTCTGATT'

	assert part1.data['sp4'] == 'CTTCGCTGCAT--ATTACGATC-ATTCTCCATGAATGATAATCTAGTTTTAGTTAAGAATATTTGCAGAAATCTCTGATT'


part0 = utils.Partition(infiles[0], name_map)
part0.indel_coder()
part1 = utils.Partition(infiles[1], name_map)
part1.indel_coder()

def test_partition_indel():
	
	assert part0.data['sp0'] == 'TATTCCTCTATTA---GTAATTGGGCTTCTACTT-TTCCAAGAGCAACTAAAAATATTCGGCGTATC---?001'

	assert part0.data['sp1'] == '-ATTCCTCTATTA---GTAATTGGGCTTCTACTT-TTCCATGAGCAACTAAAAATATACGGCGTATCT--?001'

	assert part0.data['sp2'] == '--TTCCTCTATTAA-AGGAATTGGG--TCTACATTTTCCACGAGCAACTACCCATATTCGGCGTATCTG-0110'

	assert part0.data['sp3'] == '---TCCTCTATTAA-AGGAATTGGG--TCTACATTTTCCAAGAGCAACTACCCATATTC-GCGTATCTGG0110'

	assert part0.metadata['size'] == [70, 4]

	assert part0.metadata['type'] == ['nucleid', 'indel']


def test_informative_stats():
	
	part0.informative_stats()
	part1.informative_stats()
	assert part0.metadata['informative_chars'][0] == [17, 32, 50, 51, 52]
	assert part0.metadata['informative_chars'][1] == [1, 2, 3]


def test_min_steps_char():

	assert utils.min_steps_char({"A": 1, "T": 2}, "nucleic") == 2
	assert utils.min_steps_char({"A": 1, "R": 2}, "nucleic") == 1
	assert utils.min_steps_char({"G": 1, "H": 2}, "nucleic") == 2
	assert utils.min_steps_char({"A": 1, "H": 2}, "nucleic") == 1


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
		(70, 'nucleid', [17, 32, 50, 51, 52], True), 
		(4, 'indel', [1, 2, 3], True),
		(80, 'nucleid', [5, 53, 78], True),
		(3, 'indel', [0], True)]

	assert sp3dat.partition_table == [
		(70, 'nucleid', [17, 32, 50, 51, 52], True), 
		(4, 'indel', [1, 2, 3], True),
		(80, 'nucleid', [5, 53, 78], False),
		(3, 'indel', [0], False)]

	assert sp4dat.partition_table == [
		(70, 'nucleid', [17, 32, 50, 51, 52], False), 
		(4, 'indel', [1, 2, 3], False),
		(80, 'nucleid', [5, 53, 78], True),
		(3, 'indel', [0], True)]

sp0dat = utils.Term_data('sp0')
sp3dat = utils.Term_data('sp3')
sp4dat = utils.Term_data('sp4')

def test_final_cleanup():
	os.remove('alg_test_0.fasta')
	os.remove('alg_test_1.fasta')
	os.remove(sp0dat.file)
	os.remove(sp3dat.file)
	os.remove(sp4dat.file)

