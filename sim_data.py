import os
import subprocess
import shutil
import re

seq0 = "TGCGGAAG-ATCATTGTCGAAAACC-----AGCAGAAAACCCGCGAACTCGTCTGTACTCTTGGGAAA--"
seq1 = "-GCGGTTG-ATCATTGTCGAAAACCT---AAGCAGTTTTCCCGCGAACTCGTCTGTACTCTTGGGTTTTG"

gen_content = {
    'sp0': [1, 1, 1, 1, 1, 0],
    'sp1': [1, 1, 0, 1, 1, 1],
    'sp2': [1, 1, 1, 1, 0, 1],
    'sp3': [1, 1, 1, 1, 1, 1],
    'sp4': [1, 1, 1, 0, 1, 1],
    'sp5': [1, 0, 1, 1, 1, 1],
    'sp6': [1, 1, 1, 1, 1, 1],
    'sp7': [1, 1, 0, 1, 1, 1],
    'sp8': [1, 1, 1, 1, 1, 0],
    'sp9': [1, 1, 1, 1, 0, 1],        
}

def generate_datasets(dir_name = 'temporary_directory_for_testing'):

	if not os.path.exists(dir_name):
		os.mkdir(dir_name)

	for igen, gen in enumerate(gen_content['sp0']):
		bffr = ''
		
		for isp,sp in enumerate(gen_content):

			if gen_content[sp][igen] == 1:
				bffr += f'>{sp}\n'
				
				if isp < 5:
					bffr += f'{seq0}\n'
				
				else:
					bffr += f'{seq1}\n'
		
		file_name = os.path.join(dir_name, f'gene_{igen}.fasta')

		with open(file_name, 'w') as fhandle:
			fhandle.write(bffr)


def clean(dir_name, matrix_rootname):

	if os.path.exists(dir_name):
		shutil.rmtree(dir_name)

	if os.path.exists('fasttree_datasets'):
		shutil.rmtree('fasttree_datasets')

	if os.path.exists('iqtree_datasets'):
		shutil.rmtree('iqtree_datasets')

	if os.path.exists('raxml_datasets'):
		shutil.rmtree('raxml_datasets')

	if os.path.exists('tnt_datasets'):
		shutil.rmtree('tnt_datasets')


	#for d,s,f in os.walk('.'):
	#	for filito in f:
	#		if filito.startswith(matrix_rootname):
	#			os.remove(filito)


def verify_gene_content(matrix_rootname):
	out = True
	encoded = {}

	with open(f'iqtree_datasets/{matrix_rootname}_gene_content.phy', 'r') as fh:
	
		for line in fh:
			line = line.strip()
	
			if line.startswith('s'):
				bits = re.split(r'\s+', line)
				lista = [int(x) for x in bits[1]]
				encoded[bits[0]] = lista

	for key in encoded:
		
		if not key in gen_content:
			out = False
			break

		if encoded[key] != gen_content[key]:
			out = False
			break

	return out

if __name__ == '__main__':
	folder_name = 'temporary_directory_for_gene_content_testing'
	matrix_root = 'devtest'
	generate_datasets(folder_name)
	subprocess.call(['python' ,'bad2matrix.py','-d', folder_name, '-f', '-n', matrix_root])	
	result = verify_gene_content(matrix_root)
	print(f'{result=}')
	clean(folder_name, matrix_root)
	exit()
