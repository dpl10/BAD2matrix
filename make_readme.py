from bad2matrix import documentation

accu = '# BAD2matrix\n\n'

accu += documentation.pop('DOI') + '\n\n'

for heading in documentation:
    
	if heading != 'Preamble':
		accu += f'### {heading}\n\n'			

	if heading == 'Sample input/output':

		accu += f'```bash\n{documentation[heading]}\n```\n\n'

	elif heading == 'Contact':

		for author in documentation[heading]:
			accu += f"{author}  \n"

			for dat in documentation[heading][author]:
				accu += f"{dat}  \n"

			accu += '\n'

	elif heading == 'Usage':

		accu += f'```bash\n{documentation[heading]["command"]}\n```\n\n| option | description |\n| --- | --- |\n'

		for opt in documentation[heading]['options']:

			accu += f'{opt} | {documentation[heading]["options"][opt]}\n'

		accu += '\n'

	else:

		accu += f'{documentation[heading]}\n\n'

with open('README.md', 'w') as fhandle:
	fhandle.write(accu)

exit()
