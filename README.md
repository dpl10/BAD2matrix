# BAD2matrix
merging and translating FASTA alignments into TNT, FastTree, and RAxML/IQ-Tree with indel and gene content characters

### install


### use
| option flag | description | required |
| --- | --- | --- |
|-a <i>2\|3\|4\|5\|6\|6dso\|6kgb\|6sr\|8\|10\|11\|12\|15\|18\|20</i> | Number of amino acid states (default = 20). Reduction with option ‘6dso’ follows [Dayhoff et al. (1978)](http://chagall.med.cornell.edu/BioinfoCourse/PDFs/Lecture2/Dayhoff1978.pdf); option ‘6kgb’ follows [Kosiol et al. (2004)](https://doi.org/10.1016/j.jtbi.2003.12.010); option ‘6sr’ follows [Susko and Roger (2007)](https://doi.org/10.1093/molbev/msm144); option ‘11’ follows [Buchfink et al. (2015)](https://doi.org/10.1038/nmeth.3176); and all other options follow [Murphy et al. (2000)](https://doi.org/10.1093/protein/13.3.149). | no |
| -d <i>directory–name</i> | The input directory of aligned FASTA files. Names should use the following convention: ‘>species#sequenceID’. Characters other than letters, numbers, periods, and underscores will be deleted. Use ‘-f’ for an alternate naming convention. | yes |
| -f | Use full FASTA names rather than default settings (see ‘-d’ description for default). Characters other than letters, numbers, periods, and underscores will be deleted. | no |
| -g | Do not code gene content (absence/presence). If this flag is not set, gene content is coded. | no |
| -i | Do not code indels. If this flag is not set, indels are coded. | no |
| -m <i>x</i> | Retain the upper <i>x</i> percentile of genes in the distribution of missing sequences. By default <i>x</i> = 1 (i.e. include all genes with four or more sequences). | no |
| -n <i>root–name</i> | Specify the root–name for output files. | yes |
| -o <i>list,of,terminals</i> | Specify a comma delimited list of terminals to use for rooting trees. If not specified, the first (alphabetically) is used. | no |
| -r | Specify the use of the original RAxML ‘.part’ file format (by default, the RAxML-NG/IQ-Tree format is used). | no |
<i>Italicized text following option flags should be specified by the user.</i>

### citation
If you use this software, please cite: Little et al. Submitted. BAD2matrix: better phylogenomic matrix production: concatenation; indel coding; gene content coding; reduced amino acid alphabets; occupancy filtering. [Applications in Plant Sciences.](https://doi.org/ADD_DOI_HERE)

### license
[GPL2](https://github.com/dpl10/BAD2matrix/blob/master/LICENSE)
