# BAD2matrix
merging and translating FASTA alignments into TNT, FastTree, and RAxML/IQ-Tree with indel and gene content characters

### Install
Simply clone the Github repository or download the zip version. A Python 3 interpreter is required.

### Use

```bash
python bad2matrix.py -d <directory> -n <root-name> [ -a 2|3|4|5|6|6dso|6kgb|6sr|8|10|11|12|15|18|20 ] 
[ -f ] [ -g ] [ -i ] [ -m int ] [ -o speciesA,speciesB... ] [ -t <directory> ]
```



| option flag | description | required |
| --- | --- | --- |
|-a <i>2\|3\|4\|5\|6\|6dso\|6kgb\| 6sr\|8\|10\|11\|12\|15\|18\|20</i> | Number of amino acid states (default = 20). Reduction with option ‘6dso’ follows [Dayhoff et al. (1978)](http://chagall.med.cornell.edu/BioinfoCourse/PDFs/Lecture2/Dayhoff1978.pdf); option ‘6kgb’ follows [Kosiol et al. (2004)](https://doi.org/10.1016/j.jtbi.2003.12.010); option ‘6sr’ follows [Susko and Roger (2007)](https://doi.org/10.1093/molbev/msm144); option ‘11’ follows [Buchfink et al. (2015)](https://doi.org/10.1038/nmeth.3176); and all other options follow [Murphy et al. (2000)](https://doi.org/10.1093/protein/13.3.149). | no |
| -d <i>directory–name</i> | The input directory of aligned FASTA files. The default behavior aggregates read sequences of the same species across partitions, in which case names should use the following convention: ‘>species#sequenceID’. This also implies that a species could not be represented by more than a single read within each fasta file. Characters other than letters, numbers, periods, and underscores will be deleted. Use ‘-f’ for an alternate naming convention. | yes |
| -f | Use full FASTA names rather than default settings (see ‘-d’ description for default). Sequences from the same species but different read will not be aggregated and will be considered distinct OTUs. Characters other than letters, numbers, periods, and underscores will be deleted. | no |
| -g | Do not code gene content (absence/presence). If this flag is not set, gene content is coded. | no |
| -i | Do not code indels. If this flag is not set, indels are coded using simple indel coding following [Simmons and Ochoterena (2000)](https://doi.org/10.1080/10635159950173889). | no |
| -m <i>x</i> | Retain the upper <i>x</i> percentile of genes in the distribution of missing sequences. By default <i>x</i> = 1 (i.e. include all genes with four or more sequences). | no |
| -n <i>root–name</i> | Specify the root–name for output files. | yes |
| -o <i>list,of,terminals</i> | Specify a comma delimited list of terminals to use for rooting trees. If not specified, the first (alphabetically) is used. | no |
| -r _path_ | Folder containing a morphological matrix or a set of ortholog duplication matrices. Datasets should be saved a tsv tables. Multiple states should be separated by pipes (‘\|’) in the case of a polymorphic character. | no |
<i>Italicized text following option flags should be specified by the user.</i>

### Sample input/output
```bash
./BAD2matrix.pl -d test-data -f -g -n test -o Callitris_preisii_var._preisii,Callitris_sp.,Widdringtonia_nodiflora
```

### Citation
If you use this software, please cite: Little & Salinas. 2023. BAD2matrix: better phylogenomic matrix production: concatenation; indel coding; gene content coding; reduced amino acid alphabets; occupancy filtering.

### license
[GPL2](https://github.com/dpl10/BAD2matrix/blob/master/LICENSE)

### related repositories
* [2matrix](https://github.com/nrsalinas/2matrix)