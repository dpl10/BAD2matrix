# BAD2matrix

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.10028408.svg)](https://doi.org/10.5281/zenodo.10028408)

A Python script for merging and translating FASTA alignments into TNT (extended XREAD), RAxML-NG/IQ-Tree (extended PHYLIP) and FastTree (fasta) input matrices. Optionally, it encodes indel characters using the `simple` gap coding method of Simmons and Ochoterena (2000; Gaps as characters in sequence-based phylogenetic analysis. Systematic Biology 49: 369-381. DOI 10.1080/10635159950173889), and gene content as binary characters (absence/presence). This script is slower than 2matrix.pl due to more disk use, but will not run out of RAM (hopefully).

### Installation

Simply clone the GitHub repository or download the main script (`bad2matrix.py`). A Python 3 interpreter is required.

### Usage

```bash
python bad2matrix.py -d <directory> -n <root-name> [-a 2|3|4|5|6|6dso|6kgb|6sr|8|10|11|12|15|18|20] [-f] [-g] [-i] [-m int] [-t <directory>]
```

| option | description |
| --- | --- |
-a | Number of amino acid states (default = 20). Reduction with option `6dso` follows [Dayhoff et al. (1978)](http://chagall.med.cornell.edu/BioinfoCourse/PDFs/Lecture2/Dayhoff1978.pdf); option `6kgb` follows [Kosiol et al. (2004)](https://doi.org/10.1016/j.jtbi.2003.12.010); option `6sr` follows [Susko and Roger (2007)](https://doi.org/10.1093/molbev/msm144); option `11` follows [Buchfink et al. (2015)](https://doi.org/10.1038/nmeth.3176); and all other options follow [Murphy et al. (2000)](https://doi.org/10.1093/protein/13.3.149).
-d | The input directory of aligned FASTA files. The default behavior aggregates sequences of the same species across partitions, in which case names should use the following convention: `>species#sequenceID`. This implies that a species can be represented by only one read within each FASTA file. Characters other than letters, numbers, periods, and underscores will be deleted. Use `-f` for an alternate naming convention.
-f | Use full FASTA names rather than default settings (see `-d` description for default). Sequences from the same species but different reads will not be aggregated and will be considered distinct OTUs. Characters other than letters, numbers, periods, and underscores will be deleted.
-g | Do not code gene content (absence/presence). If this flag is not set, gene content is coded.
-i | Do not code indels. If this flag is not set, indels are coded using simple indel coding following [Simmons and Ochoterena (2000)](https://doi.org/10.1080/10635159950173889).
-m | Retain the upper `x` percentile of genes in the distribution of missing sequences. By default `x` = 1 (i.e. include all genes with four or more sequences).
-n | Specify the root-name for output files.
-r | Folder containing a morphological matrix or a set of ortholog duplication matrices. Datasets should be saved a .tsv tables. Multiple states (polymorphic characters) should be separated by pipes (`\|`).

### Input specification

The program only accepts alignments in fasta format, and their file extensions should be `fa`, `fas`, or `fasta`. Alignments should also be located in a single folder, which path is parsed with the option `-d`.

Other kind of data (v. gr., a morphological matrix or multiple ortholog encoding tables) can also be concatenated. The first row of each table should contain the character names, and the first column the names of the terminals. These matrices should be saved as tab-separated values: simple text tables with tabs as column separators and `tsv` extension. These tables should be located in a single directory, which is parsed with the option `-r`. Polymorphisms should be separated by pipes (`|`): `leaves_pinnate|leaves_bipinnate`. Be aware that polymorphisms are only fully supported for the TNT output, they will be encoded as missing data in the rest of output datasets.

### Sample input/output

```bash
python bad2matrix.py -d test-data/fastas -f -g -n test
```

### Citation

Little, D. P. & N. R. Salinas. 2023. BAD2matrix: better phylogenomic matrix concatenation, indel coding, gene content coding, reduced amino acid alphabets, and occupancy filtering. Software distributed by the authors. DOI: 10.5281/zenodo.10028408.

### License

[GPL2](https://github.com/dpl10/BAD2matrix/blob/master/LICENSE)

### Related repositories

[2matrix](https://github.com/nrsalinas/2matrix)

### Contact

Nelson R. Salinas  
nrsalinas@gmail.com  
New York Botanical Garden  

Damon P. Little  
dlittle@nybg.org  
New York Botanical Garden  

