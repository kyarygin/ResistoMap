## Introduction
[ResistoMap](http://resistomap.datalaboratory.ru/) - a Web-based interactive visualization of the presence of genetic determinants conferring resistance to antibiotics, biocides and heavy metals in human gut microbiota. ResistoMap displays the data about more than 1500 published gut metagenomes of the world populations including both healthy subjects and patients.

Scripts in this repo allow to estimate resistance potential to antibiotics, biocides and heavy metals for user-defined metagenomic samples.

## Dependencies
Scripts relies on [bowtie2](http://bowtie-bio.sourceforge.net/bowtie2/index.shtml) and [diamond](https://github.com/bbuchfink/diamond) software. Please, install it and put in your executable search path (or specify full pathes to executables in `config.json`).

Also, you need to install [_Pysam_](https://pysam.readthedocs.io/en/latest/) and [_Biopython_](http://biopython.org/) modules for python.

## Usage
To calculate resistance potential for your samples, just do
```
$ python run.py sample_1.fastq sample_2.fastq sample_3.fastq
```
Script will produce two tables: _gene_table_total.tsv_ , containing abundance values for resistance genes, and _ab_table_total.tsv_, containing resistance potential for antibiotics.

You can specify output folder with `-o` key, e.g:
```
$ python run.py -o /home/user/resistomap_data sample_1.fastq sample_2.fastq sample_3.fastq
```
and number of threads for bowtie and diamond with `-n` key:
```
$ python run.py -n 20 sample_1.fastq sample_2.fastq sample_3.fastq
```
Default output folder is current dir and default number of threads is 1.

Execute `$ python run.py -h` to see help message.
