## Introduction
[ResistoMap](http://resistomap.datalaboratory.ru/) - a Web-based interactive visualization of the presence of genetic determinants conferring resistance to antibiotics, biocides and heavy metals in human gut microbiota. ResistoMap displays the data about more than 1500 published gut metagenomes of the world populations including both healthy subjects and patients.

The scripts provided in this repository allow to estimate the resistance potential to antibiotics, biocides and heavy metals for user-defined metagenomic samples.

## Dependencies
Scripts relies on [bowtie2](http://bowtie-bio.sourceforge.net/bowtie2/index.shtml) and [diamond](https://github.com/bbuchfink/diamond) software. Please install both and put them in your executable search path (or specify the full paths to executables in `config.json`).

Also, you will need to install [_Pysam_](https://pysam.readthedocs.io/en/latest/) and [_Biopython_](http://biopython.org/) modules for python.

## Usage
To calculate resistance potential for your samples, run

```
$ python run.py sample_1.fastq sample_2.fastq sample_3.fastq
```

The script will produce two tables: _gene_table_total.tsv_ - containing the abundance values for resistance genes - and _ab_table_total.tsv_ - containing the resistance potential for antibiotics.

You can specify output folder with `-o` key, e.g:
```
$ python run.py -o /home/user/resistomap_data sample_1.fastq sample_2.fastq sample_3.fastq
```
and number of threads for Bowtie2 and DIAMOND with the `-n` key - for example, to use 20 threads:
```
$ python run.py -n 20 sample_1.fastq sample_2.fastq sample_3.fastq
```
The default output folder is the current directory and the default number of threads is 1.

Execute `$ python run.py -h` to see the help message.


## Upload your data

Create a sample metadata table with following fields:
* Sample name
* Country of origin
* Host gender
* Host age
* Clinical status / diagnosis

If you miss some sample data, leave NA values in fields, e.g.:

| Sample    | Country | Gender | Age | Status             |
| --------- | ------- | ------ | --- | ------------------ |
| O2.UC20-0 | Spain   | Female | 51  | Ulcerative colitis |
| O2.UC52-2 | Spain   | Female | NA  | Healthy            |
| V1.CD15-3 | Spain   | Female | NA  | Crohns disease     |

Calculate abundance tables for samples and send us metadata table along with _ab_table_total.tsv_, _gene_table_total.tsv_ to [resistomap@gmail.com](mailto:resistomap@gmail.com).
After internal validation, the datasets will be added to the Resistomap together with the information about contributors and reference publication.
