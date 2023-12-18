# Using the CRISPRlnc prediction tool

## Preparation

Our tools are developed based on python in the Linux environment, so please make sure you can meet this requirement before using it. Please download the annotation file from our website first（http://predict.crisprlnc.cc/）. You can also use our online design web tool at http://predict.crisprlnc.cc/ .

### Additional python packages:

sys

os

biopython  (version:1.79)

imbalanced-learn  (version:0.9.1)

imblearn  (version:0.0)

importlib-metadata   (version:4.12.0)

joblib  (version:1.1.0)

numpy(version:1.21.5)

pandas (version:1.4.2)

scikit-learn(version:1.1.1)

scipy  (version:1.7.3)

sklearn(version:0.0)

ViennaRNA(version:2.6.2)

You can directly import the required packages into your current python environment with the following command:

```shell
conda env create -f CRISPRlnc.yml
```

### Required file comments:

Firstly, you need to prepare the sequence files of non-coding genes and gene annotation files (recommended to download from ensemble database), and we prepared the data of nine genomes as reference (download from http://predict.crisprlnc.cc/tool_download). If you need to design other genomes, please go to ensemble database to download them by yourself. The total number of files needed is:

exmaple (Human) :

- whole genome sequence : hg38.fa
- gtf annotation file : Homo_sapiens.GRCh38.gtf
- non-coding gene sequence file : Homo_sapiens.GRCh38.ncrna.fa

When the above file is ready, run the following command:

```shell
python embeding.py non-coding_gene_sequence_file_path gtf_annotation_file_path whole_genome_sequence_file_path species_name
```

### Off-target tool support

Our tool inherits the Cas-OFFinder tool as an off-target prediction, if you don't need the prediction results of this part, you can skip this part.

In our program, this tool has been inherited, and there is no need to download it again. You can try our prediction tool directly. If there is an error message, please check the official instructions of Cas-OFFinder. 

For detailed instructions, please refer to the official installation document.http://www.rgenome.net/cas-offinder/portable

## Instructions for use

We provide two versions of the sgRNA prediction tool: with off-target prediction and without off-target prediction.

#### no off-target prediction version:

##### input parameter(take human as an example):

sg.fa (gene sequences to be designed for sgRNAs )

genome.fa(reference genome)

distance (length of promoter sequence)

If you don't need the result of off-target prediction, or you can't use Cas-OFFinder, you can use our tool directly through the following command:

##### CRISPRko(CPU):

```shell
python CRISPRko.py sg.fa genome.fa distance 0 0
```

##### CRISPRi(CPU):

```shell
python CRISPRi.py sg.fa genome.fa distance 0 0
```

##### CRISPRko(GPU):

```shell
python CRISPRko.py sg.fa genome.fa distance 0 1
```

##### CRISPRi(GPU):

```shell
python CRISPRi.py sg.fa genome.fa distance 0 1
```

#### off-target prediction version:

##### input parameter(take human as an example):

sg.fa (gene sequences to be designed for sgRNAs )

genome.fa(reference genome)

distance (length of promoter sequence)

##### CRISPRko(CPU):

```shell
python CRISPRko.py sg.fa genome.fa distance 1 0
```

##### CRISPRi(CPU):

```
python CRISPRi.py sg.fa genome.fa distance 1 0
```

##### CRISPRko(GPU):

```shell
python CRISPRko.py sg.fa genome.fa distance 1 1
```

##### CRISPRi(GPU):

```shell
python CRISPRi.py sg.fa genome.fa distance 1 1
```

