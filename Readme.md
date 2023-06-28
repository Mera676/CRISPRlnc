# Using the CRISPRlnc prediction tool

## Preparation

Our tools are developed based on python in the Linux environment, so please make sure you can meet this requirement before using it.

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

ViennaRNA(version:2.5.0)

You can directly import the required packages into your current python environment with the following command:

```shell
conda env create -f CRISPRlnc.yml
```



### Required file comments:

First, we need the genome annotation files of the species, including the fasta sequence file of the whole genome, the gtf annotation file, and the fasta file of the non-coding gene sequence. If it is used to predict sgRNAs on **human\ Mus\Zebrafish\Drosophila melanogaster\\Paniscus\Gorilla\Macaque\Cow\Goat**  non-coding genes, **this step can be skipped**.

exmaple(hg38) :

- whole genome sequence : human_hg38.fa
- gtf annotation file : Homo_sapiens.GRCh38.107.chr.gtf
- non-coding gene sequence : Homo_sapiens.GRCh38.ncrna.fa

Then, you need to extract the location of the promoter by the following shell command:

```shell
sed 's/"/\t/g' Homo_sapiens.GRCh38.107.chr.gtf | awk 'BEGIN{OFS=FS="\t"}{if($3=="transcript") {if($7=="+") {start=$4-1000; end=$4+500;} else {if($7=="-") start=$5-500; end=$5+1000; } if(start<0) start=0; print $1,start,end,$14,$10,$7;}}' >GRCh38transcript.promoter.bed
```

### Off-target tool support

Our tool inherits the Cas-OFFinder tool as an off-target prediction, if you don't need the prediction results of this part, you can skip this part.

In our program, this tool has been inherited, and there is no need to download it again. You can try our prediction tool directly. If there is an error message, please check the official instructions of Cas-OFFinder. 

For detailed instructions, please refer to the official installation document.http://www.rgenome.net/cas-offinder/portable

## Instructions for use

We provide two versions of the sgRNA prediction tool: with off-target prediction and without off-target prediction.

#### no off-target prediction version:

##### input file(take human as an example):

Annotation/Homo_sapiens.GRCh38.ncrna.fa

Annotation/GRCh38transcript.promoter.bed

sg.fa(The only file you need to change is to replace the content with the gene sequence you want to predict)

##### If you don't need the result of off-target prediction, or you can't use Cas-OFFinder, you can use our tool directly through the following command,We prepared a total of nine genome files in our package:

###### human:

```shell
python CRISPRlnc_nooff.py Annotation/Homo_sapiens.GRCh38.ncrna.fa Annotation/GRCh38transcript.promoter.bed sg.fa
```

###### Mus:

```shell
python CRISPRlnc_nooff.py Annotation/Mus_musculus.GRCm39.ncrna.fa Annotation/GRCm39.transcript.promoter.bed sg.fa
```

###### Zebrafish:

```shell
python CRISPRlnc_nooff.py Annotation/Danio_rerio.GRCz11.ncrna.fa Annotation/GRCz11.transcript.promoter.bed sg.fa
```

###### Drosophila melanogaster:

```shell
python CRISPRlnc_nooff.py Annotation/Drosophila_melanogaster.BDGP6.32.ncrna.fa Annotation/BDGP6.transcript.promoter.bed sg.fa
```

###### Paniscus:

```shell
python CRISPRlnc_nooff.py Annotation/Pan_paniscus.panpan1.1.ncrna.fa Annotation/panpan1.1.transcript.promoter.bed sg.fa
```

###### Gorilla:

```shell
python CRISPRlnc_nooff.py Annotation/Gorilla_gorilla.gorGor4.ncrna.fa Annotation/gorGor4.transcript.promoter.bed sg.fa
```

###### Macaque:

```shell
python CRISPRlnc_nooff.py Annotation/Macaca_mulatta.Mmul_10.ncrna.fa Annotation/Mmul_10.transcript.promoter.bed sg.fa
```

###### Cow:

```shell
python CRISPRlnc_nooff.py Annotation/Bos_taurus.ARS-UCD1.2.ncrna.fa Annotation/ARS-UCD1.2.transcript.promoter.bed sg.fa
```

###### Goat:

```shell
python CRISPRlnc_nooff.py Annotation/Capra_hircus.ARS1.ncrna.fa Annotation/ARS1.transcript.promoter.bed sg.fa
```

###### other:

If you want to predict other species , please change the corresponding file location to the annotation file address of the species you need to predict, for example:

```shell
python CRISPRlnc_nooff.py ncRNA_sequence_path promoter_bedfile_path sg.fa
```



#### off-target prediction version(take human as an example):

##### input file:

Annotation/Homo_sapiens.GRCh38.ncrna.fa

Annotation/GRCh38transcript.promoter.bed

sg.fa(you need to replace the content inside with the gene sequence you want to predict)

fasta/human_hg38.fa 

##### Run order with off-target prediction:

```shell
python CRISPRlnc.py Annotation/Homo_sapiens.GRCh38.ncrna.fa Annotation/GRCh38transcript.promoter.bed sg.fa fasta/human_hg38.fa
```

If you want to predict species other than human, please change the corresponding file location to the annotation file address of the species you need to predict, for example:

```shell
python CRISPRlnc.py  ncRNA_sequence_path promoter_bedfile_path sg.fa whole_genome_path
```

