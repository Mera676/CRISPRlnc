# CRISPRlnc_SVMModel

## 1 prepare the environment

```shell
conda env create -f py38.yaml
```

## 2 Train

### 2.1 CRISPRko

```shell
python CRISPRko/train.py CRISPRko/data/CRISPRko_train_data.pkl
```

### 2.2 CRISPRi

```shell
python CRISPRi/train.py CRISPRko/data/CRISPRi_train_data.pkl
```

## 3 Predict

### 3.1 CRISPRko

 Prepare the <sgRNA_predict.txt> file according to the following format, and place the file in the <CRISPRko/data> directory:

```txt
TGAATGGCGTCGGGGTTGGA
ATCCGAGTCACGGCATTGTG
ATCCGAGTCACGGCATTGTC
AAATTCTCTTAGGAGCTCCT
GCGTTCTCTTAGGAGCTCCT
GCGGAGGACCCGGAAGGATG
AGGGCTGCAGAGGGTATGGG
AATGCTGCACTCGGTATGGG
GTAAAGATGAAGGTTAATCA
TTCCGAGTCACCGCATTGTG
```

Then use the following statement to predict the CRISPRko effectiveness of sgRNA:

```shell
python CRISPRko/predict.py CRISPRko/data/sgRNA_predict.txt
```

### 3.2 CRISPRi

 Prepare the <sgRNA_predict.txt> file according to the following format, and place the file in the <CRISPRi/data> directory:

```txt
TGAATGGCGTCGGGGTTGGA
ATCCGAGTCACGGCATTGTG
ATCCGAGTCACGGCATTGTC
AAATTCTCTTAGGAGCTCCT
GCGTTCTCTTAGGAGCTCCT
GCGGAGGACCCGGAAGGATG
AGGGCTGCAGAGGGTATGGG
AATGCTGCACTCGGTATGGG
GTAAAGATGAAGGTTAATCA
TTCCGAGTCACCGCATTGTG
```

Then use the following statement to predict the CRISPRi effectiveness of sgRNA:

```shell
python CRISPRko/predict.py CRISPRi/data/sgRNA_predict.txt
```

## 4 Train the model with your data

### 4.1 CRISPRko

Prepare the <sgRNA_train.txt> file according to the following format, and place the file in the <CRISPRko/data> directory:

```txt
TGAATGGCGTCGGGGTTGGA,1
ATCCGAGTCACGGCATTGTG,0
ATCCGAGTCACGGCATTGTC,1
AAATTCTCTTAGGAGCTCCT,1
GCGTTCTCTTAGGAGCTCCT,1
GCGGAGGACCCGGAAGGATG,0
AGGGCTGCAGAGGGTATGGG,1
AATGCTGCACTCGGTATGGG,1
GTAAAGATGAAGGTTAATCA,0
TTCCGAGTCACCGCATTGTG,1
```

Then encode the new data using the following command:

```
python CRISPRko/emb.py CRISPRko/data/sgRNA_train.txt
```

Then you can go back to the second step to continue training.

### 4.2 CRISPRi

Prepare the <sgRNA_train.txt> file according to the following format, and place the file in the <CRISPRi/data> directory:

```txt
TGAATGGCGTCGGGGTTGGA,1
ATCCGAGTCACGGCATTGTG,0
ATCCGAGTCACGGCATTGTC,1
AAATTCTCTTAGGAGCTCCT,1
GCGTTCTCTTAGGAGCTCCT,1
GCGGAGGACCCGGAAGGATG,0
AGGGCTGCAGAGGGTATGGG,1
AATGCTGCACTCGGTATGGG,1
GTAAAGATGAAGGTTAATCA,0
TTCCGAGTCACCGCATTGTG,1
```

Then encode the new data using the following command:

```
python CRISPRi/emb.py CRISPRi/data/sgRNA_train.txt
```

Then you can go back to the second step to continue training.