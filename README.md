# Welcome to the GEDI

Gedi is a software platform for working with genomic data such as sequencing reads, sequences, per-base numeric values or annotations written in Java.

Gedi is developed in the [erhard-lab](https://erhard-lab.de).

Most of the tools and methods developed in the group are part of, or a plugin to Gedi (e.g. [LFC](https://www.ncbi.nlm.nih.gov/pubmed/26160885), [PRICE](https://www.nature.com/articles/nmeth.4631), [GRAND-SLAM](https://academic.oup.com/bioinformatics/article/34/13/i218/5045735), [Peptide-PRISM](https://aacrjournals.org/cancerimmunolres/article/8/8/1018/470266/Identification-of-the-Cryptic-HLA-I)).

For documentation, see [here](https://github.com/erhard-lab/gedi/wiki)!

# News
- 7.2.2024: We now use Maven instead of Ant as our build system, and do not support Java 8 anymore

# Installation
## Prerequisites

* [Java](https://openjdk.org/) 11 or higher
* [Maven](https://maven.apache.org/) for the build script

## Installing Gedi

Follow these steps for installing Gedi:

```bash
INSTALL_DIR="$HOME/gedi"  # change this to any folder you like to install gedi to
mkdir -p $INSTALL_DIR
cd $INSTALL_DIR

git clone git@github.com:erhard-lab/gedi.git

mkdir bin
cd bin

mvn -f $INSTALL_DIR/gedi/Gedi package 
```

This will create the folder defined in $INSTALL_DIR, and inside of $INSTALL_DIR, will generate a folder named gedi (containing the source code) and a folder named bin (containing the scripts and compiled program files) The second command just created a temporary variable containing the current working directory.


## How to run
There are two options:
- Always call it with its full path (use the bin folder generated above instead of <BIN-FOLDER>):
```bash
<BIN-FOLDER>/gedi -e Version
```

- or put this folder into your [PATH variable](https://en.wikipedia.org/wiki/PATH_(variable)): 
```bash
export PATH=$PATH:<BIN-FOLDER>
gedi -e Version
```
```



