# Welcome to the GEDI

Gedi is a software platform for working with genomic data such as sequencing reads, sequences, per-base numeric values or annotations written in Java.

Gedi is developed in the [erhard-lab](https://erhard-lab.de).

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

This will create the folder defined in $INSTALL_DIR, and inside of $INSTALL_DIR, will generate a folder named gedi (containing the source code) and a folder named bin (containing the scripts and compiled program files). The second command just created a temporary variable containing the current working directory.

The bin directory should contain the following files
- bamlist2cit
- gedi
- Gedi-VERSION.jar (VERSION is the current version of gedi)
- gedi.jar (softlink)
- lib (subfolder)

**Important**: Pay attention that the folder where you execute 'mvn ...' is empty!


## How to run
There are two options:
- Always call it with its full path (use the full path to the bin folder generated above instead of <BIN-FOLDER>):
```bash
<BIN-FOLDER>/gedi -e Version
```

- or put this folder into your [PATH variable](https://en.wikipedia.org/wiki/PATH_(variable)): 
```bash
export PATH=$PATH:<BIN-FOLDER>
gedi -e Version
```
```



