# Project Organization

## Objectives:
* Create a file system for a bioinformatics project
* Identify which file types should be store in which places
* Use tools to document work

## Starting your project:
* Project organization is a critical component of a bioninformatics project
  * There are a lot of raw and processed data files
  * Can become difficult to keep track, especially if the project lasts for several months or years
  * You will find yourself needing to refer back to earlier analyses as you work on data and begin writing.
* Arriving on a project structure that works for your analyses takes some practice and trial and error.
* Operating by the rule that each project should be kept as a self-contained unit will benefit you and collaborators

## Setting up a file system:
* We will work with a simplified example of a file system
```
cd # navigate home
pwd

# Make the following directories:
dc_workshop
dc_workshop/docs
dc_workshop/data
dc_workshop/results

rm -r dc_workshop

```
* I've learned that a useful strategy is to organize each project in largely the same way
* We can facilitate this by building a special customized program that will set up our directory structure
* This will involve a new tool or program called nano
  * nano is a simple text editor
  * It can be useful for editing documents or writing short scripts, but it isn't anyone's favorite text editor.
```
nano

# Bioinformatics project template

# Instructions:
# CHANGE the PROJECT_DIR directory path in the script below:

PROJECT_DIR=/home/dcuser/dc_worshop/
cd ${PROJECT_DIR}

# DATA:
DATA_DIR=${PROJECT_DIR}/data
DATA_RAW=${DATA_DIR}/raw
DATA_PROCESSED=${DATA_DIR}/processed

# Results directories
RESULTS_DIR=${PROJECT_DIR}/results

# Documents
DOCS_DIR=${PROJECT_DIR}/docs

mkdir -p ${PROJECT_DIR} \
         ${DATA_DIR} \
         ${DATA_RAW} \
         ${DATA_PROCESSED} \
         ${RESULTS_DIR} \
         ${DOCS_DIR}
         
```
