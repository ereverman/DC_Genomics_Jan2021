## Objectives:
* Get familiar with the R environment
* Identify the key advantages of using RStudio
* Create an R project

## R Advantages:
* Developed in 1995 and is an open source and completely free program
* Packages are written by different people with diverse research interests
  * There are packages with lots of generalizable functions
  * There are packages with very specific functions and uses
* R is widely used and as a result it is often easy to get help with code and functions
  * In the form of web forums, github, and talking to colleagues
* R is powerful and runs on multiple platorms
  * Can work with very large datasets
  * A speciality of many of the commands and packages we will use is data characterization
    * Summary and subsetting functions that help you know more about your data (especially if you didn't generate it)
  
 ## RSTUDIO:
 * IDE = Integrated Development Environment
 * Graphical interface combined with a scripting text editor built in that allows you to write, run, and save code
 
 
 ## Log on to RStudio Server:
 * I am most familiar with working with R on my local computer, but we will use one that is maintained on a web server by the carpentries
 
 ec2-3-227-219-253.compute-1.amazonaws.com:8787
 
 ## Brief tour of RStudio:
 * 3 windows appear when it loads
  * Console
  * Global Environment
  * File Browser

## Create a Project:
* Allows you to save files and code related to project-specific analyses
* You have the ability to start where you left off the analysis by saving the R environment
* Makes it easier to collaborate because you are working with a self-contained project unit rather than having code spread across many directories


1. File --> Menu --> New Project
2. Select New Directory --> New Project
3. Enter dc_genomics_r for directory name (leave project subdirectory as ~)
4. Click Create Project (Look for the new project in the file browser). All R projects end with .Rproj

## Making Projects reproducible:
* We've used a bunch of tools (Fastqc, trimmomatic, samtools, etc) that were loaded into our AMI for us
* R uses a similar bunch of tools called packages that we will use to run our analyses
* These packages are frequently updated, and sometimes the analysis you are working on depends on a particular version of the package
* Some packages have specific dependency requirements (require specific versions of other packages)
* Keeping track can get complicated
* We can use packrat to help keep track of packages used for a particular project

```
install.packages("packrat")
packrat::init()
# Restarts our R session with packrat enabled

# As you add and remove packages use
packrat::snapshot()
# To save a list of the packages that are required for your analyses
```

## Create an R script:
* Allows us to keep a record of the commands and allows for comments
* File --> New File --> Rscript
* Save as genomics_r_basics





