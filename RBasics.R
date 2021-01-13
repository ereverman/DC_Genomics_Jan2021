# R Basics

# Getting to Work with R:----

# Functions have a pattern. At simplest, command()

# Show working directory similar pwd
getwd()

# to run, command return on mac, control r on windows
# or click the run button

# Exercise:

dir() # Lists files in working directory
sessionInfo() # Version of R and information about packages
date() # gives date
Sys.time() # gives current time

# Commands are case sensitive
###

# Getting help:----

round(3.14) # Rounds to nearest integer

# I want to round to a different numbers of digits:
?round() # opens a new page in the file browser window in the Help tab
args(round) # provides the syntax for implementing the function, default digits is 0
# You can also often see the syntax by hovering your mouse over the function (at least on a mac)

round(3.1456789, digits = 2) # use tab complete to help autofill the function options

# We used ? to look at a function that comes pre-installed in R
?geom.point() # not installed, we get an error
??geom.point() # Still not helpful for geompoint, but opens the help menu

help.search("chisquare test") # searches the help tab 

# Contextual help:
lm # hit tab and a menu appears that you can scroll through with arrow keys
lm() # hit tab inside the parentheses and a menu with arguments appears with additional help

###

# What we aren't going to cover:----
# all of the forms of data (matrices and lists)
# writing loops
# writing custom functions

# BUT all of these things can be done in R, and they are similar to what we've done in the shell


# Creating objects in R
# Shell does work in your filing system, and it can see all of the files that are stored there
# R works similarly to the shell except we first need to store or load objects in the R environment
# for R to "see" them and work with them

# ppt slide for syntax

a <- 1 # shortcut for arrow: option - (mac) alt - (windows)
a

b <- bat
b <- "bat"
b

B <- b # now we have two objects that store the word bat
B


# Exercise:
human_chr_number <- 23
gene_name <- 'Ccs'
ensembl_url <- 'ftp://ftp.ensemblgenomes.org/pub/bacteria/release-39/fasta/bacteria_5_collection/escherichia_coli_b_str_rel606/'
human_diploid_number <- 2 * human_chr_number


# Reassigning object names overwrites them:
gene_name <- 'tp53'

rm(gene_name)

###

# Object data types:----
# every object in R has two properties: length and classification
# Most common data types: numeric, character, logical

# Exercise:
chromosome_name <- 'chr02'
od_600_value <- 0.47
chr_position <- '1001701'
spock <- TRUE
pilot <- Earhart # What do we need to do to fix this?
pilot <- 'Earhart'

class(chromosome_name)
class(od_600_value)
class(chr_position)
class(spock)
class(pilot)

###

# Mathematical and functional operations on objects:----
# R follows the order of operations for running calculations:

(1 + (5 * 0.5))/2

# Exponents two ways:
2**2
2^2

# Numbers stored as objects can be used in operations:
human_chr_number
human_chr_number * 2

###

# Vectors:----
# Objects can store more than one element, then they become vectors.
# ppt slide

snp_genes <- c("oxtr", "ACTN3","AR","OPRM1")
class(snp_genes)
length(snp_genes)
str(snp_genes) # structure function 

###

# Creating and subsetting vectors:----
# vectors can be used to store information and parts of the vectors can be retrieved as needed

snps <- c('snp1','snp2','snp3','snp4')
snp_chr <- c('3', '11', 'X', '6')
snp_positions <- c(20000,4000,56000,10000)


# When subsetting, often use []
# get the 3rd value in the snp_genes vector:
snp_genes[3]
# Each item in a vector is indexed. If you know the position of the element you want, you can call for it with the number of the position

snp_genes[1:3] # colon indicates you want the range between and including positions 1 and 3

snp_genes[c(1,3,4)]

# We can add information to the snp_genes or any vector:
snp_genes <- c(snp_genes, "Ccs","MTF1")
snp_genes

# Logical subsetting:

snp_positions
snp_positions[snp_positions > 15000]
snp_positions[snp_positions <= 15000]

# What is actually happening:
snp_positions > 15000
which(snp_positions > 15000) # returns the index positions of the TRUE positions

# In the context of genomics analyses, we may assign a cutoff value and subset data based on the cutoff:
snp_marker_cutoff <- 15000
snp_positions[snp_positions > snp_marker_cutoff]


# Exercise:

combined <- c(snp_genes[1], snps[1], snp_chr[1], snp_positions[1])
combined
###

# Factors and dataframes:----
# Data tables can be saved in R to be used in analyses similar to assigning objects
# Keeping your data tables in a Tidy format is essential practice, but practically will make your life much easier in R
# R is installed with several ways to read in data:
read.# hit tab to see them
args(read.csv)

# Read in a variant file:
variants <- read.csv("../r_data/combined_tidy_vcf.csv")
# Appears in the environment window
# can look at the data by clicking the name, or
View(variants)
dim(variants) # rows then columns
head(variants)
tail(variants)


# Summarize and classify data frame structure:
class(variants) # in contrast to a matrix (where all elements are one data type - number or character)
summary(variants)
# returns a summary for each of the columns in the data frame

str(variants) # defines the data class for each of the columns

###
# Factors:----
# the final data structure class we will talk about
# Factors are important because we will use them to define categories in our data
# Important for subsetting and plotting data

# Let's extract one of the columns to practice with:
Samples <- variants$sample_id # $ symbol indicates you are looking for a column in the data frame

head(Samples)
str(Samples) # currently a character

# To change to a factor with categories:
Samples <- as.factor(Samples) # Uses coersion to change the data type
str(Samples) # now indicates that there are 59 levels or categories in the REF variable
plot(Samples) # kind of an odd plot, but shows how the categories are used as the x axis and the number of observations in each category is plotted

# Reorder levels:
levels(Samples)
Samples <- factor(Samples, levels = c("SRR2589044", "SRR2584863","SRR2584866"))
plot(Samples)

###

# Subsetting data frames:----
# ppt slide
variants[1,1]
variants[,2] # all contents in column 2
variants[3,] # all contents in row 3

variants[1:4, c("REF","ALT")] # Can combine index and column names

# Obtain a subset of the data frame:

SRR2584863_variants <- variants[variants$sample_id == "SRR2584863", ]
dim(variants)
dim(SRR2584863_variants)
View(SRR2584863_variants)

# the other method:
SRR2589044_variants <- subset(variants, sample_id == "SRR2589044")
View(SRR2589044_variants)

###
# Writing or saving data frames:----
write.csv(SRR2584863_variants, file = "SRR2584863_variants.csv")
###

# Importing data from excel:----
# Easiest method is to save your excel file as a .csv file
# Then you can read the file in with read.csv

# Alternatively: 
# File > Import Dataset > from excel
# requires package installation. Choose yes. Progress bar indicates how much is left
# We don't have any excel files in our instances to work with but you'd be able to browse to find the file.

###

# Assessing data with dplyr:----
# Bracket subsetting and using the subset command can get unwieldy
# the package dplyr has several functions that are useful for manipulating data frames

# Install dplyr:
install.packages("dplyr")

# load or activate the package:
library("dplyr")

# dplyr is built to work with data frames, and is especially good to very large data frames


# Select columns with dplyr:----
colnames(variants) # helpful to remind ourselves what the column names are

# SYNTAX: select(DATAFRAME, col1, col2, ...)

select(variants, sample_id, REF, ALT, DP) # select specific columns

select(variants, -CHROM) # select all columns except CHROM

# To choose rows:

# SYNTAX: filter(DATAFRAME, COLUMN_NAME == "...")
filter(variants, sample_id == "SRR2584863")

# To combine select and filter, use pipes
# these work the same as in shell, but the symbol is different
# ppt

variants %>%
  filter(sample_id == "SRR2584863") %>%
  select(REF, ALT, DP) %>%
  head()
###

# Creating data-dependent new columns:----
# mutate can be used to add a new column to a data frame that is derived from a calculation (for example) using data already in the table:
colnames(variants)

# Column QUAL has the phred-scaled confidence score that a polymorphism exists 
# Convert to a more "human-readable" number with a calculation

# Probability = 1-10^-(QUAL/10)

variants %>%
  mutate(POLPROB = 1 - (10^ -(QUAL/10))) %>%
  head()
###

# Analyzing data by group:----
# group_by() allows you to identify categories, similar to factors

variants %>%
  group_by(sample_id) %>%
  summarize(n())

# n() finds the count for each id group
# Other functions can be used here:

variants %>%
  group_by(sample_id) %>%
  summarize(max(DP))

# Here we are looking at the max depth of coverage at for each of the variants by group


