## Background and Metadata:

# Background:
* Data come from E. coli, one of the most well studied model organisms in biology.
* E. coli:
  * Rod-shaped bacteria that are part of the normal microflora of the gut
  * Survive under variable temperature, nutrient, and oxygen levels
  * Typically harmless, but some strains are associated with food-poisoning
  
 * Why study E. coli:
  * Single-celled
  * Fast generation time, doubles every 20 minutes
  * Genetic complexity is reduced compared to other models
  * A plethora of genetic tools and protocols are available to study adaptation and evolution
  
  
# The Data:
* Come from a long-term evolution experiment by Richard Lenski:
  * 12 identical populations of E. coli were initiated in February 1988
  * GOAL: assess adaptation in E. coli
  * Cultures were maintained in a glucose-limited medium supplemented with citrate
    * Under aerobic conditions, E. coli normally doesn't use citrate
  * Populations were sequenced periodically to look for adaptive changes
    * Over time the 12 populations began to diverge
    * Between generation 31,000 and 31,500, a spontaneous citrate-using variant appeared
    * Also observed hypermutability in certain regions of the genome, which may lead to faster adaptation in novel environments
  * As of 2020, the populations have undergone 73,500 generations
  * Incidently, the populations were frozen because of the COVID-19 pandemic
  
# Our Data:
* One population: Ara-3
* Three sampling periods: 5,000 15,000 and 50,000 generations
* This population changed substantially through the experiment
* OUR GOAL is to identify and characterize these changes using a variant calling pipeline.

* Here are your metadata: https://raw.githubusercontent.com/datacarpentry/wrangling-genomics/gh-pages/files/Ecoli_metadata_composite.csv

* Explanation of spreadsheet headers:
  * strain = strain name
  * generation = generation when sample was frozen
  * clade = based on parsimony-based tree
  * reference = study the samples were originally sequenced for
  * population = ancestral population group
  * mutator = hypermutability mutant status
  * facility = facility samples were sequenced at
  * run = sequence read archive sample ID
  * read_type = library type of reads
  * read_length = length of reads in sample
  * sequencing_depth = depth of sequencing
  * cit = citrate-using mutant status

## DO THIS:
1. 
