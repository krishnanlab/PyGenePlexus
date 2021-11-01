#this script filters input file, adds names from terms.tsv and writes in gmt format
#input file contains these columns curated from mygeneinfo files: 
#ensembl_id, ontology, doid_id, doid_name, evidence, pubmed
#for coexp proj I used evidence codes: EXP_IDA_IPI_IMP_IGI_TAS_IC (this was for Kaylas first project)
tic <- as.integer(as.POSIXct(Sys.time()))

# This script will hard code what doid obo files will be used

library(tidyverse)
base_dir  <- '/mnt/research/compbio/krishnanlab/projects/GenePlexus/repos/GenePlexusBackend/data_hidden2/GSCs/'
args <- commandArgs(TRUE)

FN_direct <- paste0(base_dir,'direct_annotations_doid.tsv')
#read in file as tibble with all character columns
direct_tibble <- read_delim(FN_direct, 
                                delim = "\t", 
                                col_names = TRUE,
                                col_types = cols(.default = "c")) #can do pubmed = "d" if better


colnames(direct_tibble)[1] <- "gene_id" #rename col one to gene_id so this works for ensembl or entrez file
colnames(direct_tibble)[2] <- "doid_id" # This needs to be changed becuase Chris named the files differently

#read in terms file
FN_terms <- paste0(base_dir,'doid_terms.tsv')
terms_tibble <- read_delim(FN_terms, delim = "\t", col_names = F)
colnames(terms_tibble) <- c("doid_id", "doid_name")

#read file with ontology terms and ancestors
FN_ancestor <- paste0(base_dir,'doid_ancestors.tsv')
ancestor_tibble <- read_delim(FN_ancestor, delim = "\t", col_names = F)
#give colnames to ancestor tibble
colnames(ancestor_tibble) <- c("doid_id", "ancestors")
#make ancestors tibble tidy
ancestor_tibble <- separate_rows(ancestor_tibble, ancestors, sep = ", ")
filtered_tibble <- direct_tibble[!duplicated(direct_tibble),]


# This section is the one that does the propagation
joined_tibble <- left_join(filtered_tibble, ancestor_tibble, by = "doid_id")
gene_doid_tibble <- joined_tibble %>% select(gene_id, doid_id)
gene_anc_tibble <- joined_tibble %>% select(gene_id, ancestors)
#rename cols of gene_anc so we can bind rows
colnames(gene_anc_tibble) <- c("gene_id", "doid_id")
#bind rows to get all gene associations - direct and ancestor relationships
joined_tibble <- bind_rows(gene_doid_tibble, gene_anc_tibble)
#remove duplicated rows
joined_tibble <- joined_tibble[!duplicated(joined_tibble),]
#get rid of NA rows
joined_tibble <- na.omit(joined_tibble)


#get column gene_count
# this will add the number of direct genes annotated to every GO ID (so this duplicated for every GOID for all Genes)
joined_tibble <- joined_tibble %>% group_by(doid_id) %>% mutate(gene_count = n())
#add doid_name to tibble
joined_tibble <- left_join(joined_tibble, terms_tibble, by = "doid_id")
# get gene_names in commas-space separate list
joined_tibble <- joined_tibble %>%
  group_by(doid_id, doid_name) %>%
  summarise_all(funs(paste(sort(.), collapse = ", ")))
# get only one gene count
joined_tibble <- joined_tibble %>%
    mutate(gene_count=sapply(strsplit(as.character(gene_count), ","), function(x) x[[1]][1]))
#ensure correct column order
joined_tibble <- joined_tibble[,c("doid_id", "doid_name", "gene_count", "gene_id")]
#arrange from largest to smallest gene sets
joined_tibble$gene_count <-as.integer(joined_tibble$gene_count)
joined_tibble <- joined_tibble %>% arrange(desc(gene_count))


# #write final file
joined_tibble %>%
  as.data.frame() %>%
  write_delim(paste0(base_dir,"propogated_annotations_doid.tsv"),delim = "\t",col_names = TRUE)
#script time
toc <- as.integer(as.POSIXct(Sys.time()))
print(paste("The time it took to run this script in minutes was", (toc-tic)/60, sep = " "))
