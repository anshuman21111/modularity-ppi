##This code uses PPI network data and calculates spectral entropy; all data stored in dir interactomes

library(readr)
library(DescTools)
library(ggplot2)

path_input = "interactomes/"
path_output = "adjacency_singular_values/"

filenames = list.files(path_input)
tax_ids_all = as.integer(gsub(filenames, pattern=".txt",replacement = ""))

filenames = list.files(path_output)
tax_ids_finished = as.integer(gsub(filenames, pattern=".txt",replacement = ""))

tax_ids = setdiff(tax_ids_all, tax_ids_finished)

for (current_id in tax_ids) {
  print(paste0("current_id = " , as.character(current_id)))
  current_file = paste0(path_input, current_id, ".txt")
  current_save = paste0(path_output,current_id, ".txt")
  
  file_size = file.size(current_file)
  
  if(file_size < 2e6){
    test_interactome = read.table(current_file, sep = " ", stringsAsFactors = FALSE)
    genes = as.character(unique(unlist(test_interactome)))
    N = length(genes)
    adj_matrix = matrix(0, N, N)
    row.names(adj_matrix) = genes
    colnames(adj_matrix) = genes
    
    A = as.character(test_interactome$V1)
    B = as.character(test_interactome$V2)
    adj_matrix[cbind(A,B)] = 1
    adj_matrix[cbind(B,A)] = 1
    
    try({
      tmp = svd(adj_matrix)$d;
      write.csv(tmp, file = current_save, row.names = FALSE)
    })
    #tmp = svd(adj_matrix)$d
    #write.csv(tmp, file = current_save, row.names = FALSE)
  } else{
    print(paste0("FILE TOO BIG: " ,current_file))
    
  }
  
}


#please note that all protein names have to be characters

path_input = "interactomes/"
path_output = "adjacency_singular_values/"

filenames = list.files(path_input)
tax_ids_all = as.integer(gsub(filenames, pattern=".txt",replacement = ""))

filenames = list.files(path_output)
tax_ids_finished = as.integer(gsub(filenames, pattern=".txt",replacement = ""))

tax_ids = setdiff(tax_ids_all, tax_ids_finished)


rm(list = ls())
library(readr)
library(DescTools)

path_matrix = "adjacency_singular_values/"
path_taxonomy = "taxonomy/"

## Functions

entropy <- function (sv, bins, normalizeEntropy = TRUE){
  N = length(sv)
  svn = sv/max(sv)
  
  counts = .bincode(svn,breaks = seq(0,1,1/bins),include.lowest = TRUE)
  p = counts/sum(counts)
  
  entropy = sum(-p*log(p))
  
  if(normalizeEntropy){
    entropy = entropy/log(length(p))
  }
  
  return(entropy)
}

Gini <- function(sv){
  coef = Gini(temp)
  return(coef)
  
}

## Loop


filenames = list.files(path_matrix)
tax_ids = as.integer(gsub(filenames, pattern=".txt",replacement = ""))

df = data.frame(id = NULL, species = NULL, superkingdom = NULL, phyllum = NULL,
                Gini = NULL,
                entropy_N = NULL,
                entropy_sqrtN = NULL,
                entropy_500 = NULL,
                entropy_100 = NULL,
                entropy_50 = NULL)


for (current_id in tax_ids) {
  print(paste0("current_id = " , as.character(current_id)))
  current_file = paste0(path_matrix, current_id, ".txt")
  current_taxonomy_file = paste0(path_taxonomy, current_id, ".txt")
  current_taxonomy = unlist(read_delim(current_taxonomy_file, ";", escape_double = FALSE, trim_ws = TRUE,col_names = FALSE, col_types = cols()    ))
  superkingdom = current_taxonomy[2]
  phyllum = current_taxonomy[3]
  species = (tail(current_taxonomy, n=1))
  
  name = paste0(path_matrix, current_id, ".txt")
  input = read_csv(name, col_types = cols()  )
  sv =input$x
  sv = log2(sv)
  sv[sv<0] = 0
  N = length(sv)
  
  
  stopifnot(min(sv) >= 0)
  
  
  #ginicoef = Gini(sv) ## For some reason, this is throwing an error
  ginicoef = NA
  
  entropy_N = entropy(sv, bins = N)
  entropy_sqrtN = entropy(sv, bins = floor(sqrt(N)))
  entropy_500 = entropy(sv, bins = 500)
  entropy_100 = entropy(sv, bins = 100)
  entropy_50 = entropy(sv, bins = 50)
  
  df_current = data.frame(id = current_id, species = species, superkingdom = superkingdom, phyllum = phyllum,
                          Gini = ginicoef,
                          entropy_N = entropy_N,
                          entropy_sqrtN = entropy_sqrtN,
                          entropy_500 = entropy_500,
                          entropy_100 = entropy_100,
                          entropy_50 = entropy_50)
  df = rbind(df, df_current)
  
  
}

#write.csv(df, file = "df_adj_transformed.csv", row.names = FALSE)

