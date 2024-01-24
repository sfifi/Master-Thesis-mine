library(rentrez)
library(seqinr)
library(DESeq2)
library(EnsDb.Hsapiens.v86)
library(dplyr)
library(biomaRt)
library(ggplot2)
library(readxl)


############################## Data 


# Connect to the Ensembl database
mart <- useMart("ensembl", dataset="hsapiens_gene_ensembl")

# unique_genes_RNA is a list of my genes that I created in another R script

unique_genes_RNA_S<-sample(unique_genes_RNA,50)

# Fetch protein sequences
sequences <- biomaRt::getSequence(id=unique_genes_RNA_S, type="hgnc_symbol", seqType="peptide", mart=mart)

# Filter sequences based on length
sequences <- sequences[nchar(sequences$peptide) > 100,]
sequences<-sequences$peptide
# Define the Kyte-Doolittle scale
kd_scale <- c(A=1.8, R=-4.5, N=-3.5, D=-3.5, C=2.5, Q=-3.5, E=-3.5, G=-0.4, H=-3.2, I=4.5, 
              L=3.8, K=-3.9, M=1.9, F=2.8, P=-1.6, S=-0.8, T=-0.7, W=-0.9, Y=-1.3, V=4.2)

# Sliding window function
compute_metahydropathy <- function(seq, window_size=19) {
  half_window <- floor(window_size / 2)
  len <- nchar(seq)
  values <- numeric(len)
  
  # For each position in the sequence
  for (i in 1:len) {
    start <- max(1, i - half_window)
    end <- min(len, i + half_window)
    sub_seq <- substr(seq, start, end)
    sub_values <- kd_scale[unlist(strsplit(sub_seq, NULL))]
    values[i] <- mean(sub_values, na.rm = TRUE)
  }
  
  return(values)
}

# Compute metahydropathy for each protein using a sliding window
metahydropathy_values <- sapply(sequences, compute_metahydropathy)

# Compute the mean metahydropathy value at each position
max_length <- max(sapply(metahydropathy_values, length))
padded_values <- lapply(metahydropathy_values, function(v) {
  c(v, rep(NA, max_length - length(v)))
})
mean_values <- apply(do.call(rbind, padded_values), 2, function(column) {
  mean(column, na.rm = TRUE)
})

# Plotting
data <- data.frame(x = 1:length(mean_values), y = mean_values)

# Replace leading NA values with 0
leading_na_indices <- which(!is.na(data$y))[1] - 1
data$y[1:leading_na_indices] <- 0

# Split the data into positive and negative segments
data_positive <- data
data_positive$y[data_positive$y < 0] <- 0

data_negative <- data
data_negative$y[data_negative$y > 0] <- 0

# Plot
p <- ggplot() +
  geom_line(data=data, aes(x=x, y=y)) +
  geom_hline(yintercept = 0, linetype="dashed", color = "red") +
  geom_area(data=data_positive, aes(x=x, y=y), fill="blue") +
  geom_area(data=data_negative, aes(x=x, y=y), fill="pink",alpha=0.2) +
  labs(title="Average Metahydropathy Plot", x="Residue Position", y="Metahydropathy") +
  theme_minimal() + ylim(-3,2)

print(p)



############################## random protein

library(seqinr)
library(ggplot2)
library(rentrez)

# Search UniProt for human proteins
search_results <- entrez_search(db="protein", term="human[orgn]", retmax=5000)

# Randomly select 50 protein IDs
#set.seed(123) # for reproducibility
selected_ids <- sample(search_results$ids, 50)

# Fetch the protein sequences
sequences <- entrez_fetch(db="protein", id=selected_ids, rettype="fasta", retmode="text")
fasta_list <- read.fasta(textConnection(sequences), as.string = TRUE)

# Filter sequences with length > 100 and < 600 aa
filtered_proteins <- fasta_list[sapply(fasta_list, nchar) > 600 & sapply(fasta_list, nchar) < 1000]

# Convert each sequence from string to character vector of amino acids
filtered_proteins_converted <- lapply(filtered_proteins, function(seq) {
  return(s2c(as.character(seq)))
})


#### just for elastin 

# Search UniProt for human proteins with the gene symbol ELN
search_results <- entrez_search(db="protein", term="AGO2[Gene Name] AND human[orgn]", retmax=5000)

selected_ids <- sample(search_results$ids, 5)

sequences <- entrez_fetch(db="protein", id=selected_ids, rettype="fasta", retmode="text")

fasta_list <- read.fasta(textConnection(sequences), as.string = TRUE)
# Filter sequences with length > 100 and < 600 aa
filtered_proteins <- fasta_list[sapply(fasta_list, nchar) > 200]


filtered_proteins_converted <- lapply(filtered_proteins, function(seq) {
  return(s2c(as.character(seq)))
})

##### end just for elastin #######

avg_hydropathy <- function(sequence, window_size = 19) {
  kyte_doolittle <- c(a=1.8, r=-4.5, n=-3.5, d=-3.5, c=2.5, q=-3.5, e=-3.5, g=-0.4, 
                      h=-3.2, i=4.5, l=3.8, k=-3.9, m=1.9, f=2.8, p=-1.6, s=-0.8, 
                      t=-0.7, w=-0.9, y=-1.3, v=4.2)
  
  sequence_values <- kyte_doolittle[sequence]
  seq_length <- length(sequence_values)
  
  metahydropathy <- sapply(1:(seq_length - window_size + 1), function(i) {
    mean(sequence_values[i:(i + window_size - 1)])
  })
  
  return(metahydropathy)
}

# Compute metahydropathy values for all long proteins
metahydropathy_values <- lapply(filtered_proteins_converted, avg_hydropathy)

# Find the maximum length among all metahydropathy values
max_length <- max(sapply(metahydropathy_values, length))

# Pad the values to make them of equal length
padded_values <- lapply(metahydropathy_values, function(v) {
  c(v, rep(NA, max_length - length(v)))
})

# Compute the mean metahydropathy value at each position
mean_values <- apply(do.call(rbind, padded_values), 2, mean, na.rm = TRUE)

# Plotting
data <- data.frame(x = 1:length(mean_values), y = mean_values)

# Modify the data to split it into positive and negative segments
data$y_positive <- ifelse(data$y > 0, data$y, 0)
data$y_negative <- ifelse(data$y < 0, data$y, 0)

# Plotting without annotations
p <- ggplot(data, aes(x=x)) +
  geom_line(aes(y=y)) +
  geom_hline(yintercept = 0, linetype="dashed", color = "red") +
  geom_area(aes(y=y_positive), fill="blue") +
  geom_area(aes(y=y_negative), fill="pink",alpha=0.2) +
  labs(title="Average Metahydropathy Plot", x="Residue Position", y="Metahydropathy") +
  theme_minimal() + ylim(-3,1) +
  theme(axis.title = element_text(size = 14),   # Adjust the size of axis titles
        axis.text = element_text(size = 12))

print(p)




########### looking for SP in radom ###########



library(seqinr)
library(ggplot2)
library(rentrez)

# Search UniProt for human proteins
search_results <- entrez_search(db="protein", term="human[orgn]", retmax=5000)

# Randomly select 200 protein IDs
selected_ids <- sample(search_results$ids, 100)

# Fetch the protein sequences
sequences <- entrez_fetch(db="protein", id=selected_ids, rettype="fasta", retmode="text")
fasta_list <- read.fasta(textConnection(sequences), as.string = TRUE)

# Filter sequences with length > 100
filtered_proteins <- fasta_list[sapply(fasta_list, nchar) > 100]

# Convert each sequence from string to character vector of amino acids
converted_proteins <- lapply(filtered_proteins, function(seq) {
  return(s2c(as.character(seq)))
})

# Truncate each sequence to ensure at least 50 metahydropathy values
window_size = 19
truncated_proteins <- lapply(converted_proteins, function(seq) {
  return(seq[1:min(100 + window_size - 1, length(seq))])
})

# Compute metahydropathy values function
avg_hydropathy <- function(sequence, window_size = 9) {
  kyte_doolittle <- c(a=1.8, r=-4.5, n=-3.5, d=-3.5, c=2.5, q=-3.5, e=-3.5, g=-0.4, 
                      h=-3.2, i=4.5, l=3.8, k=-3.9, m=1.9, f=2.8, p=-1.6, s=-0.8, 
                      t=-0.7, w=-0.9, y=-1.3, v=4.2)
  
  sequence_values <- kyte_doolittle[sequence]
  seq_length <- length(sequence_values)
  
  metahydropathy <- sapply(1:(seq_length - window_size + 1), function(i) {
    mean(sequence_values[i:(i + window_size - 1)])
  })
  
  return(metahydropathy)
}

# Compute metahydropathy values for all truncated proteins
metahydropathy_values <- lapply(truncated_proteins, avg_hydropathy)

# Ensure all proteins are of the same length
max_length <- max(sapply(metahydropathy_values, length))
padded_values <- lapply(metahydropathy_values, function(v) {
  c(v, rep(NA, max_length - length(v)))
})

# Compute the mean metahydropathy value at each position
mean_values <- apply(do.call(rbind, padded_values), 2, mean, na.rm = TRUE)

# Plotting
data <- data.frame(x = 1:length(mean_values), y = mean_values)
data$y_positive <- ifelse(data$y > 0, data$y, 0)
data$y_negative <- ifelse(data$y < 0, data$y, 0)

# Plot
p <- ggplot(data, aes(x=x)) +
  geom_line(aes(y=y)) +
  geom_hline(yintercept = 0, linetype="dashed", color = "red") +
  geom_area(aes(y=y_positive), fill="blue") +
  geom_area(aes(y=y_negative), fill="pink",alpha=0.2) +
  labs(title="Average Metahydropathy Plot", x="Residue Position", y="Metahydropathy") +
  theme_minimal() + ylim(-1,0.5) + xlim(0,100)

print(p)





########## looking for SP in mine ############

# Assuming you have your protein sequences stored in a variable named 'sequences'
# If not, you might need to read or fetch them first.

# Connect to the Ensembl database
mart <- biomaRt::useMart("ensembl", dataset="hsapiens_gene_ensembl")

unique_genes_RNA_S<-sample(unique_genes_RNA,19)

# Fetch protein sequences
sequences <- biomaRt::getSequence(id=unique_genes_RNA_S, type="hgnc_symbol", seqType="peptide", mart=mart)

# Filter sequences based on length
sequences <- sequences[nchar(sequences$peptide) > 1000,]

sequences<-sequences$peptide

# Convert each sequence from string to character vector of amino acids
sequences_p <- lapply(sequences, function(seq) {
  return(s2c(as.character(seq)))
})

# Truncate each sequence to ensure at least 50 metahydropathy values
window_size = 19
sequences_p <- lapply(sequences_p, function(seq) {
  return(seq[1:min(50 + window_size - 1, length(seq))])
})

# Compute metahydropathy values function
avg_hydropathy <- function(sequence, window_size = 19) {
  kyte_doolittle <- c(A=1.8, R=-4.5, N=-3.5, D=-3.5, C=2.5, Q=-3.5, E=-3.5, G=-0.4, H=-3.2, I=4.5, 
                      L=3.8, K=-3.9, M=1.9, F=2.8, P=-1.6, S=-0.8, T=-0.7, W=-0.9, Y=-1.3, V=4.2)
  
  sequence_values <- kyte_doolittle[sequence]
  seq_length <- length(sequence_values)
  
  metahydropathy <- sapply(1:(seq_length - window_size + 1), function(i) {
    mean(sequence_values[i:(i + window_size - 1)])
  })
  
  return(metahydropathy)
}

# Compute metahydropathy for each protein using a sliding window
metahydropathy_values <- lapply(sequences_p, avg_hydropathy)

# Compute the mean metahydropathy value at each position
max_length <- max(sapply(metahydropathy_values, length))
padded_values <- lapply(metahydropathy_values, function(v) {
  c(v, rep(NA, max_length - length(v)))
})
mean_values <- apply(do.call(rbind, padded_values), 2, function(column) {
  mean(column, na.rm = TRUE)
})

# Plotting
data <- data.frame(x = 1:length(mean_values), y = mean_values)

# Replace leading NA values with 0
leading_na_indices <- which(!is.na(data$y))[1] - 1
data$y[1:leading_na_indices] <- 0

# Split the data into positive and negative segments
data_positive <- data
data_positive$y[data_positive$y < 0] <- 0

data_negative <- data
data_negative$y[data_negative$y > 0] <- 0

# Plot
p <- ggplot() +
  geom_line(data=data, aes(x=x, y=y)) +
  geom_hline(yintercept = 0, linetype="dashed", color = "red") +
  geom_area(data=data_positive, aes(x=x, y=y), fill="blue") +
  geom_area(data=data_negative, aes(x=x, y=y), fill="pink",alpha=0.2) +
  labs(title="Average Metahydropathy Plot", x="Residue Position", y="Metahydropathy") +
  theme_minimal() + ylim(-0.8,0.6)

print(p)

############ From c-terminus #############



# Assuming you have already fetched the sequences and filtered them based on length
# Connect to the Ensembl database
mart <- useMart("ensembl", dataset="hsapiens_gene_ensembl")

unique_genes_RNA_S<-sample(unique_genes_RNA,20)

# Fetch protein sequences
sequences <- biomaRt::getSequence(id=unique_genes_RNA_S, type="hgnc_symbol", seqType="peptide", mart=mart)

# Filter sequences based on length
sequences <- sequences[nchar(sequences$peptide) > 1000,]

# Reverse each sequence
reversed_proteins <- lapply(sequences, rev)

# Truncate each sequence to keep the last 50 amino acids
truncated_proteins <- lapply(reversed_proteins, function(seq) {
  start_index <- max(1, length(seq) - 299)
  return(seq[start_index:length(seq)])
})

# Convert each sequence from string to character vector of amino acids
truncated_proteins <- lapply(truncated_proteins, function(seq) {
  return(s2c(as.character(seq)))
})

# Truncate each sequence to ensure at least 50 metahydropathy values
window_size = 19
truncated_proteins <- lapply(truncated_proteins, function(seq) {
  return(seq[1:min(300 + window_size - 1, length(seq))])
})

# Compute metahydropathy values function
avg_hydropathy <- function(sequence, window_size = 19) {
  kyte_doolittle <- c(A=1.8, R=-4.5, N=-3.5, D=-3.5, C=2.5, Q=-3.5, E=-3.5, G=-0.4, H=-3.2, I=4.5, 
                      L=3.8, K=-3.9, M=1.9, F=2.8, P=-1.6, S=-0.8, T=-0.7, W=-0.9, Y=-1.3, V=4.2)
  
  sequence_values <- kyte_doolittle[sequence]
  seq_length <- length(sequence_values)
  
  metahydropathy <- sapply(1:(seq_length - window_size + 1), function(i) {
    mean(sequence_values[i:(i + window_size - 1)])
  })
  
  return(metahydropathy)
}

# Compute metahydropathy for each protein using a sliding window
metahydropathy_values <- lapply(truncated_proteins , avg_hydropathy)

# Compute the mean metahydropathy value at each position
max_length <- max(sapply(metahydropathy_values, length))
padded_values <- lapply(metahydropathy_values, function(v) {
  c(v, rep(NA, max_length - length(v)))
})
mean_values <- apply(do.call(rbind, padded_values), 2, function(column) {
  mean(column, na.rm = TRUE)
})

# Plotting
data <- data.frame(x = 1:length(mean_values), y = mean_values)

# Replace leading NA values with 0
leading_na_indices <- which(!is.na(data$y))[1] - 1
data$y[1:leading_na_indices] <- 0

# Split the data into positive and negative segments
data_positive <- data
data_positive$y[data_positive$y < 0] <- 0

data_negative <- data
data_negative$y[data_negative$y > 0] <- 0

# Plot
p <- ggplot() +
  geom_line(data=data, aes(x=x, y=y)) +
  geom_hline(yintercept = 0, linetype="dashed", color = "red") +
  geom_area(data=data_positive, aes(x=x, y=y), fill="blue") +
  geom_area(data=data_negative, aes(x=x, y=y), fill="pink",alpha=0.1) +
  labs(title="Average Metahydropathy Plot", x="Residue Position", y="Metahydropathy") +
  theme_minimal() + ylim(-0.8,0.5)

print(p)

