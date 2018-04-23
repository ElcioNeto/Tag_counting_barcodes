cat("\n")
cat("Program version 1.0   Access at https://github.com/kdanielmorais/Tag_counting_barcodes ")
cat("\n")
cat("\n")

##Checking packages
list.of.packages = c("optparse", "stringi")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages,  repos = "http://cran.us.r-project.org")

library("optparse")

option_list1 = list(
  make_option(c("-f", "--fast"), type="character",  
        help="A joined .fasta file [Required]", metavar="character"),
  make_option(c("-b", "--bar"), type="character", 
        help="A barcode file [Required]", metavar="character"),
  make_option(c("-t", "--tag"), type="character", default= "FALSE", action= "store",
        help="write TRUE if you want to see the header of the sequence with a tag pair [default= %default] "),
  make_option(c("-n", "--name"), type="character", default="FALSE", action= "store",
        help="write TRUE if you want the final sequences named and the barcode tags kept in the sequences, if FALSE, barcode tags will be trimmed off [default= %default] ", metavar="character")
  )

opt_parser = OptionParser(option_list=option_list1);
opt = parse_args(opt_parser)

if (is.null(opt$fast)){
  print_help(opt_parser)
  stop("At least two input files must be supplied. \n--fast [sequences.fasta] and --bar [barcodes.txt] )\n", call.=FALSE)
}

g1 <- opt$tag
g2 <- opt$name
a1 <- opt$fast
a2 <- opt$bar

cat(timestamp())
cat("\n")
#Entering data
#Get sequences fasta and barcodes file
fasta_df = read.delim(a1, header = F, stringsAsFactors = F) # took 30 seconds
barcodes = read.delim(a2, header = F, stringsAsFactors = F)
numb_forward_bc = unique(barcodes[,2])
numb_reverse_bc = unique(barcodes[,6])

#Generating reverse complement barcodes
revcomp = function(DNAstr) {
  step1 = chartr("ACGT","TGCA",DNAstr)
  step2 = unlist(strsplit(step1, split=""))
  step3 = rev(step2)
  step4 = paste(step3, collapse="")
  return(step4)
}
rev_comp_bc=sapply(barcodes[,7], revcomp)

#Generating base for search
df_base_fasta = data.frame(header=fasta_df[seq(1,nrow(fasta_df),2),], bases=as.character(fasta_df[seq(2,nrow(fasta_df),2),]), stringsAsFactors = F)
for (j in numb_forward_bc){
  df_base_fasta[, ncol(df_base_fasta) +1] <- substr(df_base_fasta[,2], 1, j)
  names(df_base_fasta)[ncol(df_base_fasta)] <- paste0("fwd_bc_length_", j)
}
fwd_cols = ncol(df_base_fasta)

##finding fwd and rev barcodes and appending names to the headers
#it would be good to check if multiple primers were found in the same seq. Check if any line of "hit" objects have more than one number.

# finding fwd barc first.
hit_fwd_bc = sapply(df_base_fasta[3:fwd_cols], barcodes[,3], FUN = match) 

cat("\n")

#Changing orientation of reads and checking fwd primer issues
no_fwd_hits = rowSums(!is.na(hit_fwd_bc))
test_if_more_than_one_fwd_barcode_hit_a_sequence = which(no_fwd_hits >1)

#checking if multiple forward barcodes are hitting the same sample
if(length(test_if_more_than_one_fwd_barcode_hit_a_sequence)>0) {
  cat("error!\nMultiple barcodes hitting the same sequence.")
} else { cat("Forward barcodes are fine")
}
position_of_no_fwd_hits = which(no_fwd_hits == 0)

df_base_fasta[position_of_no_fwd_hits,2] <- sapply(df_base_fasta[position_of_no_fwd_hits,2], revcomp)

#Recreating the df_base_fasta with the reversed possible barcodes
df_base_fasta = df_base_fasta[,1:2]

for (j in numb_forward_bc){
  df_base_fasta[, ncol(df_base_fasta) +1] <- substr(df_base_fasta[,2], 1, j)
  names(df_base_fasta)[ncol(df_base_fasta)] <- paste0("fwd_bc_length_", j)
} 

for (j in (numb_reverse_bc)-1){
  df_base_fasta[, ncol(df_base_fasta) +1] <- substr(df_base_fasta[,2], nchar(df_base_fasta[,2])-j, nchar(df_base_fasta[,2]))
  names(df_base_fasta)[ncol(df_base_fasta)] <- paste0("rev_bc_length_", (j)+1)
}
rev_cols = ncol(df_base_fasta)

# finding fwd_barc again but having all seqs in the same orientation
hit_fwd_bc = sapply(df_base_fasta[3:fwd_cols], barcodes[,3], FUN = match)

#appending read names with fwd tag ID
for (i in (1:nrow(barcodes))){
  for ( j in 1:length(numb_forward_bc)) {
    df_base_fasta[c(which(hit_fwd_bc[,j] == i)),1] <- paste(df_base_fasta[c(which(hit_fwd_bc[,j] == i)),1], barcodes[i,4], sep = "|")
  }
}

#appending read names with rev tag ID
hits_rev_bc = as.matrix(sapply(df_base_fasta[,((fwd_cols)+1):rev_cols], rev_comp_bc, FUN = match)) #finding rev_barc, having all seqs in the same frame (orientation)

for (i in (1:nrow(barcodes))){
  for (j in (1:length(numb_reverse_bc))) {
    df_base_fasta[c(which(hits_rev_bc[,j] == i)),1] <- paste(df_base_fasta[c(which(hits_rev_bc[,j] == i)),1], barcodes[i,8], sep = "|")
  }
}

#Filtering sequences with tags
isolating = strsplit(df_base_fasta[,1], "\\|")
library(stringi)
isolating2 = stri_list2matrix(isolating)
headers_with2_tags = which(!is.na(isolating2)[3,]) #getting only the sequences that have two tags. 
#isolating2[,which(!is.na(isolating2)[3,])]
two_tags_headers_and_seq = df_base_fasta[headers_with2_tags,c(1,2)] 

if (g1 == "TRUE") {
  taged_seqs = as.vector(t(two_tags_headers_and_seq))
  write(taged_seqs,"tagged_fasta.fasta")
  cat("\n")
  }

#Counting tags pairs
isolated_tags = isolating2[,headers_with2_tags]
transposed_tags = t(isolated_tags[2:3,])
all_tags = apply(transposed_tags, 1, paste, collapse = "|")
counted_tags = as.data.frame(table(all_tags))
sorted_counted_tags = counted_tags[order(counted_tags[,2],decreasing = T),]

#Changing headers
sample_name_bc = data.frame(samples = barcodes[,1], barcodes = apply(barcodes[,c(4,8)], 1, paste, collapse = "|"), fwd_size = barcodes[,2], rev_size = barcodes[,6], stringsAsFactors = F)
sample_ids_with_tag_position = sample_name_bc[match(all_tags, sample_name_bc$barcodes),c(1,3,4)]
headers = c(isolated_tags[1,])
renaming_tags = data.frame(headers, sample_ids_with_tag_position, stringsAsFactors = F)
appending_seq = cbind(renaming_tags, two_tags_headers_and_seq)
removed_NA_id_seqs = appending_seq[which(!is.na(appending_seq[,2])),]

name = sample_name_bc[match( sorted_counted_tags[,1], sample_name_bc[,2]),1]
sorted_counted_tags$tags_name=name

#Write the counts of each tag pair and the sample names for each tag pair, based on the barcodes file.
write.table(sorted_counted_tags, "read_counts2.tsv", row.names = F, sep = "\t", quote = F)

#write file sequences with name but still with tags
final_renamed_sequences_with_tags_df = cbind(apply(removed_NA_id_seqs[,c(1,2)], 1,paste, collapse = "|"), removed_NA_id_seqs[,6])

if (g2 == "TRUE"){
  named_seqs_with_tags = as.vector(t(final_renamed_sequences_with_tags_df))
  write(named_seqs_with_tags,"renamed_fasta_with_tags.fasta")
  cat("\n")
  }

#write sequences with name and without tags
final_renamed_sequences_no_tags_df = cbind(final_renamed_sequences_with_tags_df[,1], substr(final_renamed_sequences_with_tags_df[,2], (removed_NA_id_seqs$fwd_size)+1, (nchar(removed_NA_id_seqs$bases))-(removed_NA_id_seqs$rev_size)))
named_seqs_without_tags = as.vector(t(final_renamed_sequences_no_tags_df))

write(named_seqs_without_tags,"renamed_fasta_without_tags.fasta")

cat(timestamp())
cat("\n")
