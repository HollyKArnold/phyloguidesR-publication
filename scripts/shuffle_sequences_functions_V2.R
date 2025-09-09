# write_shuffles_mix
# Randomly mix of sequences from two different fastas for mixed simulations. 
# in_fasta_list: a vector of fasta file paths.
# N the number of permutations to output mixed trees
# out_dir: the output directory
write_shuffles_mix = function(in_fasta_list, N, out_fasta_prefix, out_dir, guides = FALSE){
  
  # Checks - list should have at least two elements
  if(length(in_fasta_list) <= 1){stop("FAIL: Must provide more than two fastas for a mix.")}
  
  # Check that all the fastas exist. 
  for(f in 1:length(in_fasta_list)){
    # get the current fasta file
    cur_fasta = in_fasta_list[f]
    
    # Check file exists
    if(!file.exists(cur_fasta)){stop("FAIL: File does not exist at that path")}
  }
  
  if(guides){
    
    
    # Read in the first file
    cur_name = names(in_fasta_list)[1]
    mix = as_tibble(as.data.frame(read.fasta(in_fasta_list[1]))) %>%
      rename("ID" = seq.name) %>%
      rename({{cur_name}} := seq.text)
    
    #Keep the ref seqs for later
    mix_ref = mix %>%
      filter(stringr::str_detect(string = ID, pattern = "REF")) 
    mix_ref = as.data.frame(mix_ref)
    colnames(mix_ref) = c("ID", "SEQ")
    
    # Now, filter out the ref sequences for the mix matrix
    mix = 
      mix %>%
      filter(stringr::str_detect(string = ID, pattern = "REF", negate = TRUE))
    
    
    for(f in 2:length(in_fasta_list)){
      cur_name = names(in_fasta_list)[f]
      cur_fasta = as_tibble(as.data.frame(read.fasta(in_fasta_list[f]))) %>%
        rename("ID" = seq.name) %>%
        rename({{cur_name}} := seq.text) %>%
        filter(stringr::str_detect(string = "ID", pattern = "REF", negate = TRUE))
      
      mix = 
        mix %>%
        left_join(cur_fasta, join_by(ID))
    }
    
    N_mix = ncol(mix) - 1
    
    for(i in 1:N){
      
      # Get the random mix of columns.
      choose_HVR = sample(x = seq(from = 1, to = N_mix, by = 1), size = nrow(mix), replace = TRUE, prob = rep(x = 1/N_mix, length = N_mix)) + 1
      
      #Format the data matrix for the mixing.
      cur_shuffle = as.data.frame(mix)
      cur_shuffle$choose_HVR = choose_HVR
      cur_shuffle$SEQ = NA
      
      for(j in 1:nrow(cur_shuffle)){
        cur_shuffle[j,"SEQ"] = cur_shuffle[j, as.vector(cur_shuffle[j, "choose_HVR"])]
      }
      cur_shuffle = rbind(cur_shuffle[,c("ID", "SEQ")],
                          mix_ref[,c("ID", "SEQ")])
      
      # Mix the cur suffled rows.
      cur_shuffle = cur_shuffle[sample(x = nrow(cur_shuffle), replace = FALSE), ]
      
      # Write the output.
      filename = paste0(c(out_fasta_prefix, i, "_shuffle.align"), collapse = "", sep = "")
      filename = file.path(out_dir, filename)
      write_fasta(df = cur_shuffle, file = filename)
      
    }
    
    
    
  }else{
    
    # Read in the first file
    cur_name = names(in_fasta_list)[1]
    mix = as_tibble(as.data.frame(read.fasta(in_fasta_list[1]))) %>%
      rename("ID" = seq.name) %>%
      rename({{cur_name}} := seq.text)
    
    for(f in 2:length(in_fasta_list)){
      cur_name = names(in_fasta_list)[f]
      cur_fasta = as_tibble(as.data.frame(read.fasta(in_fasta_list[f]))) %>%
        rename("ID" = seq.name) %>%
        rename({{cur_name}} := seq.text)
      
      mix = 
        mix %>%
        left_join(cur_fasta, join_by(ID))
    }
    
    N_mix = ncol(mix) - 1
    
    for(i in 1:N){
      
      # Get the random mix of columns.
      choose_HVR = sample(x = seq(from = 1, to = N_mix, by = 1), size = nrow(mix), replace = TRUE, prob = rep(x = 1/N_mix, length = N_mix)) + 1
      
      #Format the data matrix for the mixing.
      cur_shuffle = as.data.frame(mix)
      cur_shuffle$choose_HVR = choose_HVR
      cur_shuffle$SEQ = NA
      
      for(j in 1:nrow(cur_shuffle)){
        cur_shuffle[j,"SEQ"] = cur_shuffle[j, as.vector(cur_shuffle[j, "choose_HVR"])]
      }
      cur_shuffle = cur_shuffle[,c("ID", "SEQ")]
      
      # Mix the cur suffle rows.
      cur_shuffle = cur_shuffle[sample(x = nrow(cur_shuffle), replace = FALSE), ]
      
      # Write the output.
      filename = paste0(c(out_fasta_prefix, i, "_shuffle.align"), collapse = "", sep = "")
      filename = file.path(out_dir, filename)
      write_fasta(df = cur_shuffle, file = filename)
      
    }
    
  }
}

# Randomly choose a sequence from N different fastas which have no guides. 
# in_fasta_list: a vector of fasta file paths.
# N the number of permutations to output mixed trees
# out_dir: the output directory
write_shuffles_mix_no_guides = function(in_fasta_list, N, out_fasta_prefix, out_dir){
  
  # Checks - list should have at least two elements
  if(length(in_fasta_list) <= 1){stop("FAIL: Must provide more than two fastas for a mix.")}
  
  # Check that all the fastas exist. 
  for(f in 1:length(in_fasta_list)){
    # get the current fasta file
    cur_fasta = in_fasta_list[f]
    
    # Check file exists
    if(!file.exists(cur_fasta)){stop("FAIL: File does not exist at that path")}
  }
  
  
  # Read in the first file
  cur_name = names(in_fasta_list)[1]
  mix = as_tibble(as.data.frame(read.fasta(in_fasta_list[1]))) %>%
    rename("ID" = seq.name) %>%
    rename({{cur_name}} := seq.text)
  
  for(f in 2:length(in_fasta_list)){
    cur_name = names(in_fasta_list)[f]
    cur_fasta = as_tibble(as.data.frame(read.fasta(in_fasta_list[f]))) %>%
      rename("ID" = seq.name) %>%
      rename({{cur_name}} := seq.text)
    
    mix = 
      mix %>%
      left_join(cur_fasta, join_by(ID))
  }
  
  N_mix = ncol(mix) - 1
  
  for(i in 1:N){
    
    # Get the random mix of columns.
    choose_HVR = sample(x = seq(from = 1, to = N_mix, by = 1), size = nrow(mix), replace = TRUE, prob = rep(x = 1/N_mix, length = N_mix)) + 1
    
    #Format the data matrix for the mixing.
    cur_shuffle = as.data.frame(mix)
    cur_shuffle$choose_HVR = choose_HVR
    cur_shuffle$SEQ = NA
    
    for(j in 1:nrow(cur_shuffle)){
      cur_shuffle[j,"SEQ"] = cur_shuffle[j, as.vector(cur_shuffle[j, "choose_HVR"])]
    }
    cur_shuffle = cur_shuffle[,c("ID", "SEQ")]
    
    # Mix the cur suffled rows.
    cur_shuffle = cur_shuffle[sample(x = nrow(cur_shuffle), replace = FALSE), ]
    
    # Write the output.
    filename = paste0(c(out_fasta_prefix, i, "_shuffle.fasta"), collapse = "", sep = "")
    filename = file.path(out_dir, filename)
    write_fasta(df = cur_shuffle, file = filename)
    
  }
  
  
}

# Writes a shuffled version of an input fasta
# in_fasta: Input fasta sequence file to shuffle
# N: the number of shuffles
# out_fasta_prefix: Output fasta file prefix
# out_dir: output directory for which to files out to.
write_shuffles = function(in_fasta, N, out_fasta_prefix, out_dir){
  if(file.exists(in_fasta)){
    in_fasta = as.data.frame(read.fasta(in_fasta))
    colnames(in_fasta) = c("ID", "SEQ")
  }else{
    stop("Your input fasta file does not exist. ")
  }
  
  for(i in 1:N){
    cur_shuffle = in_fasta[sample(x = nrow(in_fasta), replace = FALSE), ]
    filename = paste0(c(out_fasta_prefix, i, "_shuffle.align"), collapse = "", sep = "")
    filename = file.path(out_dir, filename)
    write_fasta(df = cur_shuffle, file = filename)
  }
  
  
}

# Writes out a fasta file. 
# df: This is a dataframe where colnames are ID for the identifier and SEQ for the sequence.
# File the file path to write out to and filename.
write_fasta = function(df, file) {
  
  sink(file)
  for (i in 1:nrow(df)) {
    cur = paste0(c(">", as.vector(df[i, "ID"]), "\n", 
                   as.character(df[i, "SEQ"]), "\n"), sep = "", 
                 collapse = "")
    cat(cur)
  }
  closeAllConnections()
  
}
