# AUTHOR: ARNOLD
# Functions for maketree_sequences_V2.R

# Write commands for making trees to output to file
make_fasttrees = function(in_dir, N, file_prefix, out_dir, path.fasttree = "/local/cqls/opt/x86_64/bin/FastTree"){
  cmd = vector()
  
  # Check that all files exist
  all_exist = TRUE
  for(i in 1:N){
    
    filename = paste0(c(file_prefix, i, "_shuffle.align"), sep = "", collapse = "")
    filename = file.path(in_dir, filename)
    if(!file.exists(filename)){
      
      stop(paste0(c("FAIL: File doesn't exist", filename), sep = "", collapse = ""))
    }
  }
  if(all_exist){
    print("PASS: All files exist")
  }  
  
  # Now, we will run fasttree for each
  for(i in 1:N){
    filename = paste0(c(file_prefix, i, "_shuffle.align"), sep = "", collapse = "")
    filename = file.path(in_dir, filename)
    cmd = c(cmd, make_fasttree(path.to.alignment.file = filename, path.fasttree = path.fasttree, output.file = paste0(c(filename, ".tre"), sep = "", collapse = "")))
  }
  return(cmd)
}


# Make the fast tree command itself
make_fasttree = function(path.to.alignment.file,
                         fasttree.flags = "-nt -gtr -gamma",
                         path.fasttree,
                         output.file){
  
  cmd = paste0(c("hpcman queue submit \'",
                 path.fasttree,
                 " ",
                 fasttree.flags,
                 " ",
                 path.to.alignment.file,
                 " > ",
                 output.file,
                 "\' -r tree_sim_",
                 basename(output.file), 
                 " -p 1 -q sharpton"), sep = "", collapse = "")
  
  #my_cat(cmd)
  #system(cmd)
  
  return(cmd)
}

# Modification to increase the amount of memory for building larger fasttrees. 
make_fasttree_large = function(path.to.alignment.file,
                         fasttree.flags = "-nt -gtr -gamma",
                         path.fasttree,
                         output.file){
  
  cmd = paste0(c("hpcman queue submit \'",
                 path.fasttree,
                 " ",
                 fasttree.flags,
                 " ",
                 path.to.alignment.file,
                 " > ",
                 output.file,
                 "\' -r tree_sim_",
                 basename(output.file), 
                 " -p 1 -q sharpton -m 32G"), sep = "", collapse = "")
  
  #my_cat(cmd)
  #system(cmd)
  
  return(cmd)
}

# Modification to command to increase the amount of memory for building larger fasttrees. 
make_fasttrees_large = function(in_dir, N, file_prefix, out_dir, path.fasttree = "/local/cqls/opt/x86_64/bin/FastTree"){
  cmd = vector()
  
  # Check that all files exist
  all_exist = TRUE
  for(i in 1:N){
    # recreate filename
    filename = paste0(c(file_prefix, i, "_shuffle.align"), sep = "", collapse = "")
    filename = file.path(in_dir, filename)
    if(!file.exists(filename)){
      
      stop(paste0(c("FAIL: File doesn't exist", filename), sep = "", collapse = ""))
    }
  }
  if(all_exist){
    print("PASS: All files exist")
  }  
  
  # Now, we will run fasttree for each
  for(i in 1:N){
    filename = paste0(c(file_prefix, i, "_shuffle.align"), sep = "", collapse = "")
    filename = file.path(in_dir, filename)
    cmd = c(cmd, make_fasttree_large(path.to.alignment.file = filename, path.fasttree = path.fasttree, output.file = paste0(c(filename, ".tre"), sep = "", collapse = "")))
  }
  return(cmd)
}

# Make fasttree command file and use ruby for building
make_fasttree_ruby = function(path.to.alignment.file,
                               fasttree.flags = "-nt -gtr -gamma",
                               path.fasttree,
                               output.file){
  
  cmd = paste0(c("hpcman queue submit \'",
                 path.fasttree,
                 " ",
                 fasttree.flags,
                 " ",
                 path.to.alignment.file,
                 " > ",
                 output.file,
                 "\' -r tree_sim_",
                 basename(output.file), 
                 " -p 1 -q biomed"), sep = "", collapse = "")
  
  #my_cat(cmd)
  #system(cmd)
  
  return(cmd)
}

# Make fasttree command and use ruby for building
make_fasttrees_ruby = function(in_dir, N, file_prefix, out_dir, path.fasttree = "/local/cqls/opt/x86_64/bin/FastTree"){
  cmd = vector()
  
  # Check that all files exist
  all_exist = TRUE
  for(i in 1:N){
    # recreate filename
    filename = paste0(c(file_prefix, i, "_shuffle.align"), sep = "", collapse = "")
    filename = file.path(in_dir, filename)
    if(!file.exists(filename)){
      
      stop(paste0(c("FAIL: File doesn't exist", filename), sep = "", collapse = ""))
    }
  }
  if(all_exist){
    print("PASS: All files exist")
  }  
  
  # Now, we will run fasttree for each
  for(i in 1:N){
    filename = paste0(c(file_prefix, i, "_shuffle.align"), sep = "", collapse = "")
    filename = file.path(in_dir, filename)
    cmd = c(cmd, make_fasttree_ruby(path.to.alignment.file = filename, path.fasttree = path.fasttree, output.file = paste0(c(filename, ".tre"), sep = "", collapse = "")))
  }
  return(cmd)
}
