# in_dir
# N: how many file iterations to look for
# out_dir: where to write alginemnt file to.

mothur_align_sim = function(in_dir, N, file_prefix, out_dir, template){
  
  # Check that all files exist
  all_exist = TRUE
  for(i in 1:N){
    # recreate filename
    filename = paste0(c(file_prefix, i, "_shuffle.fasta"), sep = "", collapse = "")
    filename = file.path(in_dir, filename)
    if(!file.exists(filename)){
      stop("FAIL: File doesn't exist")
    }
  }
  if(all_exist){
    print("PASS: All files exist")
  }  
  
  # Now, we will run mothur for each
  for(i in 1:N){
    filename = paste0(c(file_prefix, i, "_shuffle.fasta"), sep = "", collapse = "")
    filename = file.path(in_dir, filename)
    mothur_align_seqs(candidate = filename, output.directory = out_dir, template = template)
  }
}


# candidate: Candidate sequences to align
# output.directory: Output directory for alignment
# template: template seed alignment.
mothur_align_seqs = function(candidate,
                             output.directory,
                             template,
                             search = "kmer",
                             ksize = 8,
                             align = "needleman",
                             match = 1,
                             mismatch = -1,
                             gapopen = -2,
                             gapextend = -1,
                             flip = TRUE,
                             threshold = 0.50,
                             processors = 100){
  if(!file.exists(candidate)){
    stop("Candidate alignment file doesn't exist.")
  }
  if(!file.exists(template)){
    stop("Template alignment file doesn't exist.")
  }
  
  cmd = paste0(c("/local/cqls/opt/x86_64/bin/mothur ",
                 "\"#align.seqs(candidate=",
                 candidate,
                 ", template=",
                 template,
                 ", search=",
                 search,
                 ", ksize=",
                 ksize,
                 ", align=",
                 align,
                 ", match=",
                 match,
                 ", mismatch=",
                 mismatch,
                 ", gapopen=",
                 gapopen,
                 ", gapextend=",
                 gapextend,
                 ", flip=",
                 flip,
                 ", threshold=",
                 threshold,
                 ", processors=",
                 processors,
                 ")\""),
               sep = "", collapse = "")
  
  
  my_cat(cmd)
  print(cmd)
  withr::with_dir(new = output.directory, system(cmd))
  
}
