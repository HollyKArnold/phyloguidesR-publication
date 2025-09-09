
# AUTHOR: Arnold
# DATE 1/10/25
# PURPOSE: Curate silva seed file for purpose of simulations.
SCRIPTS="/nfs3/Sharpton_Lab/prod/prod_restructure/projects/arnoldhk/2025_HVR_Guide_Phylogenetic_Integration/2025_hvr_guide_phylogenetic_integration/scripts/"

#echo "Downloading files..."
#wget -P  /nfs3/Sharpton_Lab/prod/prod_restructure/projects/arnoldhk/2025_HVR_Guide_Phylogenetic_Integration/reference_data/ https://mothur.s3.us-east-2.amazonaws.com/wiki/silva.seed_v138_2.tgz

#echo "Unpacking files..."
#tar zxvf silva.seed_v138_2.tgz

#echo "Getting some stats on the silva seed files..."
#echo "The silva seed 138 alignment file has this many sequences: "
#cat /nfs3/Sharpton_Lab/prod/prod_restructure/projects/arnoldhk/2025_HVR_Guide_Phylogenetic_Integration/reference_data/silva.seed_v138_2.align | grep ">" | wc -l

#echo "The silva seed 138 taxonomy file has this many sequences: "
#cat /nfs3/Sharpton_Lab/prod/prod_restructure/projects/arnoldhk/2025_HVR_Guide_Phylogenetic_Integration/reference_data/silva.seed_v138_2.tax | wc -l

#echo "The silva seed 138 alignment file has this many Eukaroytes:"
#cat /nfs3/Sharpton_Lab/prod/prod_restructure/projects/arnoldhk/2025_HVR_Guide_Phylogenetic_Integration/reference_data/silva.seed_v138_2.align | grep "Eukaryota;" | wc -l

#echo "The silva seed 138 alignment file has this many Bacteria:"
#cat /nfs3/Sharpton_Lab/prod/prod_restructure/projects/arnoldhk/2025_HVR_Guide_Phylogenetic_Integration/reference_data/silva.seed_v138_2.align  | grep "Bacteria;" | wc -l

#echo "The silva seed 138 alignment file has this many Archaea:"
#cat /nfs3/Sharpton_Lab/prod/prod_restructure/projects/arnoldhk/2025_HVR_Guide_Phylogenetic_Integration/reference_data/silva.seed_v138_2.align  | grep "Archaea;" | wc -l

#echo "Removing sequences labeled as eukaryotes from the files..." 
#perl $SCRIPTS/curate_silva_seed_remove_eukaryotes.pl

#echo "Getting some stats on the silva seed file which has eukaryotes removed..."
#echo "The silva seed 138 alignment file has this many sequences: "
#cat /nfs3/Sharpton_Lab/prod/prod_restructure/projects/arnoldhk/2025_HVR_Guide_Phylogenetic_Integration/reference_data/silva.seed_v138_2_Eukaryota_excluded.align | grep ">" | wc -l

#echo "The silva seed 138 alignment file has this many Eukaroytes:"
#cat /nfs3/Sharpton_Lab/prod/prod_restructure/projects/arnoldhk/2025_HVR_Guide_Phylogenetic_Integration/reference_data/silva.seed_v138_2_Eukaryota_excluded.align  | grep "Eukaryota;" | wc -l

#echo "The silva seed 138 alignment file has this many Bacteria:"
#cat /nfs3/Sharpton_Lab/prod/prod_restructure/projects/arnoldhk/2025_HVR_Guide_Phylogenetic_Integration/reference_data/silva.seed_v138_2_Eukaryota_excluded.align  | grep "Bacteria;" | wc -l

#echo "The silva seed 138 alignment file has this many Archaea:"
#cat /nfs3/Sharpton_Lab/prod/prod_restructure/projects/arnoldhk/2025_HVR_Guide_Phylogenetic_Integration/reference_data/silva.seed_v138_2_Eukaryota_excluded.align  | grep "Archaea;" | wc -l
# 
# echo "The number of sequences with the current V1 8F: "
# cat /nfs3/Sharpton_Lab/prod/prod_restructure/projects/arnoldhk/2025_HVR_Guide_Phylogenetic_Integration/reference_data/silva.seed_v138_2_Eukaryota_excluded.align | sed 's/-//g' | sed 's/\.//g' | sed 's/U/T/g' | grep "AGAGTTTGATCATGGCTCA" -B1 | grep ">" | wc -l
# 
# echo "The number of sequences with the current V1 68F: "
# cat /nfs3/Sharpton_Lab/prod/prod_restructure/projects/arnoldhk/2025_HVR_Guide_Phylogenetic_Integration/reference_data/silva.seed_v138_2_Eukaryota_excluded.align | sed 's/-//g' | sed 's/\.//g' | sed 's/U/T/g' | grep "T[ACTG]A[ACTG]ACATGCAAGTCG[AG][AG]CG" -B1 | grep ">" | wc -l                                                                     
# 
# echo "The number of sequences with the current V2: "
# cat /nfs3/Sharpton_Lab/prod/prod_restructure/projects/arnoldhk/2025_HVR_Guide_Phylogenetic_Integration/reference_data/silva.seed_v138_2_Eukaryota_excluded.align | sed 's/-//g' | sed 's/\.//g' | sed 's/U/T/g' | grep "GGCG[ACG]ACGGGTGAGTAA" -B1 | grep ">" | wc -l
# 
# echo "The number of sequences with the current V3: "
# cat /nfs3/Sharpton_Lab/prod/prod_restructure/projects/arnoldhk/2025_HVR_Guide_Phylogenetic_Integration/reference_data/silva.seed_v138_2_Eukaryota_excluded.align | sed 's/-//g' | sed 's/\.//g' | sed 's/U/T/g' | grep "CCTACGGGAGGCAGCAG" -B1 | grep ">" | wc -l
# 
# echo "The number of sequences with the current V4: "
# cat /nfs3/Sharpton_Lab/prod/prod_restructure/projects/arnoldhk/2025_HVR_Guide_Phylogenetic_Integration/reference_data/silva.seed_v138_2_Eukaryota_excluded.align | sed 's/-//g' | sed 's/\.//g' | sed 's/U/T/g' | grep "GTG[CT]CAGC[AC]GCCGCGGTAA" -B1 | grep ">" | wc -l
# 
# echo "The number of sequences with the current V5: "
# cat /nfs3/Sharpton_Lab/prod/prod_restructure/projects/arnoldhk/2025_HVR_Guide_Phylogenetic_Integration/reference_data/silva.seed_v138_2_Eukaryota_excluded.align | sed 's/-//g' | sed 's/\.//g' | sed 's/U/T/g' | grep "ATTAGA[AT]ACCC[CGT][ACTG]GTAGTCC" -B1 | grep ">" | wc -l
# 
# echo "The number of sequences with the current V6: "
# cat /nfs3/Sharpton_Lab/prod/prod_restructure/projects/arnoldhk/2025_HVR_Guide_Phylogenetic_Integration/reference_data/silva.seed_v138_2_Eukaryota_excluded.align | sed 's/-//g' | sed 's/\.//g' | sed 's/U/T/g' | grep "ACTCAAA[GT]GAATTGACGGGG[AG]C" -B1 | grep ">" | wc -l
# 
# echo "The number of sequences with the current V7: "
# cat /nfs3/Sharpton_Lab/prod/prod_restructure/projects/arnoldhk/2025_HVR_Guide_Phylogenetic_Integration/reference_data/silva.seed_v138_2_Eukaryota_excluded.align | sed 's/-//g' | sed 's/\.//g' | sed 's/U/T/g' | grep "GTG[CG]TGCATGG[CT]TGTCGTCAGCT" -B1 | grep ">" | wc -l
# 
# echo "The number of sequences with the current V8: "
# cat /nfs3/Sharpton_Lab/prod/prod_restructure/projects/arnoldhk/2025_HVR_Guide_Phylogenetic_Integration/reference_data/silva.seed_v138_2_Eukaryota_excluded.align | sed 's/-//g' | sed 's/\.//g' | sed 's/U/T/g' | grep "GGAAGG[CT]GGGGA[CT]GACGTCAA" -B1 | grep ">" | wc -l
# 
# echo "The number of sequences with the current V9: "
# cat /nfs3/Sharpton_Lab/prod/prod_restructure/projects/arnoldhk/2025_HVR_Guide_Phylogenetic_Integration/reference_data/silva.seed_v138_2_Eukaryota_excluded.align | sed 's/-//g' | sed 's/\.//g' | sed 's/U/T/g' | grep "TG[CT]AC[AT]CACCGCCCGTC" -B1 | grep ">" | wc -l
# 
# echo "The number of sequences with the current VEnd: "
# cat /nfs3/Sharpton_Lab/prod/prod_restructure/projects/arnoldhk/2025_HVR_Guide_Phylogenetic_Integration/reference_data/silva.seed_v138_2_Eukaryota_excluded.align | sed 's/-//g' | sed 's/\.//g' | sed 's/U/T/g' | grep "AAGTCGTAACAAGGTAACCGTA" -B1 | grep ">" | wc -l
# 
# echo "Finding sequences which have all the hyperconserved regions..."
# cat /nfs3/Sharpton_Lab/prod/prod_restructure/projects/arnoldhk/2025_HVR_Guide_Phylogenetic_Integration/reference_data/silva.seed_v138_2_Eukaryota_excluded.align | sed 's/-//g' | sed 's/\.//g' | sed 's/U/T/g' | grep "T[ACTG]A[ACTG]ACATGCAAGTCG[AG][AG]CG" -B1 | grep "GGCG[ACG]ACGGGTGAGTAA" -B1 | grep "CCTACGGGAGGCAGCAG" -B1 | grep "GTG[CT]CAGC[AC]GCCGCGGTAA" -B1 | grep "ATTAGA[AT]ACCC[CGT][ACTG]GTAGTCC" -B1 | grep "ACTCAAA[GT]GAATTGACGGGG[AG]C" -B1 | grep "GTG[CG]TGCATGG[CT]TGTCGTCAGCT" -B1 | grep "GGAAGG[CT]GGGGA[CT]GACGTCAA" -B1 | grep "TG[CT]AC[AT]CACCGCCCGTC" -B1 > /nfs3/Sharpton_Lab/prod/prod_restructure/projects/arnoldhk/2025_HVR_Guide_Phylogenetic_Integration/reference_data/silva.seed_v138_2_Eukaryota_excluded_all_hyperconserved.fasta

#echo "Counting the number of characters in each line of the alignment..."
#cat /nfs3/Sharpton_Lab/prod/prod_restructure/projects/arnoldhk/2025_HVR_Guide_Phylogenetic_Integration/reference_data/silva.seed_v138_2_Eukaryota_excluded.align | grep -v ">" | awk '{print length}' | sort -n | uniq

#echo "Filtering the alignment file to remove columns of gaps..."
#/local/cluster/bin/mothur "#filter.seqs(fasta=/nfs3/Sharpton_Lab/prod/prod_restructure/projects/arnoldhk/2025_HVR_Guide_Phylogenetic_Integration/reference_data/silva.seed_v138_2_Eukaryota_excluded.align, vertical=T, processors=20)"

#echo "Counting the number of characters in the filtered alignment"
#cat /nfs3/Sharpton_Lab/prod/prod_restructure/projects/arnoldhk/2025_HVR_Guide_Phylogenetic_Integration/reference_data/silva.seed_v138_2_Eukaryota_excluded.filter.fasta | grep -v ">" | awk '{print length}' | sort -n | uniq

#echo "Getting the SILVA fasta for full length sequences..."
#cat /nfs3/Sharpton_Lab/prod/prod_restructure/projects/arnoldhk/2025_HVR_Guide_Phylogenetic_Integration/reference_data/silva.seed_v138_2_Eukaryota_excluded.filter.fasta | sed 's/\.//g' | sed 's/-//g' | sed 's/>/>REF./g' > /nfs3/Sharpton_Lab/prod/prod_restructure/projects/arnoldhk/2025_HVR_Guide_Phylogenetic_Integration/reference_data/silva.seed_v138_2_Eukaryota_excluded.filter.refseqs.fasta

#echo "Renaming files to something more reasonable"
#ln -s silva.seed_v138_2_Eukaryota_excluded.filter.fasta ref.align
#ln -s silva.seed_v138_2_Eukaryota_excluded.filter.refseqs.fasta full_length.fasta

echo "Finding the reference E. coli sequence"
cat ref.align | grep -n "AE005174.Esch1125"
head -n 808 ref.align | tail -n 2 > e_coli_ref_AE005174_Esch1125.fasta
ln -s e_coli_ref_AE005174_Esch1125.fasta e_coli_ref.fasta

