#!/bin/bash
# AUTHOR: ARNOLD
# DATE 1/10/25
# PURPOSE: To simulate sequences from silva seed. 



echo "Simulating regions V1 - V9 at constant length..."
echo "Simulating V1..."
perl ../2025_hvr_guide_phylogenetic_integration/scripts/simulate_sequences_constant_length.pl -p V1.fasta -a ref.align -m e_coli_ref.fasta -l 50 -o seqs_V1_L50.fasta -e 50
perl ../2025_hvr_guide_phylogenetic_integration/scripts/simulate_sequences_constant_length.pl -p V1.fasta -a ref.align -m e_coli_ref.fasta -l 100 -o seqs_V1_L100.fasta -e 100
perl ../2025_hvr_guide_phylogenetic_integration/scripts/simulate_sequences_constant_length.pl -p V1.fasta -a ref.align -m e_coli_ref.fasta -l 200 -o seqs_V1_L200.fasta -e 200
perl ../2025_hvr_guide_phylogenetic_integration/scripts/simulate_sequences_constant_length.pl -p V1.fasta -a ref.align -m e_coli_ref.fasta -l 300 -o seqs_V1_L300.fasta -e 300
perl ../2025_hvr_guide_phylogenetic_integration/scripts/simulate_sequences_constant_length.pl -p V1.fasta -a ref.align -m e_coli_ref.fasta -l 400 -o seqs_V1_L400.fasta -e 400
perl ../2025_hvr_guide_phylogenetic_integration/scripts/simulate_sequences_constant_length.pl -p V1.fasta -a ref.align -m e_coli_ref.fasta -l 500 -o seqs_V1_L500.fasta -e 500

echo "Simulating V2..."
perl ../2025_hvr_guide_phylogenetic_integration/scripts/simulate_sequences_constant_length.pl -p V2.fasta -a ref.align -m e_coli_ref.fasta -l 50 -o seqs_V2_L50.fasta -e 50
perl ../2025_hvr_guide_phylogenetic_integration/scripts/simulate_sequences_constant_length.pl -p V2.fasta -a ref.align -m e_coli_ref.fasta -l 100 -o seqs_V2_L100.fasta -e 100
perl ../2025_hvr_guide_phylogenetic_integration/scripts/simulate_sequences_constant_length.pl -p V2.fasta -a ref.align -m e_coli_ref.fasta -l 200 -o seqs_V2_L200.fasta -e 200
perl ../2025_hvr_guide_phylogenetic_integration/scripts/simulate_sequences_constant_length.pl -p V2.fasta -a ref.align -m e_coli_ref.fasta -l 300 -o seqs_V2_L300.fasta -e 300
perl ../2025_hvr_guide_phylogenetic_integration/scripts/simulate_sequences_constant_length.pl -p V2.fasta -a ref.align -m e_coli_ref.fasta -l 400 -o seqs_V2_L400.fasta -e 400
perl ../2025_hvr_guide_phylogenetic_integration/scripts/simulate_sequences_constant_length.pl -p V2.fasta -a ref.align -m e_coli_ref.fasta -l 500 -o seqs_V2_L500.fasta -e 500

echo "Simulating V3..."
perl ../2025_hvr_guide_phylogenetic_integration/scripts/simulate_sequences_constant_length.pl -p V3.fasta -a ref.align -m e_coli_ref.fasta -l 50 -o seqs_V3_L50.fasta -e 50
perl ../2025_hvr_guide_phylogenetic_integration/scripts/simulate_sequences_constant_length.pl -p V3.fasta -a ref.align -m e_coli_ref.fasta -l 100 -o seqs_V3_L100.fasta -e 100
perl ../2025_hvr_guide_phylogenetic_integration/scripts/simulate_sequences_constant_length.pl -p V3.fasta -a ref.align -m e_coli_ref.fasta -l 200 -o seqs_V3_L200.fasta -e 200
perl ../2025_hvr_guide_phylogenetic_integration/scripts/simulate_sequences_constant_length.pl -p V3.fasta -a ref.align -m e_coli_ref.fasta -l 300 -o seqs_V3_L300.fasta -e 300
perl ../2025_hvr_guide_phylogenetic_integration/scripts/simulate_sequences_constant_length.pl -p V3.fasta -a ref.align -m e_coli_ref.fasta -l 400 -o seqs_V3_L400.fasta -e 400
perl ../2025_hvr_guide_phylogenetic_integration/scripts/simulate_sequences_constant_length.pl -p V3.fasta -a ref.align -m e_coli_ref.fasta -l 500 -o seqs_V3_L500.fasta -e 500

echo "Simulating V4..."
perl ../2025_hvr_guide_phylogenetic_integration/scripts/simulate_sequences_constant_length.pl -p V4.fasta -a ref.align -m e_coli_ref.fasta -l 50 -o seqs_V4_L50.fasta -e 50
perl ../2025_hvr_guide_phylogenetic_integration/scripts/simulate_sequences_constant_length.pl -p V4.fasta -a ref.align -m e_coli_ref.fasta -l 100 -o seqs_V4_L100.fasta -e 100
perl ../2025_hvr_guide_phylogenetic_integration/scripts/simulate_sequences_constant_length.pl -p V4.fasta -a ref.align -m e_coli_ref.fasta -l 200 -o seqs_V4_L200.fasta -e 200
perl ../2025_hvr_guide_phylogenetic_integration/scripts/simulate_sequences_constant_length.pl -p V4.fasta -a ref.align -m e_coli_ref.fasta -l 300 -o seqs_V4_L300.fasta -e 300
perl ../2025_hvr_guide_phylogenetic_integration/scripts/simulate_sequences_constant_length.pl -p V4.fasta -a ref.align -m e_coli_ref.fasta -l 400 -o seqs_V4_L400.fasta -e 400
perl ../2025_hvr_guide_phylogenetic_integration/scripts/simulate_sequences_constant_length.pl -p V4.fasta -a ref.align -m e_coli_ref.fasta -l 500 -o seqs_V4_L500.fasta -e 500

echo "Simulating V5..."
perl ../2025_hvr_guide_phylogenetic_integration/scripts/simulate_sequences_constant_length.pl -p V5.fasta -a ref.align -m e_coli_ref.fasta -l 50 -o seqs_V5_L50.fasta -e 50
perl ../2025_hvr_guide_phylogenetic_integration/scripts/simulate_sequences_constant_length.pl -p V5.fasta -a ref.align -m e_coli_ref.fasta -l 100 -o seqs_V5_L100.fasta -e 100
perl ../2025_hvr_guide_phylogenetic_integration/scripts/simulate_sequences_constant_length.pl -p V5.fasta -a ref.align -m e_coli_ref.fasta -l 200 -o seqs_V5_L200.fasta -e 200
perl ../2025_hvr_guide_phylogenetic_integration/scripts/simulate_sequences_constant_length.pl -p V5.fasta -a ref.align -m e_coli_ref.fasta -l 300 -o seqs_V5_L300.fasta -e 300
perl ../2025_hvr_guide_phylogenetic_integration/scripts/simulate_sequences_constant_length.pl -p V5.fasta -a ref.align -m e_coli_ref.fasta -l 400 -o seqs_V5_L400.fasta -e 400
perl ../2025_hvr_guide_phylogenetic_integration/scripts/simulate_sequences_constant_length.pl -p V5.fasta -a ref.align -m e_coli_ref.fasta -l 500 -o seqs_V5_L500.fasta -e 500

echo "Simulating V6..."
perl ../2025_hvr_guide_phylogenetic_integration/scripts/simulate_sequences_constant_length.pl -p V6.fasta -a ref.align -m e_coli_ref.fasta -l 50 -o seqs_V6_L50.fasta -e 50
perl ../2025_hvr_guide_phylogenetic_integration/scripts/simulate_sequences_constant_length.pl -p V6.fasta -a ref.align -m e_coli_ref.fasta -l 100 -o seqs_V6_L100.fasta -e 100
perl ../2025_hvr_guide_phylogenetic_integration/scripts/simulate_sequences_constant_length.pl -p V6.fasta -a ref.align -m e_coli_ref.fasta -l 200 -o seqs_V6_L200.fasta -e 200
perl ../2025_hvr_guide_phylogenetic_integration/scripts/simulate_sequences_constant_length.pl -p V6.fasta -a ref.align -m e_coli_ref.fasta -l 300 -o seqs_V6_L300.fasta -e 300
perl ../2025_hvr_guide_phylogenetic_integration/scripts/simulate_sequences_constant_length.pl -p V6.fasta -a ref.align -m e_coli_ref.fasta -l 400 -o seqs_V6_L400.fasta -e 400
perl ../2025_hvr_guide_phylogenetic_integration/scripts/simulate_sequences_constant_length.pl -p V6.fasta -a ref.align -m e_coli_ref.fasta -l 500 -o seqs_V6_L500.fasta -e 500

echo "Simulating V7..."
perl ../2025_hvr_guide_phylogenetic_integration/scripts/simulate_sequences_constant_length.pl -p V7.fasta -a ref.align -m e_coli_ref.fasta -l 50 -o seqs_V7_L50.fasta -e 50
perl ../2025_hvr_guide_phylogenetic_integration/scripts/simulate_sequences_constant_length.pl -p V7.fasta -a ref.align -m e_coli_ref.fasta -l 100 -o seqs_V7_L100.fasta -e 100
perl ../2025_hvr_guide_phylogenetic_integration/scripts/simulate_sequences_constant_length.pl -p V7.fasta -a ref.align -m e_coli_ref.fasta -l 200 -o seqs_V7_L200.fasta -e 200
perl ../2025_hvr_guide_phylogenetic_integration/scripts/simulate_sequences_constant_length.pl -p V7.fasta -a ref.align -m e_coli_ref.fasta -l 300 -o seqs_V7_L300.fasta -e 300
perl ../2025_hvr_guide_phylogenetic_integration/scripts/simulate_sequences_constant_length.pl -p V7.fasta -a ref.align -m e_coli_ref.fasta -l 400 -o seqs_V7_L400.fasta -e 400
perl ../2025_hvr_guide_phylogenetic_integration/scripts/simulate_sequences_constant_length.pl -p V7.fasta -a ref.align -m e_coli_ref.fasta -l 500 -o seqs_V7_L500.fasta -e 500

echo "Simulating V8..."
perl ../2025_hvr_guide_phylogenetic_integration/scripts/simulate_sequences_constant_length.pl -p V8.fasta -a ref.align -m e_coli_ref.fasta -l 50 -o seqs_V8_L50.fasta -e 50
perl ../2025_hvr_guide_phylogenetic_integration/scripts/simulate_sequences_constant_length.pl -p V8.fasta -a ref.align -m e_coli_ref.fasta -l 100 -o seqs_V8_L100.fasta -e 100
perl ../2025_hvr_guide_phylogenetic_integration/scripts/simulate_sequences_constant_length.pl -p V8.fasta -a ref.align -m e_coli_ref.fasta -l 200 -o seqs_V8_L200.fasta -e 200
perl ../2025_hvr_guide_phylogenetic_integration/scripts/simulate_sequences_constant_length.pl -p V8.fasta -a ref.align -m e_coli_ref.fasta -l 300 -o seqs_V8_L300.fasta -e 300
perl ../2025_hvr_guide_phylogenetic_integration/scripts/simulate_sequences_constant_length.pl -p V8.fasta -a ref.align -m e_coli_ref.fasta -l 400 -o seqs_V8_L400.fasta -e 400
perl ../2025_hvr_guide_phylogenetic_integration/scripts/simulate_sequences_constant_length.pl -p V8.fasta -a ref.align -m e_coli_ref.fasta -l 500 -o seqs_V8_L500.fasta -e 500

echo "Simulating V9..."
perl ../2025_hvr_guide_phylogenetic_integration/scripts/simulate_sequences_constant_length.pl -p V9.fasta -a ref.align -m e_coli_ref.fasta -l 50 -o seqs_V9_L50.fasta -e 50
perl ../2025_hvr_guide_phylogenetic_integration/scripts/simulate_sequences_constant_length.pl -p V9.fasta -a ref.align -m e_coli_ref.fasta -l 100 -o seqs_V9_L100.fasta -e 100
perl ../2025_hvr_guide_phylogenetic_integration/scripts/simulate_sequences_constant_length.pl -p V9.fasta -a ref.align -m e_coli_ref.fasta -l 200 -o seqs_V9_L200.fasta -e 200
perl ../2025_hvr_guide_phylogenetic_integration/scripts/simulate_sequences_constant_length.pl -p V9.fasta -a ref.align -m e_coli_ref.fasta -l 300 -o seqs_V9_L300.fasta -e 300
perl ../2025_hvr_guide_phylogenetic_integration/scripts/simulate_sequences_constant_length.pl -p V9.fasta -a ref.align -m e_coli_ref.fasta -l 400 -o seqs_V9_L400.fasta -e 400
perl ../2025_hvr_guide_phylogenetic_integration/scripts/simulate_sequences_constant_length.pl -p V9.fasta -a ref.align -m e_coli_ref.fasta -l 500 -o seqs_V9_L500.fasta -e 500


echo "Simulating regions V1 - V9 by hvr..."
perl  ../2025_hvr_guide_phylogenetic_integration/scripts/simulate_sequences_hvr.pl -a ref.align -o seqs_V1_LHVR.fasta -s 105 -t 632 
perl  ../2025_hvr_guide_phylogenetic_integration/scripts/simulate_sequences_hvr.pl -a ref.align -o seqs_V2_LHVR.fasta -s 663 -t 1815 
perl  ../2025_hvr_guide_phylogenetic_integration/scripts/simulate_sequences_hvr.pl -a ref.align -o seqs_V3_LHVR.fasta -s 1816 -t 2364 
perl  ../2025_hvr_guide_phylogenetic_integration/scripts/simulate_sequences_hvr.pl -a ref.align -o seqs_V4_LHVR.fasta -s 2365 -t 3131
perl  ../2025_hvr_guide_phylogenetic_integration/scripts/simulate_sequences_hvr.pl -a ref.align -o seqs_V5_LHVR.fasta -s 3132 -t 3595 
perl  ../2025_hvr_guide_phylogenetic_integration/scripts/simulate_sequences_hvr.pl -a ref.align -o seqs_V6_LHVR.fasta -s 3596 -t 4010 
perl  ../2025_hvr_guide_phylogenetic_integration/scripts/simulate_sequences_hvr.pl -a ref.align -o seqs_V7_LHVR.fasta -s 4011 -t  4969
perl  ../2025_hvr_guide_phylogenetic_integration/scripts/simulate_sequences_hvr.pl -a ref.align -o seqs_V8_LHVR.fasta -s 4970 -t 5372 
perl  ../2025_hvr_guide_phylogenetic_integration/scripts/simulate_sequences_hvr.pl -a ref.align -o seqs_V9_LHVR.fasta -s 5373 -t  5806
perl  ../2025_hvr_guide_phylogenetic_integration/scripts/simulate_sequences_hvr.pl -a ref.align -o seqs_VFull.fasta -s 1 -t  5806


echo "Simulating regions V1 - V9 by hvr...DONE"

echo "Simulating common experimental regions by hvr"
perl  ../2025_hvr_guide_phylogenetic_integration/scripts/simulate_sequences_hvr.pl -a ref.align -o seqs_V1V2_LHVR.fasta -s 105 -t 1815
perl  ../2025_hvr_guide_phylogenetic_integration/scripts/simulate_sequences_hvr.pl -a ref.align -o seqs_V1V3_LHVR.fasta -s 105 -t 2364
perl  ../2025_hvr_guide_phylogenetic_integration/scripts/simulate_sequences_hvr.pl -a ref.align -o seqs_V3V4_LHVR.fasta -s 1816 -t 3131
perl  ../2025_hvr_guide_phylogenetic_integration/scripts/simulate_sequences_hvr.pl -a ref.align -o seqs_V3V5_LHVR.fasta -s 1816 -t 3595
perl  ../2025_hvr_guide_phylogenetic_integration/scripts/simulate_sequences_hvr.pl -a ref.align -o seqs_V4V5_LHVR.fasta -s 2365 -t 3595
perl  ../2025_hvr_guide_phylogenetic_integration/scripts/simulate_sequences_hvr.pl -a ref.align -o seqs_V5V7_LHVR.fasta -s 3132 -t 4969
perl  ../2025_hvr_guide_phylogenetic_integration/scripts/simulate_sequences_hvr.pl -a ref.align -o seqs_V7V9_LHVR.fasta -s 4011 -t 5806

cat full_length.fasta seqs_V1V2_LHVR.fasta > seqs_V1V2_LHVR_guides.fasta 
cat full_length.fasta seqs_V1V3_LHVR.fasta > seqs_V1V3_LHVR_guides.fasta 
cat full_length.fasta seqs_V3V4_LHVR.fasta > seqs_V3V4_LHVR_guides.fasta 
cat full_length.fasta seqs_V3V5_LHVR.fasta > seqs_V3V5_LHVR_guides.fasta 
cat full_length.fasta seqs_V4V5_LHVR.fasta > seqs_V4V5_LHVR_guides.fasta 
cat full_length.fasta seqs_V5V7_LHVR.fasta > seqs_V5V7_LHVR_guides.fasta 
cat full_length.fasta seqs_V7V9_LHVR.fasta > seqs_V7V9_LHVR_guides.fasta 

cat full_length.fasta seqs_VFull.fasta > seqs_vfull_guides.fasta 

echo "Adding guides to expeirmental data"
cat experimental_demo_guides_few.fasta experimental_demo_short.fasta > experimental_demo_short_guides_few.fasta
cat experimental_demo_guides_many.fasta experimental_demo_short.fasta > experimental_demo_short_guides_many.fasta
cat experimental_demo_guides_blast_99.fasta experimental_demo_short.fasta > experimental_demo_short_guides_blast_99.fasta



