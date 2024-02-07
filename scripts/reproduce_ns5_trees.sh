#!/bin/bash

##########################################################################################
# Mapping Glycoprotein Structure Reveals Defining Events in the Evolution of the Flaviviridae
# ----------------------------------------------------------------------------------------
# Authors: Jonathon C.O. Mifsud, Spyros Lytras, Michael Oliver, Kamilla Toon, Vincenzo A. Costa,
#          Edward C. Holmes, Joe Grove
# -----------------------------------------------------------------------------------------
# This script is dedicated to creating the NS5b flaviviridae tree variations
# It was run on a .pbs cluster but the script has been adapted here to run on a local machine
##########################################################################################

# Load necessary modules
#module load trimal/1.4.1
#module load mafft/7.402
#module load muscle/5.1
#module load clustal-omega/1.2.4
#module load iqtree/1.6.12

# Function to perform sequence alignment using MAFFT, MUSCLE, and Clustal Omega
function align {
  local sequences="$1"
  local filename="$(basename -- "$sequences")"
  local mafft_outfile="${filename}_untrimmed_MAFFT_$(date '+%Y%m%d').fasta"
  local muscle_outfile="${filename}_untrimmed_MUSCLE_$(date '+%Y%m%d').fasta"
  local clustalo_outfile="${filename}_untrimmed_ClustalOmega_$(date '+%Y%m%d').fasta"

  echo "Performing MAFFT alignment..."
  mafft --ep 0 --genafpair --maxiterate 1000 --thread 6 --threadit 6 "$sequences" > "$mafft_outfile" || { echo "MAFFT alignment failed"; exit 1; }

  echo "Performing MUSCLE alignment..."
  muscle -in "$sequences" -out "$muscle_outfile" || { echo "MUSCLE alignment failed"; exit 1; }
  
  echo "Performing Clustal Omega alignment..."
  clustalo -i "$sequences" --threads=4 --auto -o "$clustalo_outfile" || { echo "Clustal Omega alignment failed"; exit 1; }
}

# Function to trim alignments using Trimal for different gap thresholds and conservation levels
function trim {
  local filename="$1"
  
  # Iterate through aligners and perform trimming
  for aligner in MAFFT MUSCLE ClustalOmega; do
    local untrimmed="${filename}_untrimmed_${aligner}_$(date '+%Y%m%d').fasta"
    local trimmed_prefix="${filename}_trimmed_${aligner}"

    # Iterate through gap threshold and conservation values
    for gap_threshold in 0.9 0.8 0.7; do
      for cons in 2.5 5 7.5 10 12.5 15 17.5 20; do
        local trimal_outfile="${trimmed_prefix}_cons${cons}_gt${gap_threshold}_$(date '+%Y%m%d').fasta"
        echo "Trimming with Trimal: $trimal_outfile"
        trimal -in "$untrimmed" -out "$trimal_outfile" -gt "$gap_threshold" -cons "$cons" || { echo "Trimal trimming failed"; exit 1; }
      done
    done

    # Perform gappyout trimming
    local trimal_gappyout_outfile="${trimmed_prefix}_gappyout_$(date '+%Y%m%d').fasta"
    echo "Gappyout trimming with Trimal: $trimal_gappyout_outfile"
    trimal -in "$untrimmed" -out "$trimal_gappyout_outfile" -gappyout || { echo "Trimal gappyout trimming failed"; exit 1; }
  done
}

# Function to execute phylogenetic analysis on trimmed alignments using IQ-TREE
function phylogeny {
  local filename="$1"
  local wd=$(dirname "$1")  # Assuming $1 is a path to the sequences file
  local input_files=($(find "$wd" -name "${filename}_trimmed_*_$(date '+%Y%m%d').fasta"))

  # Models for phylogenetic analysis
  local models=("LG" "LG+F+R10" "FLAVI")

  for file in "${input_files[@]}"; do
    for model in "${models[@]}"; do
      echo "Running IQ-TREE for $file with model $model"
      iqtree -s "$file" -m "$model" -bb 1000 -alrt 1000 -nt AUTO || { echo "IQ-TREE analysis failed"; exit 1; }
    done
  done
}

# Main execution block
sequences="/sequences/sequences_for_alignments/ns5b/flaviviridae_ns5.fasta"

# Ensure sequences variable is provided
if [[ -z "$sequences" ]]; then
  echo "Sequences variable is not set. Exiting."
  exit 1
fi

# Execute functions with the correct arguments
basename=$(basename -- "$sequences")
align "$sequences"
trim "$basename"
phylogeny "$basename"
