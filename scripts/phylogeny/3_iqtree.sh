#!/bin/bash
# runtime: ~3 min per marker
start=$(date +%s)
DATE="2023-06-24"

# Declare an array of string with type-- note NO COMMAS between elements in bash array
declare -a StringArray=("18sgappy")
# "its2" rbcl" "18sgappy" "18sfull", 18sb" "18sa" "cat" 

# Iterate the string array using for loop
for MARKER in ${StringArray[@]}; do
  echo $MARKER
  FILEPATH=~/ownCloud/proj/single_cell/data/phylo/${DATE}_$MARKER/aligns.fasta
  iqtree -s $FILEPATH -bb 1000 -nt AUTO
done
end=$(date +%s)
elapsed=$(($end-$start))
echo "Elapsed Time: $((elapsed/60)) min"