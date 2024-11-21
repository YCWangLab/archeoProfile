# 1. break WGS data (Rice_320) into 53bp kmer by FastK
for file1 in *_1.fastq.gz
do
  file2=$(echo "$file1" | sed 's/_1/_2/g')
  bname=`echo $file1|cut -d_ -f1`
  FastK -T16 -NDB/$bname -P. -M32 -k53 -t1 $file1 $file2
done

# 2. get high copy kmer from DB, merge all high copy kmer to OUT/HighCopyMers
input_file=kmer_count.txt  # including 2 columns, the first is sample ID, the second is copy number
sh Logex_merge.sh $input_file
########Logex_merge.sh##########
#!/bin/bash

input_file=$1

line_number=0
lines=()

while IFS= read -r line || [[ -n "$line" ]]; do
    lines+=("$line")  
    ((line_number++)) 
    # process every eight lines
    if (( line_number % 8 == 0 )); then
        echo "Processing lines $((line_number - 7)) to $line_number:"
        params=()
        for current_line in "${lines[@]}"; do
            column1=$(echo "$current_line" | awk '{print $1}')
            column2=$(echo "$current_line" | awk '{print $2}')
            params+=("$column1" "$column2")
        done

        for param in "${params[@]}"
        do
          echo "$param"
          Logex -T16 -h 'OUT/Merge_'${params[0]}' = A['${params[1]}'-] |+ B['${params[3]}'-] |+ C['${params[5]}'-] |+ D['${params[7]}'-] |+ E['${params[9]}'-] |+ F['${params[11]}'-] |+ G['${params[13]}'-] |+ H['${params[15]}'-]' DB/${params[0]} DB/${params[2]} DB/${params[4]} DB/${params[6]} DB/${params[8]} DB/${params[10]} DB/${params[12]} DB/${params[14]}
          break
        done       
       lines=()
    fi
done < $input_file

# process the remaining lines
if (( ${#lines[@]} > 0 )); then
    params=()
    for current_line in "${lines[@]}"; do
        column1=$(echo "$current_line" | awk '{print $1}')
        column2=$(echo "$current_line" | awk '{print $2}')

        params+=("$column1" "$column2")
    done

    for param in "${params[@]}"; do
        echo "$param"
        Logex -T16 -h 'OUT/Merge_'${params[0]}' = A['${params[1]}'-] |+ B['${params[3]}'-] |+ C['${params[5]}'-] |+ D['${params[7]}'-] |+ E['${params[9]}'-] |+ F['${params[11]}'-] |+ G['${params[13]}'-]' DB/${params[0]} DB/${params[2]} DB/${params[4]} DB/${params[6]} DB/${params[8]} DB/${params[10]} DB/${params[12]}
	break
    done

fi
############################################
 
 Logex -T16 -h 'OUT/merge_1 = A |+ B |+ C |+ D |+ E |+ F |+ G |+ H' OUT/Merge_R0 OUT/Merge_S105 OUT/Merge_S112 OUT/Merge_S127 OUT/Merge_S12 OUT/Merge_S134 OUT/Merge_S141 OUT/Merge_S149
 Logex -T16 -h 'OUT/merge_2 = A |+ B |+ C |+ D |+ E |+ F |+ G |+ H' OUT/Merge_S156 OUT/Merge_S163 OUT/Merge_S170 OUT/Merge_S178 OUT/Merge_S186 OUT/Merge_S193 OUT/Merge_S207 OUT/Merge_S20
 Logex -T16 -h 'OUT/merge_3 = A |+ B |+ C |+ D |+ E |+ F |+ G |+ H' OUT/Merge_S214 OUT/Merge_S221 OUT/Merge_S229 OUT/Merge_S236 OUT/Merge_S243 OUT/Merge_S250 OUT/Merge_S258 OUT/Merge_S266
 Logex -T16 -h 'OUT/merge_4 = A |+ B |+ C |+ D |+ E |+ F |+ G |+ H' OUT/Merge_S273 OUT/Merge_S280 OUT/Merge_S288 OUT/Merge_S295 OUT/Merge_S301 OUT/Merge_S309 OUT/Merge_S316 OUT/Merge_S35
 Logex -T16 -h 'OUT/merge_5 = A |+ B |+ C |+ D |+ E |+ F |+ G |+ H' OUT/Merge_S42 OUT/Merge_S57 OUT/Merge_S5 OUT/Merge_S64 OUT/Merge_S71 OUT/Merge_S79 OUT/Merge_S86 OUT/Merge_S93
 Logex -T16 -h 'OUT/HighCopyMers = A |+ B |+ C |+ D |+ E' OUT/merge_1 OUT/merge_2 OUT/merge_3 OUT/merge_4 OUT/merge_5

# 3. extract HighCopyMers from every WGS data
for file in ./DB/*ktab
do
 bname=`basename $file|cut -d. -f1`
 # extract kmer from 'A' that also exist in 'B',  '&.' sign take kmer count in A for the output
 Logex -T16 -h 'OUT/'${bname}' = A &. B' DB/$bname OUT/HighCopyMers
done

 # convert kmer table (*ktab) to fasta
for file in ./OUT/S*ktab
do
  bname=`basename $file|cut -d. -f1`
  Tabex $file LIST > $bname.fa
  awk 'NR > 1 {print $1, $2, $4}' $bname.fa > tmp && mv tmp $bname.fa
done

for file in *fa
do
  awk '{print $1 $3,$2}' $file | tr " " "\n" | sed '1~2s/^/>/' > tmp && mv tmp $file
done

# 4. get population-specific kmer
 
 # Merge individuals belonging to the same population, *id file including  individuals' sample ID of specific population 
for file in *_id
do
  bname=`echo $file|cut -d_ -f1`
  cat $(< "$file") > ${bname}_out
  seqkit rmdup ${bname}_out -s -D ${bname}_dup > ${bname}.fa
done

 # select kmers that exist in more than 80%(85%,90%,95%) individuals
mkdir 80percent_copy
for file in *_out
do 
  bname=`echo $file|cut -d_ -f1`
  awk 'NR % 2 == 0' $file | sort | uniq -c > ./80percent_copy/${bname}_tmp
done

file=individual_percent  # the first line is population, the second line is 80% individual number of the population
tail -n +2 $file|while read line; do 
  pop=`echo $line|cut -d" " -f1`
  per1=`echo $line|cut -d" " -f2`
  cat ./80percent_copy/${pop}_tmp | awk '$1 >= '$per1'' | sed 's/^[ \t]*//; s/[ \t]\+/\n/g' > ./${pop}_80_seq
done

for file in ./80percent_copy/*seq
do
  bname=`basename $file|cut -d_ -f1-2`
  ID=`basename $file|cut -d_ -f1`
  awk 'NR % 2 == 1 {print ">'$ID'" $0 "_" NR; next} {print}' $file > ./80percent_copy/$bname.fa
done

python extract_specific.py aus_80%.fa rayada_80%.fa aromatic_80%.fa nivara_80%.fa indica_80%.fa temperate_80%.fa tropical_80%.fa rufipogon_80%.fa

###########extract_specific.py#####################
import sys
from Bio import SeqIO

if len(sys.argv) < 2:
    print("Usage: python extract_specific.py <file1.fa> <file2.fa> ... <fileN.fa>")
    sys.exit(1)

fasta_files = sys.argv[1:] 
all_sequences = {}

for index, fasta_file in enumerate(fasta_files):
    current_sequences = {str(record.seq): record.id for record in SeqIO.parse(fasta_file, "fasta")}
    all_sequences[index] = current_sequences
unique_sequences = {}

for i in range(len(fasta_files)):
    current_sequences = set(all_sequences[i].keys())
    for j in range(len(fasta_files)):
        if i != j:
            current_sequences -= set(all_sequences[j].keys())
    unique_sequences[i] = {seq: all_sequences[i][seq] for seq in current_sequences}

for i, unique in unique_sequences.items():
    if unique: 
        with open(f"specific_to_{fasta_files[i]}", "w") as output_file:
            for seq in unique:
                output_file.write(f">{unique[seq]}\n{seq}\n") 

print("unique sequences saved to the corresponding file")
#######################################################

# 5. mapping to pan-genome
bowtie2-build -f ref_66_genome.fa.gz --threads 16 ref_66_genome
DB=/home/usr/huangm/reference/66_rice_pangenome/ref_66_genome
for file in $fa_path/specific_to*fa
do
  bname=`basename $file | cut -d'_' -f3- | cut -d'.' -f1`
  sbatch -J $bname -o $bname.bowtie2.log ../bowtie2_fa_mapping.sh $bname $DB $file
done
##########bowtie2_fa_mapping.sh####################
#!/bin/bash
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=40
#SBATCH --mem=170000

module load bowtie2
module load samtools

bname=${1}
DB=${2}
file=${3}

echo "mapping $bname againest $DB"
bowtie2 --threads 35 -k 1000 -x $DB -f $file --no-unal | samtools view -@ 5 -b > $bname.bam
echo "done"



