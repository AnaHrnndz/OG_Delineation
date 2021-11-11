#!/bin/bash

TREE=$1
ALG=$2
OUT_DIR=$3
CPU=$4
taxonomy=$5
reftre=$6
NAME=$(basename $TREE)
echo $NAME


echo '#create groups and get fastas'
python /data/projects/find_dups/findDups_countLoss.py --tree $TREE --raw_fasta $ALG \
   --output_path $OUT_DIR --midpoint no  --taxonomy $taxonomy --reftree $reftre

cd $OUT_DIR

grep -c '^>' *.faa >initial_size.txt

echo '#Create aln for each group'
#for line in `find *_*.faa`; do mafft --auto $line > $line.aln; done
snakemake --snakefile /data/projects/find_dups/mafft.snakefile --cores $CPU

echo '#create hmm for each group'
#for line in `find *.aln`; do hmmbuild $line.hmm $line; done
snakemake --snakefile /data/projects/find_dups/hmmbuild.snakefile --cores $CPU

echo '#create hmm_groups database'
cat *.hmm >$NAME.hmm
hmmpress $NAME.hmm

hmmscan  --max --tblout $NAME.out.tsv $NAME.hmm *.notOG.faa
awk '!/^ *#/ {print}' $NAME.out.tsv > $NAME.clean.tsv

python3.6 /data/projects/find_dups/parse_hmmscan.py $NAME.clean.tsv result_hmmscan.tsv

echo '#expand hmm_groups with hmmscan results'
python3.6 /data/projects/find_dups/expand_og.py result_hmmscan.tsv  *.notOG.faa
grep -c '^>' *.faa >final_size.txt

# echo '#draw final tree to see reasingments'
# QT_QPA_PLATFORM=offscreen xvfb-run python3.6 ../distribution_hmm.py ./ $TREE $GNAMES P53 ./FINALtree.$NAME.pdf ./not_og-$NAME.faa 7742 midpoint


