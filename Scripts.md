# Scripts
## Sequences downloading
```shell script
linfasta allMers_20_to_50.fa /dev/stdout | process_ncbi_names.sh /dev/stdin sanitized_all_mers.fna
```
## Filtering the 2.3 million covid genomes into a set with good names (suitable for phylogeny)
```shell script
linfasta covid_2m_secondary.fasta /dev/stdout | process_ncbi_names.sh /dev/stdin sanitized_covid.fna
```
## This cleans up the genomes by removing extraneous or erroneous sequences and tails
```shell script
awk '{if (NR % 2 == 0) {gsub("[BDEFHIJKLMNOPQRSUVWXYZ]+","N"); sub("^N",""); sub("N$",""); sub("AAAAA+$",""); sub("AAAAAAAAAAAAA.*",""); sub("TTTTTTTTTTTTTTT+",""); print $0;} else print $0;}' sanitized_covid.fna > captain_covid.fna
```
## Sort the genomes by length so that we can operate on a small subset first
```shell script
sortgenome captain_covid.fna sanitized_covid.srt LEN
head -200000 sanitized_covid.srt > captain100k.fa
```
## Cluster the 100,000 longest genomes to find reps
```shell script
numactl -i all canolax4b -db captain100k.fa -o cap-100k.fa -k 14 -t 120 -fitC .001
```
## Those 100k genomes turned into only 339 reps! Now align all sequences to that set
```shell script
canolax4b -db cap-100k.fa -q sanitized_covid.srt -o map_captain100k_all.tsv -k 14 -t 120 -local -min2
```
## Extract the non-mapping set (only 46,426 mapped at > .001)
```shell script
awk -F'\t' '$3 > .001' map_captain100k_all.tsv | cut -f1 | sed 's/^/>/' > toClust.post100k.txt
grep -F -f toClust.post100k.txt -A1 --no-group-separator captain_covid.fna > toClust.post100k.fna
```
## Cluster that remaining set and merge the cluster sets
### (Note, because we partitioned by sequence length, this partitioning is guaranteed lossless)
```shell script
canolax4b -db toClust.post100k.fna -o post-100k-001.fa -k 14 -t 120 -fitC .001
cat post-100k-001.fa cap-100k.fa > final_covid_clust.fna
```
## Get the final mapping statistics to the final cluster set (1534 genomes)
```shell script
canolax4b -db final_covid_clust.fna -q sanitized_covid.srt -o map_captain_all_final.tsv -k 14 -t 120 -local -min2

cut -f2 map_captain_all_final.tsv | sort | uniq -c | sort -hr | sed 's/ *//' | sed 's/ /\t/' > final_covid_cluster_memb.tsv
```

## Mapping all MERS against all covid to find highest homology matches and focus on those
```shell script
time makeblastdb -dbtype nucl -in final_covid_clust.fna

time blastn -db final_covid_clust.fna -query clustered001.fa -outfmt "6 std qcovs" -num_threads 62 -subject_besthit -max_target_seqs 1 -max_hsps 1 -task blastn -reward 2 -penalty -7 -gapopen 2 -gapextend 4 -no_greedy -mt_mode 0 > blastAll-2-7.b6
```
## A more comprehensive list of high-identity matches is given by
```shell script
time blastn -db final_covid_clust.fna -query clustered001.fa -outfmt "6 std qcovs" -num_threads 62 -subject_besthit -max_target_seqs 1000000 -max_hsps 1000000 -task blastn -reward 1 -penalty -5 -gapopen 3 -gapextend 3 -no_greedy -mt_mode 0 > blastAll-1-5-1Mseqs.b6
```
## Now we do an MSA of the 99.8% "all corona" reps, 618 (including 1 outlier): 
```shell script
/virus $ export MAFFT_BINARIES=$HOME/build/mafft-main/core/; $HOME/build/mafft-main/core/mafft --globalpair --maxiterate 10 --thread 64 --linelength -1 --adjustdirection FullCorona17_covid_mers-002.fa > fullCorona_imp620.msa
```
## We created a multiple sequence alignment of 618 diverse coronaviruses.
```shell script
export MAFFT_BINARIES=$HOME/build/mafft-main/core/; $HOME/build/mafft-main/core/mafft --thread 32 --linelength -1 --nomemsave --adjustdirection --maxiterate 10 --add toAdd.fa "fullCorona_imp620 (copy).msa" > added_fullCorona_to_skeleton_10iter.msa

time raxml-ng --threads 48 --msa added_fullCorona_to_skeleton_10iter.msa --prefix FC1 --model GTR+G --all --seed 123 --bs-metric fbp,tbe
```

