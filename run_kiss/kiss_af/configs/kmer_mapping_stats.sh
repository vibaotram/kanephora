## summarize kmer mapping results 

## count number of MQ0 kmers
stats=$(ls /scratch/most_kmer/kiss_af/output/8.MAPPING/output_file.*_vs_CC1.8_v2_pseudomolecule_cat_sorted.bam.stats)

mq0=0
for s in $stats
do
	mq0_s=$(head -n 20 $s | tail -n 1 | cut -f3)
	mq0=$(($mq0+$mq0_s))
done

## count number of unmapped kmers
um=0
for s in $stats
do
	um_s=$(head -n 16 $s | tail -n 1 | cut -f3)
	um=$(($um+$um_s))
done

## count number of mapped kmers on each chromosome
pos=$(ls /scratch/most_kmer/kiss_af/output/9.KMERPOSITION/output_file.*_vs_CC1.8_v2_pseudomolecule_cat_KMERPOSITION.txt)
for c in "Chr01" "Chr02" "Chr03" "Chr04" "Chr05" "Chr06" "Chr07" "Chr08" "Chr09" "Chr10" "Chr11" "Contig"
do
	echo $c >> count.log
	chr=0
	for p in $pos; do echo $p >> count.log; chr_s=$(grep $c $p | wc -l); chr=$(($chr+$chr_s)); done
	echo -e "$c\t$chr" >> count_kmer_chr.txt
done &

