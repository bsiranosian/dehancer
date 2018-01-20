# my working direcrory for these files
# wdir=/home/ben/ak_local/enhancer_conservation/encode_data/broad
wdir=/users/bsiranos/analysis/enhancer_conservation/encode_data/broad

# do this manually to fix and standardize some filenames
cd $wdir
mkdir bed
cut -f 1,2,3 wgEncodeBroadHistoneDnd41H3k04me3Pk.broadPeak > bed/Dnd41_H3k04me3.bed
cut -f 1,2,3 wgEncodeBroadHistoneDnd41H3k27acPk.broadPeak > bed/Dnd41_H3k27ac.bed
cut -f 1,2,3 wgEncodeBroadHistoneGm12878H3k04me3StdPkV2.broadPeak > bed/Gm12878_H3k04me3.bed
cut -f 1,2,3 wgEncodeBroadHistoneGm12878H3k27acStdPk.broadPeak > bed/Gm12878_H3k27ac.bed
cut -f 1,2,3 wgEncodeBroadHistoneH1hescH3k27acStdPk.broadPeak > bed/H1hesc_H3k27ac.bed
cut -f 1,2,3 wgEncodeBroadHistoneH1hescH3k4me3StdPk.broadPeak > bed/H1hesc_H3k04me3.bed
cut -f 1,2,3 wgEncodeBroadHistoneHelas3H3k27acStdPk.broadPeak > bed/Helas3_H3k27ac.bed
cut -f 1,2,3 wgEncodeBroadHistoneHelas3H3k4me3StdPk.broadPeak > bed/Helas3_H3k04me3.bed
cut -f 1,2,3 wgEncodeBroadHistoneHepg2H3k27acStdPk.broadPeak > bed/Hepg2_H3k27ac.bed
cut -f 1,2,3 wgEncodeBroadHistoneHepg2H3k4me3StdPk.broadPeak > bed/Hepg2_H3k04me3.bed
cut -f 1,2,3 wgEncodeBroadHistoneHmecH3k27acStdPk.broadPeak > bed/Hmec_H3k27ac.bed
cut -f 1,2,3 wgEncodeBroadHistoneHmecH3k4me3StdPk.broadPeak > bed/Hmec_H3k04me3.bed
cut -f 1,2,3 wgEncodeBroadHistoneHsmmH3k27acStdPk.broadPeak > bed/Hsmm_H3k27ac.bed
cut -f 1,2,3 wgEncodeBroadHistoneHsmmH3k4me3StdPk.broadPeak > bed/Hsmm_H3k04me3.bed
cut -f 1,2,3 wgEncodeBroadHistoneHsmmtH3k27acStdPk.broadPeak > bed/Hsmmt_H3k27ac.bed
cut -f 1,2,3 wgEncodeBroadHistoneHsmmtH3k4me3StdPk.broadPeak > bed/Hsmmt_H3k04me3.bed
cut -f 1,2,3 wgEncodeBroadHistoneHuvecH3k27acStdPk.broadPeak > bed/Huvec_H3k27ac.bed
cut -f 1,2,3 wgEncodeBroadHistoneHuvecH3k4me3StdPk.broadPeak > bed/Huvec_H3k04me3.bed
cut -f 1,2,3 wgEncodeBroadHistoneK562H3k27acStdPk.broadPeak > bed/K562_H3k27ac.bed
cut -f 1,2,3 wgEncodeBroadHistoneK562H3k4me3StdPk.broadPeak > bed/K562_H3k04me3.bed
cut -f 1,2,3 wgEncodeBroadHistoneMonocd14ro1746H3k04me3Pk.broadPeak > bed/Monocd14ro1746_H3k04me3.bed
cut -f 1,2,3 wgEncodeBroadHistoneMonocd14ro1746H3k27acPk.broadPeak > bed/Monocd14ro1746_H3k27ac.bed
cut -f 1,2,3 wgEncodeBroadHistoneNhaH3k27acStdPk.broadPeak > bed/Nha_H3k27ac.bed
cut -f 1,2,3 wgEncodeBroadHistoneNhaH3k4me3StdPk.broadPeak > bed/Nha_H3k04me3.bed
cut -f 1,2,3 wgEncodeBroadHistoneNhdfadH3k27acStdPk.broadPeak > bed/Nhdfad_H3k27ac.bed
cut -f 1,2,3 wgEncodeBroadHistoneNhdfadH3k4me3StdPk.broadPeak > bed/Nhdfad_H3k04me3.bed
cut -f 1,2,3 wgEncodeBroadHistoneNhekH3k27acStdPk.broadPeak > bed/Nhek_H3k27ac.bed
cut -f 1,2,3 wgEncodeBroadHistoneNhekH3k4me3StdPk.broadPeak > bed/Nhek_H3k04me3.bed
cut -f 1,2,3 wgEncodeBroadHistoneNhlfH3k27acStdPk.broadPeak > bed/Nhlf_H3k27ac.bed
cut -f 1,2,3 wgEncodeBroadHistoneNhlfH3k4me3StdPk.broadPeak > bed/Nhlf_H3k04me3.bed
cut -f 1,2,3 wgEncodeBroadHistoneOsteoblH3k27acStdPk.broadPeak > bed/Osteobl_H3k27ac.bed
cut -f 1,2,3 wgEncodeBroadHistoneOsteoH3k04me3Pk.broadPeak > bed/Osteobl_H3k04me3.bed

mkdir original 
mv *.broadPeak original

# sort bedfiles
for i in $(ls *.bed); do bedtools sort -i $i > $i.sorted; mv -f $i.sorted $i; done

# make intersection set for enhancers and promoters
mkdir intersect
# enhancers have K3k27ac AND NOT H3k04me3
# promoters have K3k27ac AND H3k04me3
# the merge of these sets will be used to make a complement for the negative set 

# extend positive peaks by this amount
slop_positive=1000

# location of the genome bedfile
genome_bed=/users/bsiranos/analysis/enhancer_conservation/hg19_sort.genome
# cells=("Dnd41" "Gm12878" "H1hesc" "Helas3" "Hepg2" "Hmec" "Hsmm" "Hsmmt" "Huvec" "K562" "Monocd14ro1746" "Nha" "Nhdfad" "Nhek" "Nhlf" "Osteobl" )
cells=($(ls bed/*H3k27ac* | cut -f1 -d'_'| cut -f 2 -d'/'))
for cell in "${cells[@]}"; do
	# with this method we miss enhancers and promoters that are less than 2kb in length
	# in addition, we also want regions that only have 1000bp inside the interval 
	# this can be fixed by extending enhancer and promoter regions by 100bp on either side
	# then any interval will have at least one window in it
	# and the end 1000 bp of a region will have a chance to get a window that extends
	# past the end of the interval by 1000bp
	# run this through a merge to concat any linked intervals after
	echo "creating positive and negative set: " $cell
	bedtools intersect -a bed/$cell"_H3k27ac.bed" -b bed/$cell"_H3k04me3.bed" -v | bedtools slop -g $genome_bed -b $slop_positive | bedtools sort | bedtools merge > intersect/$cell"_enhancers.bed"
	bedtools intersect -a bed/$cell"_H3k27ac.bed" -b bed/$cell"_H3k04me3.bed" | bedtools slop -g $genome_bed -b $slop_positive | bedtools sort | bedtools merge > intersect/$cell"_promoters.bed"
	cat bed/$cell"_H3k27ac.bed" bed/$cell"_H3k04me3.bed" | bedtools sort | bedtools merge  > intersect/$cell"_mergepeak.bed"
	# might want some more intelligent way to make the negative set of intervals
	# only on mapped, non repeat DNA or something? 
	# for now we just do a complement to get the DNA intervals
	bedtools complement -i intersect/$cell"_mergepeak.bed" -g $genome_bed > intersect/$cell"_negative.bed"
done

# and now we can use bedtools to make sliding windows of these
# set window and step size
window_size=2000
step_size=500

mkdir windows
for cell in "${cells[@]}"; do
	# for enhancers, promoters and negative set
	# but need to check that each interval is large enough - filter with that awk command
	echo "making windows: " $cell
	bedtools makewindows -b intersect/$cell"_enhancers.bed" -w $window_size -s $step_size | awk '(($3-$2)==2000) { print }' > windows/$cell"_enhancers_windows.bed"
	bedtools makewindows -b intersect/$cell"_promoters.bed" -w $window_size -s $step_size | awk '(($3-$2)==2000) { print }'> windows/$cell"_promoters_windows.bed"
	bedtools makewindows -b intersect/$cell"_negative.bed" -w $window_size -s $step_size | awk '(($3-$2)==2000) { print }'> windows/$cell"_negative_windows.bed"
done

# extract each of these from the genome
hg19_genome=/mnt/data/annotations/by_organism/human/hg19.GRCh37/hg19.genome.fa

mkdir fasta
for cell in "${cells[@]}"; do
	# for enhancers, promoters and negative set
	echo "extracing fasta: " $cell
	bedtools getfasta -fi $hg19_genome -bed windows/$cell"_enhancers_windows.bed" > fasta/$cell"_enhancers_windows.fa"
	bedtools getfasta -fi $hg19_genome -bed windows/$cell"_promoters_windows.bed" > fasta/$cell"_promoters_windows.fa"
	bedtools getfasta -fi $hg19_genome -bed windows/$cell"_negative_windows.bed" > fasta/$cell"_negative_windows.fa"
done


