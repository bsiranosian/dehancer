# my working direcrory for these files
wdir=/home/ben/ak_local/enhancer_conservation/encode_data/broad

# do this manually to fix and standardize some filenames
cd $wdir
cut -f 1,2,3 wgEncodeBroadHistoneDnd41H3k04me3Pk.broadPeak > Dnd41_H3k04me3.bed
cut -f 1,2,3 wgEncodeBroadHistoneDnd41H3k27acPk.broadPeak > Dnd41_H3k27ac.bed
cut -f 1,2,3 wgEncodeBroadHistoneGm12878H3k04me3StdPkV2.broadPeak > Gm12878_H3k04me3.bed
cut -f 1,2,3 wgEncodeBroadHistoneGm12878H3k27acStdPk.broadPeak > Gm12878_H3k27ac.bed
cut -f 1,2,3 wgEncodeBroadHistoneH1hescH3k27acStdPk.broadPeak > H1hesc_H3k27ac.bed
cut -f 1,2,3 wgEncodeBroadHistoneH1hescH3k4me3StdPk.broadPeak > H1hesc_H3k04me3.bed
cut -f 1,2,3 wgEncodeBroadHistoneHelas3H3k27acStdPk.broadPeak > Helas3_H3k27ac.bed
cut -f 1,2,3 wgEncodeBroadHistoneHelas3H3k4me3StdPk.broadPeak > Helas3_H3k04me3.bed
cut -f 1,2,3 wgEncodeBroadHistoneHepg2H3k27acStdPk.broadPeak > Hepg2_H3k27ac.bed
cut -f 1,2,3 wgEncodeBroadHistoneHepg2H3k4me3StdPk.broadPeak > Hepg2_H3k04me3.bed
cut -f 1,2,3 wgEncodeBroadHistoneHmecH3k27acStdPk.broadPeak > Hmec_H3k27ac.bed
cut -f 1,2,3 wgEncodeBroadHistoneHmecH3k4me3StdPk.broadPeak > Hmec_H3k04me3.bed
cut -f 1,2,3 wgEncodeBroadHistoneHsmmH3k27acStdPk.broadPeak > Hsmm_H3k27ac.bed
cut -f 1,2,3 wgEncodeBroadHistoneHsmmH3k4me3StdPk.broadPeak > Hsmm_H3k04me3.bed
cut -f 1,2,3 wgEncodeBroadHistoneHsmmtH3k27acStdPk.broadPeak > Hsmmt_H3k27ac.bed
cut -f 1,2,3 wgEncodeBroadHistoneHsmmtH3k4me3StdPk.broadPeak > Hsmmt_H3k04me3.bed
cut -f 1,2,3 wgEncodeBroadHistoneHuvecH3k27acStdPk.broadPeak > Huvec_H3k27ac.bed
cut -f 1,2,3 wgEncodeBroadHistoneHuvecH3k4me3StdPk.broadPeak > Huvec_H3k04me3.bed
cut -f 1,2,3 wgEncodeBroadHistoneK562H3k27acStdPk.broadPeak > K562_H3k27ac.bed
cut -f 1,2,3 wgEncodeBroadHistoneK562H3k4me3StdPk.broadPeak > K562_H3k04me3.bed
cut -f 1,2,3 wgEncodeBroadHistoneMonocd14ro1746H3k04me3Pk.broadPeak > Monocd14ro1746_H3k04me3.bed
cut -f 1,2,3 wgEncodeBroadHistoneMonocd14ro1746H3k27acPk.broadPeak > Monocd14ro1746_H3k27ac.bed
cut -f 1,2,3 wgEncodeBroadHistoneNhaH3k27acStdPk.broadPeak > Nha_H3k27ac.bed
cut -f 1,2,3 wgEncodeBroadHistoneNhaH3k4me3StdPk.broadPeak > Nha_H3k04me3.bed
cut -f 1,2,3 wgEncodeBroadHistoneNhdfadH3k27acStdPk.broadPeak > Nhdfad_H3k27ac.bed
cut -f 1,2,3 wgEncodeBroadHistoneNhdfadH3k4me3StdPk.broadPeak > Nhdfad_H3k04me3.bed
cut -f 1,2,3 wgEncodeBroadHistoneNhekH3k27acStdPk.broadPeak > Nhek_H3k27ac.bed
cut -f 1,2,3 wgEncodeBroadHistoneNhekH3k4me3StdPk.broadPeak > Nhek_H3k04me3.bed
cut -f 1,2,3 wgEncodeBroadHistoneNhlfH3k27acStdPk.broadPeak > Nhlf_H3k27ac.bed
cut -f 1,2,3 wgEncodeBroadHistoneNhlfH3k4me3StdPk.broadPeak > Nhlf_H3k04me3.bed
cut -f 1,2,3 wgEncodeBroadHistoneOsteoblH3k27acStdPk.broadPeak > Osteobl_H3k27ac.bed
cut -f 1,2,3 wgEncodeBroadHistoneOsteoH3k04me3Pk.broadPeak > Osteobl_H3k04me3.bed

mkdir original 
mv *.broadPeak original

# sort bedfiles
for i in $(ls *.bed); do bedtools sort -i $i > $i.sorted; mv -f $i.sorted $i; done

# make intersection set for enhancers and promoters
mkdir intersect
# enhancers have K3k27ac AND NOT H3k04me3
# promoters have K3k27ac AND H3k04me3
# the merge of these sets will be used to make a complement for the negative set 
# peaks closer than merge_dist will be merged int one
merge_dist=500
# and for the negative set they wll be extended by slop
slop=500

# location of the genome bedfile
genome_bed=/home/ben/ak_local/enhancer_conservation/hg19.genome
# cells=("Dnd41" "Gm12878" "H1hesc" "Helas3" "Hepg2" "Hmec" "Hsmm" "Hsmmt" "Huvec" "K562" "Monocd14ro1746" "Nha" "Nhdfad" "Nhek" "Nhlf" "Osteobl" )
cells=($(ls *H3k27ac* | cut -f1 -d'_'))
for cell in "${cells[@]}"; do
	# bedtools intersect -a $cell"_H3k27ac.bed" -b $cell"_H3k04me3.bed" -v > intersect/$cell"_enhancers.bed"
	# bedtools intersect -a $cell"_H3k27ac.bed" -b $cell"_H3k04me3.bed" > intersect/$cell"_promoters.bed"
	cat $cell"_H3k27ac.bed" $cell"_H3k04me3.bed" | bedtools sort | bedtools merge -d $merge_dist | bedtools slop -g $genome_bed -b $slop > intersect/$cell"_mergepeak.bed"
done

# might want some more intelligent way to make the negative set of intervals
# only on mapped, non repeat DNA or something? 
# for now we just do a completed 