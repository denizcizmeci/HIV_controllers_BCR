# Extended somatic hypermutation and clonal evolution in B-cells of HIV controllers with broadly neutralizing antibodies
# FASTA to clone, version - 3 June 2020

samples=(
	donor149812
	donor201441
	donor211774
	donor280008
	donor315504
	donor330183
	donor386576
	donor444154
	donor534694
	donor622800
	donor628655
	donor756587
	donor785360
	)

# for sample in "${samples[@]}"
# do
#     VAR1="L001_R1_001.fastq.gz"
# 	VAR2="L001_R2_001.fastq.gz"
# 	VAR3="./FASTQ/"

# 	f1="$VAR3$sample$VAR1"
# 	f2="$VAR3$sample$VAR2"
# 	fa=".fasta"
# 	fastapath="./FASTA/"
# 	outfile="$fastapath$sample$fa"

# 	pandaseq -f $f1 -r $f2 > $outfile
# done


for sample in "${samples[@]}"
do
    VAR1=".fasta"
	fastapath="./FASTA/"
	f1="$fastapath$sample$VAR1"

	AssignGenes.py igblast -s $f1 -b ~/share/igblast \
    --organism human --loci ig --format blast

done


for sample in "${samples[@]}"
do
    VAR1="_igblast.fmt7"
    VAR2="./FASTA/"
	f1="$VAR2$sample$VAR1"

	VAR3=".fasta"
    VAR4="./FASTA/"
	f2="$VAR4$sample$VAR3"

	VAR5="_igblast_db-pass.tsv"
	f3="$VAR2$sample$VAR5"

	VAR6="_igblast_db-pass_clone-pass.tsv"
	f4="$VAR2$sample$VAR6"

	VAR7="_igblast_db-pass_clone-pass_germ-pass.tsv"
	f5="$VAR2$sample$VAR7"

	VAR8="./heavy/"
	f6="$VAR8$sample$VAR7"


	MakeDb.py igblast -i $f1 -s $f2 \
    -r ~/share/germlines/imgt/human/vdj/imgt_human_IGHV.fasta ~/share/germlines/imgt/human/vdj/imgt_human_IGHD.fasta ~/share/germlines/imgt/human/vdj/imgt_human_IGHJ.fasta 

	DefineClones.py -d $f3 --act set --model ham \
    --norm len --dist 0.16

    CreateGermlines.py -d $f4 -g dmask --cloned \
    -r ~/share/germlines/imgt/human/vdj/imgt_human_IGHV.fasta ~/share/germlines/imgt/human/vdj/imgt_human_IGHD.fasta ~/share/germlines/imgt/human/vdj/imgt_human_IGHJ.fasta

    mv $f5 $f6
done


# docker run -it --workdir /data -v $(pwd):/data:z kleinstein/immcantation:devel bash
# BuildTrees.py -d donor_igblast_db-pass_clone-pass_germ-pass.tsv --collapse --nproc 2 --ncdr3 --clean all --minseq 3 --igphyml --omega e,ce






