#/!/bin/bash
qig "conda run -n snakemake snakemake ../../../../../tmp/LE1.CS.bowtie2.bam \
	-s align.smk --use-conda --conda-frontend conda --cores 16" \
	-zslot 1 \
	-zmem 40G \
	-zo snake.log -ze snake.err \
	-zmailaddrs lucy.mahony@partners.basf.com \
	-zasync 



# conda run -n snakemake snakemake -j1  ../../../../../tmp/CS.star.index -s align.smk --use-conda --conda-frontend conda --cores 16

