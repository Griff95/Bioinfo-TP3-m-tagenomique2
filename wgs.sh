#!/usr/bin/bash 


# indexation bowtie2
bowtie2-build databases/all_genome.fasta databases/all_genome.fasta
bowtie2-build databases/resfinder.fna databases/resfinder.fna



# question 1: identifier et quantifier les bactéries présentes dans votre échantillon
bowtie2 -p 8 -x databases/all_genome.fasta -1 fastq/EchD_R1.fastq.gz -2 fastq/EchD_R2.fastq.gz -S res/bowtie_echD.sam



# question 2: quantifier l'abondance de chaque bactérie en analysant le fichier .sam
# 	a) convertir fichier sam en bam
# 		l'option -@ 8 permet de définir le nombre de thread
# 		l'option -b permet d'obtenir un output au format bam
# 		l'option -S Ignored compatibility with previous samtools versions
#		l'option -1 permet la fast BAM compression
samtools view -@ 8 -bS1 res/bowtie_echD.sam > res/bowtie_echD.bam

# 	b) trier le fichier bam
# 		l'option -l 1 permet de définir le niveau de compression (0: uncompressed, 1: fastest but minimal compression, 9: best compression but slowest to write)
samtools sort -@ 8 -l 1 res/bowtie_echD.bam > res/bowtie_echD_sorted.bam

#	c) indexer le fichier bam
samtools index res/bowtie_echD_sorted.bam

#	d) extraction du comptage
samtools idxstats res/bowtie_echD_sorted.bam > res/bowtie_echD_sorted_idxstats

#	f) association gi --> annotation
grep ">" databases/all_genome.fasta | cut -f 2 -d ">" > res/association.tsv


# question 3: assembler le génome des bactéries présentes
#	l'option -t défini le nombre de thread
#	l'option -m défini la mémoire utilisée
megahit -1 fastq/EchD_R1.fastq.gz -2 fastq/EchD_R2.fastq.gz --k-list 21 -t 4 -m 0.5 -o res/megahit_result


# question 4: prédire les gènes présents sur les contigs avec prodigal
prodigal -i res/megahit_result/final.contigs.fa -d res/prodigal_result.fna


# question 5: sélectionner les gènes "complets"
sed "s:>:*\n>:g" res/prodigal_result.fna | sed -n "/partial=00/,/*/p"|grep -v "*" > res/genes_full.fna


# question 6: annoter les gènes "complets" controle la banque resfinder avec blastn
makeblastdb -in databases/resfinder.fna -dbtype nucl

blastn -query res/genes_full.fna -db databases/resfinder.fna -out res/blastn_results.out -outfmt "6 qseqid sseqid evalue pident qcovhsp" -perc_identity 0.8 -qcov_hsp_perc 0.8 -evalue 0.003 -num_threads 8
