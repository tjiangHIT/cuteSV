# Generation of simulated datasets

---
### Preparation
	$ gzip -d sim_del.bed.gz
	$ gzip -d sim_ins.bed.gz
	$ gzip -d sim_inv.bed.gz
	$ gzip -d sim_dup.bed.gz
	$ gzip -d sim_tra.bed.gz
	$ mkdir workdir && cd workdir
	$ mv ../sim_*.bed ../LASeR.bed ./
	$ git clone https://github.com/tjiangHIT/VISOR.git
	$ cd VISOR && python3 setup.py install && cd ..

### Genome Modification
	$ VISOR HACk -g reference_genome.fasta -bed sim_del.bed -o donor_genome_del
	$ VISOR HACk -g reference_genome.fasta -bed sim_ins.bed -o donor_genome_ins
	$ VISOR HACk -g reference_genome.fasta -bed sim_inv.bed -o donor_genome_inv
	$ VISOR HACk -g reference_genome.fasta -bed sim_dup.bed -o donor_genome_dup
	$ VISOR HACk -g reference_genome.fasta -bed sim_tra.bed -o donor_genome_tra

### Simulated Alignments Generation 
	$ VISOR LASeR -g reference_genome.fasta -s donor_genome_del -bed LASeR.bed -o data_del_5x -c 5 --threads 16 --noaddtag
	$ VISOR LASeR -g reference_genome.fasta -s donor_genome_del -bed LASeR.bed -o data_del_10x -c 10 --threads 16 --noaddtag
	$ VISOR LASeR -g reference_genome.fasta -s donor_genome_del -bed LASeR.bed -o data_del_20x -c 20 --threads 16 --noaddtag
	$ VISOR LASeR -g reference_genome.fasta -s donor_genome_del -bed LASeR.bed -o data_del_30x -c 30 --threads 16 --noaddtag

	$ VISOR LASeR -g reference_genome.fasta -s donor_genome_ins -bed LASeR.bed -o data_ins_5x -c 5 --threads 16 --noaddtag
	$ VISOR LASeR -g reference_genome.fasta -s donor_genome_ins -bed LASeR.bed -o data_ins_10x -c 10 --threads 16 --noaddtag
	$ VISOR LASeR -g reference_genome.fasta -s donor_genome_ins -bed LASeR.bed -o data_ins_20x -c 20 --threads 16 --noaddtag
	$ VISOR LASeR -g reference_genome.fasta -s donor_genome_ins -bed LASeR.bed -o data_ins_30x -c 30 --threads 16 --noaddtag

	$ VISOR LASeR -g reference_genome.fasta -s donor_genome_inv -bed LASeR.bed -o data_inv_5x -c 5 --threads 16 --noaddtag
	$ VISOR LASeR -g reference_genome.fasta -s donor_genome_inv -bed LASeR.bed -o data_inv_10x -c 10 --threads 16 --noaddtag
	$ VISOR LASeR -g reference_genome.fasta -s donor_genome_inv -bed LASeR.bed -o data_inv_20x -c 20 --threads 16 --noaddtag
	$ VISOR LASeR -g reference_genome.fasta -s donor_genome_inv -bed LASeR.bed -o data_inv_30x -c 30 --threads 16 --noaddtag

	$ VISOR LASeR -g reference_genome.fasta -s donor_genome_dup -bed LASeR.bed -o data_dup_5x -c 5 --threads 16 --noaddtag
	$ VISOR LASeR -g reference_genome.fasta -s donor_genome_dup -bed LASeR.bed -o data_dup_10x -c 10 --threads 16 --noaddtag
	$ VISOR LASeR -g reference_genome.fasta -s donor_genome_dup -bed LASeR.bed -o data_dup_20x -c 20 --threads 16 --noaddtag
	$ VISOR LASeR -g reference_genome.fasta -s donor_genome_dup -bed LASeR.bed -o data_dup_30x -c 30 --threads 16 --noaddtag

	$ VISOR LASeR -g reference_genome.fasta -s donor_genome_tra -bed LASeR.bed -o data_tra_5x -c 5 --threads 16 --noaddtag
	$ VISOR LASeR -g reference_genome.fasta -s donor_genome_tra -bed LASeR.bed -o data_tra_10x -c 10 --threads 16 --noaddtag
	$ VISOR LASeR -g reference_genome.fasta -s donor_genome_tra -bed LASeR.bed -o data_tra_20x -c 20 --threads 16 --noaddtag
	$ VISOR LASeR -g reference_genome.fasta -s donor_genome_tra -bed LASeR.bed -o data_tra_30x -c 30 --threads 16 --noaddtag
	
