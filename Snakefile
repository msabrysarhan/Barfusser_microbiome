import glob, os



configfile: "general_config.yaml"

#report: "report.html"
#Retrieve all samples in data directory
INPUT_DIR=config["INPUT_DIR"]
samples = []
for file in glob.glob(INPUT_DIR+"/*fastq.gz"):
    new_sample = file.split("/")[-1].split(".")[0]
    if new_sample not in samples:
        samples.append(new_sample)

##Snakemake rules
#hallo
rule all:
    input:
        fastp_out = expand(["out/qc_reads/{sample}/{sample}.merged.qc.fastq.gz", "out/qc_reads/{sample}/{sample}.R{num}.unmerged.qc.fastq.gz", "out/qc_reads/{sample}/{sample}.{mode}.qc_report.html"], sample = samples, num = [1,2], mode = ["PE","merged"]),
        lib_complex_out = expand(["out/qc_reads/{sample}/{sample}.uniqeness.hist", "out/qc_reads/{sample}/{sample}.uniqeness.pdf"], sample = samples),
        rmdup_out = expand(["out/qc_reads/{sample}/{sample}.rmdup.merged.qc.fastq.gz", "out/qc_reads/{sample}/{sample}.R{num}.rmdup.unmerged.qc.fastq.gz"], sample = samples, num = [1,2]),
        spades_out = expand(["out/assembly/meta_spades/{sample}/{sample}.contigs"], sample = samples),
        megahit_out = expand(["out/assembly/megahit/{sample}/{sample}.contigs/contigs.fasta","out/assembly/megahit/{sample}/{sample}.contigs"], sample = samples),
        filter_assembly_out = expand(["out/assembly/{assembler}/{sample}/{sample}.binning/scaffolds.1000.{end}"], sample = samples, end = ["fasta", "sorted.bam", "sorted.depth", "sorted.depth.raw", "fasta.gz"], assembler = ["megahit", "meta_spades"]),
        MetaBAT_out = expand(["out/assembly/{assembler}/{sample}/{sample}.binning/metabat"], sample = samples, assembler = ["megahit", "meta_spades"]),
        MaxBin_out = expand(["out/assembly/{assembler}/{sample}/{sample}.binning/maxbin"], sample = samples, assembler = ["megahit", "meta_spades"]),
        concoct_out = expand(["out/assembly/{assembler}/{sample}/{sample}.binning/concoct.coverage_table.tsv", "out/assembly/{assembler}/{sample}/{sample}.binning/concoct"], sample = samples, assembler = ["megahit", "meta_spades"]),
        dastool_out = expand(["out/assembly/{assembler}/{sample}/{sample}.binning/{sample}_DASTool_bins", "out/assembly/{assembler}/{sample}/{sample}.binning/{binner}.contigs2bins.tsv"], sample = samples, assembler = ["megahit", "meta_spades"], binner = ["metabat", "maxbin", "concoct"]),
        checkM_out = expand(["out/assembly/{assembler}/{sample}/{sample}.checkm", "out/assembly/{assembler}/{sample}/{sample}.checkm/{sample}.checkm.summary"], sample = samples, assembler = ["megahit", "meta_spades"]),
        gtdb_out = expand(["out/assembly/{assembler}/{sample}/{sample}.GTDB"], sample = samples, assembler = ["megahit", "meta_spades"]),
        human_rCRS_out = expand(["out/human/{sample}/{sample}.mito.sorted.bam"], sample = samples),
        human_hg19_out = expand(["out/human/{sample}/{sample}.auto.sorted.bam"], sample = samples),
        human_qualimap_out = expand(["out/human/{sample}/{sample}.{mode}.qualimap"], sample = samples, mode = ["auto", "mito"]),
        human_mapDamage_out = expand(["out/human/{sample}/{sample}.{mode}.mapDamage",
        "out/human/{sample}/{sample}.auto.sorted.rescaled.bam", 
        "out/human/{sample}/{sample}.mito.sorted.rescaled.bam",
        "out/human/{sample}/{sample}.{mode}.mapDamage/Fragmisincorporation_plot.pdf"], sample = samples, mode = ["auto", "mito"]),
        sex_out = expand(["out/human/{sample}/{sample}.skglnd.sex", "out/human/{sample}/{sample}.mtnk.sex"], sample = samples),
        schmutzi_out = expand(["out/human/{sample}/{sample}.mito.sorted.MD.bam", "out/human/{sample}/{sample}.schmutzi", "out/human/{sample}/{sample}.mito.fasta"], sample = samples),
        y_haplo_out = expand(["out/human/{sample}/Y_haplogroup/{sample}.y.vcf", "out/human/{sample}/Y_haplogroup/{sample}.yhaplo"], sample = samples),
        #"out/human/{sample}/Y_haplogroup/{sample}.Y_chr.hg", 
        mito_haplogroup_out = expand(["out/human/{sample}/{sample}.mito.sorted.rescaled.vcf", "out/human/{sample}/{sample}.mito.hg"], sample = samples),
        diamond_out = expand(["out/taxonomic_classifications/diamond/{sample}/{sample}.nr.tab", "out/taxonomic_classifications/diamond/{sample}/{sample}.discarded.nr.tab"], sample =  samples),
        megan_out = expand(["out/taxonomic_classifications/diamond/{sample}/{sample}.nr.rma6", "out/taxonomic_classifications/diamond/{sample}/{sample}.discarded.nr.rma6"], sample = samples),
        metaphlan_out = expand(["out/taxonomic_classifications/metaphlan/{sample}/{sample}.def.metaphlan", "out/taxonomic_classifications/metaphlan/{sample}/{sample}.F.def.metaphlan", "out/taxonomic_classifications/metaphlan/{sample}/{sample}.R.def.metaphlan", "out/taxonomic_classifications/metaphlan/{sample}/merged_abundance_table_species.txt", "out/taxonomic_classifications/metaphlan/{sample}/{sample}.heatmap.pdf"], sample = samples),
        kraken_out = expand(["out/taxonomic_classifications/kraken/{sample}/{sample}.kraken.out",
        "out/taxonomic_classifications/kraken/{sample}/{sample}.krona.html",
        "out/taxonomic_classifications/kraken/{sample}/{sample}.kraken.report",
        "out/taxonomic_classifications/kraken/{sample}/{sample}.kraken.bracken",
        "out/taxonomic_classifications/kraken/{sample}/{sample}.discarded.kraken.out",
        "out/taxonomic_classifications/kraken/{sample}/{sample}.discarded.krona.html",
        "out/taxonomic_classifications/kraken/{sample}/{sample}.discarded.kraken.report",
        "out/taxonomic_classifications/kraken/{sample}/{sample}.discarded.kraken.bracken"], sample = samples),
        busco_out = expand(["out/assembly/{assembler}/{sample}/{sample}.busco"], sample = samples, assembler = ["megahit", "meta_spades"]),
        chrX_contam_out = expand(["out/human/{sample}/angsd.chrX.contam", "out/human/{sample}/angsd.chrX.R.contam"], sample = samples),
        stats_out = expand(["out/qc_reads/{sample}/{sample}.seq_stats.txt", "out/qc_reads/{sample}/{sample}.seq_stats.txt.pdf"], sample = samples),
        

#----------------------------------------------------#

threads_num = {}
for i in samples:
    threads_num[i]=max(round(os.path.getsize("data/"+ str(i) +".R1.fastq.gz")/10**9), 4)

print(threads_num)
#-------------General Quality Control------------------#
#QC_check and trimming using Fastp
rule fastp_trim:
    """Quality trimming and PE mergeing"""
    input:
        reads1 = "data/{sample}.R1.fastq.gz",
        reads2 = "data/{sample}.R2.fastq.gz"
    output:
        reads_qc1 = "out/qc_reads/{sample}/{sample}.R1.qc.fastq.gz",
        reads_qc2 = "out/qc_reads/{sample}/{sample}.R2.qc.fastq.gz",
        report_PE = "out/qc_reads/{sample}/{sample}.PE.qc_report.html",
        json_PE = "out/qc_reads/{sample}/{sample}.PE.qc_report.json"
    envmodules:
        "fastp/0.20.1", "gcc"
    threads: lambda wildcards: threads_num[wildcards.sample]
    resources:
        mem_mb=lambda wildcards, input, attempt: max(2*input.size_mb, 10000) ** 1 * attempt
    shell:
        """
        fastp -i {input.reads1} -I {input.reads2} \
        -o {output.reads_qc1} -O {output.reads_qc2} \
        -h {output.report_PE} -j {output.json_PE} \
        --thread {threads} --cut_mean_quality 30
        """

#Estimation of the library complexity using bbmap
rule lib_complexity:
    input:
        reads_1 = "out/qc_reads/{sample}/{sample}.R1.qc.fastq.gz",
        reads_2 = "out/qc_reads/{sample}/{sample}.R2.qc.fastq.gz"
    output:
        hist = "out/qc_reads/{sample}/{sample}.uniqeness.hist",
        hist_plot = report("out/qc_reads/{sample}/{sample}.uniqeness.pdf", category = "Quality Control")
    envmodules:
        "bbmap", "R"
    shell:
        """
        bbcountunique.sh in={input.reads_1} in2={input.reads_2} out={output.hist}
        Rscript src/scripts/LibraryComplexityPlot.R {output.hist} {output.hist_plot}
        """

#Merging pair-end reads using Fastp
rule fastp_merge:
    """Quality trimming and PE mergeing"""
    input:
        reads1 = "data/{sample}.R1.fastq.gz",
        reads2 = "data/{sample}.R2.fastq.gz"
    output:
        reads_merged = "out/qc_reads/{sample}/{sample}.merged.qc.fastq.gz",
        reads_unmerged1 = "out/qc_reads/{sample}/{sample}.R1.unmerged.qc.fastq.gz",
        reads_unmerged2 = "out/qc_reads/{sample}/{sample}.R2.unmerged.qc.fastq.gz",
        report_merged = "out/qc_reads/{sample}/{sample}.merged.qc_report.html",
        json_merged = "out/qc_reads/{sample}/{sample}.merged.qc_report.json"
    envmodules:
        "fastp/0.20.1", "gcc"
    threads: lambda wildcards: threads_num[wildcards.sample]
    resources:
        mem_mb=lambda wildcards, input, attempt: max(2*input.size_mb, 10000) ** 1 * attempt
    shell:
        """
        fastp -i {input.reads1} -I {input.reads2} \
        -o {output.reads_unmerged1} -O {output.reads_unmerged2} \
        --merge --merged_out {output.reads_merged} \
        --overlap_len_require 10 -h {output.report_merged} \
        -j {output.json_merged} --length_required 25 --thread {threads}
        """

#Removing PCR duplicates using SeqKit
rule rmdup:
    """Removing PCR duplicates from the PE merged data"""
    input:
        reads_merged = "out/qc_reads/{sample}/{sample}.merged.qc.fastq.gz",
        reads1_unmerged = "out/qc_reads/{sample}/{sample}.R1.unmerged.qc.fastq.gz",
        reads2_unmerged = "out/qc_reads/{sample}/{sample}.R2.unmerged.qc.fastq.gz"
    output:
        dedup = "out/qc_reads/{sample}/{sample}.rmdup.merged.qc.fastq.gz",
        reads1_unmerged = "out/qc_reads/{sample}/{sample}.R1.rmdup.unmerged.qc.fastq.gz",
        reads2_unmerged = "out/qc_reads/{sample}/{sample}.R2.rmdup.unmerged.qc.fastq.gz"
    threads: lambda wildcards: threads_num[wildcards.sample]
    envmodules:
        "seqkit/2.0.0", "gcc"
    resources:
        mem_mb=lambda wildcards, input, attempt: max(2*input.size_mb, 10000) ** 1 * attempt
    shell:
        """
        seqkit rmdup -s -j {threads} {input.reads_merged} | gzip > {output.dedup}
        seqkit rmdup -s -j {threads} {input.reads1_unmerged} | gzip > {output.reads1_unmerged}
        seqkit rmdup -s -j {threads} {input.reads2_unmerged} | gzip > {output.reads2_unmerged}
        """

#----------------------------------------------------#

#-------Metagenomic Assembly and MAGs Taxonomy-------#



#print(os.path.getsize("out/qc_reads/{sample}/{sample}.R1.qc.fastq.gz"))

#num_threads = os.path.getsize("out/qc_reads/{sample}/{sample}.R1.qc.fastq.gz")
#print(num_threads)

#def num_threads(wildcards):
#    return threads * 150


rule MetaSPAdes:
    input:
        reads1 = "out/qc_reads/{sample}/{sample}.R1.qc.fastq.gz",
        reads2 = "out/qc_reads/{sample}/{sample}.R2.qc.fastq.gz"
    output:
        assembly_dir = directory("out/assembly/meta_spades/{sample}/{sample}.contigs")
    threads: lambda wildcards: threads_num[wildcards.sample]*2
    resources: 
        mem_mb = lambda wildcards, input, attempt: max(2*input.size_mb, 10000) ** 1.1 * attempt,
        disk_mb = lambda wildcards, input, attempt: max(2*input.size_mb, 10000) ** 1.1 * attempt
    envmodules:
        "spades", "gcc"
    shell:
        """
        spades.py -1 {input.reads1} -2 {input.reads2} -m {resources.mem_mb} -t {threads} -o {output.assembly_dir} --meta
        rm -r out/assembly/meta_spades/{wildcards.sample}/{wildcards.sample}.contigs/tmp
        rm -r out/assembly/meta_spades/{wildcards.sample}/{wildcards.sample}.contigs/corrected
        rm -r out/assembly/meta_spades/{wildcards.sample}/{wildcards.sample}.contigs/K55
        rm -r out/assembly/meta_spades/{wildcards.sample}/{wildcards.sample}.contigs/K33
        rm -r out/assembly/meta_spades/{wildcards.sample}/{wildcards.sample}.contigs/K21
        rm -r out/assembly/meta_spades/{wildcards.sample}/{wildcards.sample}.contigs/misc
        rm out/assembly/meta_spades/{wildcards.sample}/{wildcards.sample}.contigs/params.txt
        """

rule megahit:
    input:
        reads1 = "out/qc_reads/{sample}/{sample}.R1.qc.fastq.gz",
        reads2 = "out/qc_reads/{sample}/{sample}.R2.qc.fastq.gz"
    output:
        contigs = "out/assembly/megahit/{sample}/{sample}.contigs/contigs.fasta",
        assembly_dir = directory("out/assembly/megahit/{sample}/{sample}.contigs")
    threads: lambda wildcards: threads_num[wildcards.sample]*2
    resources:
        mem_mb=lambda wildcards, input, attempt: max(2*input.size_mb, 10000) ** 1.1 * attempt,
        disk_mb=lambda wildcards, input, attempt: max(2*input.size_mb, 10000) ** 1.1 * attempt
    envmodules:
        "megahit", "gcc"
    shell:
        """
        rm -rf {output.assembly_dir}
        megahit -1 {input.reads1} -2 {input.reads2} -t {threads} -o {output.assembly_dir}
        rm -r out/assembly/megahit/{wildcards.sample}/{wildcards.sample}.contigs/intermediate_contigs
        mv out/assembly/megahit/{wildcards.sample}/{wildcards.sample}.contigs/final.contigs.fa {output.contigs}
        """

rule filter_assembly:
    input:
        reads1 = "out/qc_reads/{sample}/{sample}.R1.qc.fastq.gz",
        reads2 = "out/qc_reads/{sample}/{sample}.R2.qc.fastq.gz",
        contigs = "out/assembly/{assembler}/{sample}/{sample}.contigs"
    output:
        contigs = "out/assembly/{assembler}/{sample}/{sample}.binning/scaffolds.1000.fasta",
        bam = "out/assembly/{assembler}/{sample}/{sample}.binning/scaffolds.1000.sorted.bam",
        depth = "out/assembly/{assembler}/{sample}/{sample}.binning/scaffolds.1000.sorted.depth",
        depth_raw = "out/assembly/{assembler}/{sample}/{sample}.binning/scaffolds.1000.sorted.depth.raw",
        zip_contigs = "out/assembly/{assembler}/{sample}/{sample}.binning/scaffolds.1000.fasta.gz"
    threads: lambda wildcards: threads_num[wildcards.sample]
    envmodules:
        "seqtk", "bowtie2", "samtools", "metabat", "gcc"
    resources:
        mem_mb=lambda wildcards, input, attempt: max(2*input.size_mb, 10000) ** 1 * attempt,
        disk_mb=lambda wildcards, input, attempt: max(2*input.size_mb, 10000) ** 1 * attempt
    shell:
        """
        seqtk seq -L 1000 {input.contigs}/contigs.fasta > {output.contigs}
        seqtk seq -L 1000 {input.contigs}/contigs.fasta | gzip > {output.zip_contigs}
        bowtie2-build {output.contigs} out/assembly/{wildcards.assembler}/{wildcards.sample}/{wildcards.sample}.binning/scaffolds.1000
        bowtie2 -x out/assembly/{wildcards.assembler}/{wildcards.sample}/{wildcards.sample}.binning/scaffolds.1000 -1 {input.reads1} -2 {input.reads2} -S out/assembly/{wildcards.assembler}/{wildcards.sample}/{wildcards.sample}.binning/scaffolds.1000.sam --no-unal --threads {threads}
        samtools view -@ {threads} -Sbq 30 out/assembly/{wildcards.assembler}/{wildcards.sample}/{wildcards.sample}.binning/scaffolds.1000.sam > out/assembly/{wildcards.assembler}/{wildcards.sample}/{wildcards.sample}.binning/scaffolds.1000.bam
        samtools sort -@ {threads} out/assembly/{wildcards.assembler}/{wildcards.sample}/{wildcards.sample}.binning/scaffolds.1000.bam > {output.bam}
        samtools index -@ {threads} {output.bam}
        jgi_summarize_bam_contig_depths {output.bam} --outputDepth {output.depth_raw}
        cat {output.depth_raw} | cut -f1,3 > {output.depth}
        rm out/assembly/{wildcards.assembler}/{wildcards.sample}/{wildcards.sample}.binning/scaffolds.1000.sam
        rm out/assembly/{wildcards.assembler}/{wildcards.sample}/{wildcards.sample}.binning/scaffolds.1000.bam
        rm out/assembly/{wildcards.assembler}/{wildcards.sample}/{wildcards.sample}.binning/scaffolds.1000.rev.2.bt2
        rm out/assembly/{wildcards.assembler}/{wildcards.sample}/{wildcards.sample}.binning/scaffolds.1000.rev.1.bt2
        rm out/assembly/{wildcards.assembler}/{wildcards.sample}/{wildcards.sample}.binning/scaffolds.1000.2.bt2
        rm out/assembly/{wildcards.assembler}/{wildcards.sample}/{wildcards.sample}.binning/scaffolds.1000.1.bt2
        rm out/assembly/{wildcards.assembler}/{wildcards.sample}/{wildcards.sample}.binning/scaffolds.1000.4.bt2
        rm out/assembly/{wildcards.assembler}/{wildcards.sample}/{wildcards.sample}.binning/scaffolds.1000.3.bt2
        """


rule MetaBAT:
    input:
        contigs = "out/assembly/{assembler}/{sample}/{sample}.binning/scaffolds.1000.fasta.gz",
        depth = "out/assembly/{assembler}/{sample}/{sample}.binning/scaffolds.1000.sorted.depth.raw"
    output:
        base = directory("out/assembly/{assembler}/{sample}/{sample}.binning/metabat")
    envmodules:
        "metabat", "gcc"
    resources:
        mem_mb=lambda wildcards, input, attempt: max(2*input.size_mb, 10000) ** 1 * attempt
    threads:
        4
    shell:
        """
        metabat -i {input.contigs} -a {input.depth} -o {output.base}/metabat -m 1500 --minClsSize 10000
        """


rule MaxBin:
    input:
        contigs = "out/assembly/{assembler}/{sample}/{sample}.binning/scaffolds.1000.fasta",
        depth = "out/assembly/{assembler}/{sample}/{sample}.binning/scaffolds.1000.sorted.depth"
    output:
        base = directory("out/assembly/{assembler}/{sample}/{sample}.binning/maxbin")
    envmodules:
        "maxbin","perl/5.36.0", "gcc", "bowtie2"
    threads:
        4
    resources:
        mem_mb=lambda wildcards, input, attempt: max(2*input.size_mb, 10000) ** 1 * attempt
    shell:
        """
        run_MaxBin.pl -contig {input.contigs} -abund {input.depth} -out {output.base}
        mkdir out/assembly/{wildcards.assembler}/{wildcards.sample}/{wildcards.sample}.binning/maxbin
        cd out/assembly/{wildcards.assembler}/{wildcards.sample}/{wildcards.sample}.binning
        for i in $(ls maxbin*.fasta | rev | cut -c 7- | rev); do mv $i.fasta maxbin/$i.fa; done
        """

rule CONCOCT:
    input:
        bam = "out/assembly/{assembler}/{sample}/{sample}.binning/scaffolds.1000.sorted.bam",
        contigs = "out/assembly/{assembler}/{sample}/{sample}.binning/scaffolds.1000.fasta",
        depth = "out/assembly/{assembler}/{sample}/{sample}.binning/scaffolds.1000.sorted.depth"
    output:
        cov_tab = "out/assembly/{assembler}/{sample}/{sample}.binning/concoct.coverage_table.tsv",
        concoct_dir = directory("out/assembly/{assembler}/{sample}/{sample}.binning/concoct")
    threads:
        4
    resources:
        mem_mb=lambda wildcards, input, attempt: max(2*input.size_mb, 10000) ** 1 * attempt
    envmodules:
        "samtools", "python3/3.7.0", "concoct/1.1.0-py37", "gcc", "compat"
    shell:
        """
        samtools index {input.bam}
        cut_up_fasta.py {input.contigs} -c 10000 -o 0 --merge_last -b out/assembly/{wildcards.assembler}/{wildcards.sample}/{wildcards.sample}.binning/contigs_10K.bed > out/assembly/{wildcards.assembler}/{wildcards.sample}/{wildcards.sample}.binning/contigs_10K.fna
        concoct_coverage_table.py out/assembly/{wildcards.assembler}/{wildcards.sample}/{wildcards.sample}.binning/contigs_10K.bed {input.bam} > {output.cov_tab}
        concoct --composition_file out/assembly/{wildcards.assembler}/{wildcards.sample}/{wildcards.sample}.binning/contigs_10K.fna --coverage_file {output.cov_tab} --threads 8 -b out/assembly/{wildcards.assembler}/{wildcards.sample}/{wildcards.sample}.binning/
        merge_cutup_clustering.py out/assembly/{wildcards.assembler}/{wildcards.sample}/{wildcards.sample}.binning/clustering_gt1000.csv > out/assembly/{wildcards.assembler}/{wildcards.sample}/{wildcards.sample}.binning/clustering_merged.csv
        mkdir {output.concoct_dir}
        extract_fasta_bins.py {input.contigs} out/assembly/{wildcards.assembler}/{wildcards.sample}/{wildcards.sample}.binning/clustering_merged.csv --output_path {output.concoct_dir}
        cd {output.concoct_dir}
        for file in $(ls *.fa)
        do
        minimumsize=50000
        actualsize=$(wc -c < $file)
        if [ $actualsize -le $minimumsize ]; then
            rm $file
        fi
        done
        for i in $(ls * | rev | cut -c 4- | rev)
        do
        mv $i.fa concoct.$i.fa
        done
        """


rule DAS_Tool:
    input:
        metabat = "out/assembly/{assembler}/{sample}/{sample}.binning/metabat",
        maxbin = "out/assembly/{assembler}/{sample}/{sample}.binning/maxbin",
        concoct = "out/assembly/{assembler}/{sample}/{sample}.binning/concoct",
        contigs = "out/assembly/{assembler}/{sample}/{sample}.binning/scaffolds.1000.fasta"
    output:
        metabat = "out/assembly/{assembler}/{sample}/{sample}.binning/metabat.contigs2bins.tsv",
        maxbin = "out/assembly/{assembler}/{sample}/{sample}.binning/maxbin.contigs2bins.tsv",
        concoct = "out/assembly/{assembler}/{sample}/{sample}.binning/concoct.contigs2bins.tsv",
        DAS_Tool = directory("out/assembly/{assembler}/{sample}/{sample}.binning/{sample}_DASTool_bins")
    envmodules:
        "dastool", "usearch", "R/3.5.1", "gcc"
    threads:
        4
    resources:
        mem_mb=lambda wildcards, input, attempt: max(2*input.size_mb, 10000) * 1 * attempt,
        disk_mb=lambda wildcards, input, attempt: max(2*input.size_mb, 10000) * 1 * attempt
    shell:
        """
        mkdir {output.DAS_Tool}
        src/scripts/Fasta_to_Scaffolds2Bin.sh -i {input.metabat} -e fa > {output.metabat}
        src/scripts/Fasta_to_Scaffolds2Bin.sh -i {input.maxbin} -e fa > {output.maxbin}
        src/scripts/Fasta_to_Scaffolds2Bin.sh -i {input.concoct} -e fa > {output.concoct}
        DAS_Tool -i {output.metabat},{output.maxbin},{output.concoct} -l metabat,maxbin,concoct  -c {input.contigs} -o out/assembly/{wildcards.assembler}/{wildcards.sample}/{wildcards.sample}.binning/{wildcards.sample} -t {threads} --write_bins 1 --score_threshold 0.0
        """


rule checkM:
    input:
        bins = "out/assembly/{assembler}/{sample}/{sample}.binning/{sample}_DASTool_bins"
    output:
        checkm = directory("out/assembly/{assembler}/{sample}/{sample}.checkm"),
        summary = "out/assembly/{assembler}/{sample}/{sample}.checkm/{sample}.checkm.summary"
    threads:
        8
    resources:
        mem_mb=lambda wildcards, input, attempt: max(2*input.size_mb, 10000) ** 1 * attempt,
        disk_mb=lambda wildcards, input, attempt: max(2*input.size_mb, 10000) ** 1 * attempt
    envmodules:
        "checkm", "gcc", "python3"
    shell:
        """
        checkm lineage_wf {input.bins} {output.checkm} -x ".fa" -t {threads} --pplacer_threads {threads} --tab_table > {output.summary}
        
        rm -rf {output.checkm}/storage
        rm -rf {output.checkm}/bins
        """


rule BUSCO:
    input:
        bins = "out/assembly/{assembler}/{sample}/{sample}.binning/{sample}_DASTool_bins"
    output:
        busco = directory("out/assembly/{assembler}/{sample}/{sample}.busco"),
        #summary = "out/assembly/{assembler}/{sample}/{sample}.busco/{sample}.checkm.summary"
    threads:
        8
    resources:
        mem_mb=lambda wildcards, input, attempt: max(2*input.size_mb, 10000) ** 1 * attempt,
        disk_mb=lambda wildcards, input, attempt: max(2*input.size_mb, 10000) ** 1 * attempt
    shell:
        """
        export PATH="/apps/augustus/3.4.0/bin:$PATH"
        export PATH="/apps/augustus/3.4.0/scripts:$PATH"
        export AUGUSTUS_CONFIG_PATH="/apps/augustus/3.4.0/config/"
        export BUSCO_CONFIG_FILE="busco_config.ini"
        ./src/scripts/BUSCO.sh {input.bins} {output.busco} {wildcards.sample}_busco
        #busco -i {input.bins} -f -o {wildcards.sample}_busco -m genome --auto-lineage --out_path {output.busco}
        rm -rf out/assembly/{wildcards.assembler}/{wildcards.sample}/{wildcards.sample}.busco/{wildcards.sample}_busco/*fa
        rm -rf out/assembly/{wildcards.assembler}/{wildcards.sample}/{wildcards.sample}.busco/{wildcards.sample}_busco/logs
        rm busco*log
        """

rule GTDB:
    input:
        bins = "out/assembly/{assembler}/{sample}/{sample}.binning/{sample}_DASTool_bins"
    output:
        GTDB = directory("out/assembly/{assembler}/{sample}/{sample}.GTDB")
    threads:
        8
    envmodules:
        "gtdbtk", "gcc"
    resources:
        mem_mb=lambda wildcards, input, attempt: max(2*input.size_mb, 10000) ** 1 * attempt,
        disk_mb=lambda wildcards, input, attempt: max(2*input.size_mb, 10000) ** 1 * attempt
    shell:
        """
        gtdbtk classify_wf --genome_dir {input.bins} --cpus {threads} -x fa --out_dir {output.GTDB} --force --min_perc_aa 5 --pplacer_cpus {threads}
        #rm -r {output.GTDB}/*/intermediate_results
        """

#-------------------------------------------------#

#-------------Human DNA Analysis------------------#
rule human_rCRS:
    input:
        reads = "out/qc_reads/{sample}/{sample}.rmdup.merged.qc.fastq.gz",
        ref = config["rCRS"]
    output:
        mito = "out/human/{sample}/{sample}.mito.sorted.bam"
    envmodules:
        "bwa", "samtools", "gcc"
    threads: lambda wildcards: threads_num[wildcards.sample]
    resources:
        mem_mb=lambda wildcards, input, attempt: max(2*input.size_mb, 10000) ** 1 * attempt,
        disk_mb=lambda wildcards, input, attempt: max(2*input.size_mb, 10000) ** 1 * attempt
    shell:
        """
        bwa aln -t {threads} {input.ref} {input.reads} > out/human/{wildcards.sample}/{wildcards.sample}.mito.sai
        bwa samse {input.ref} out/human/{wildcards.sample}/{wildcards.sample}.mito.sai {input.reads} | samtools view -q 30 -bSh > out/human/{wildcards.sample}/{wildcards.sample}.mito.bam
        samtools sort -@ {threads} out/human/{wildcards.sample}/{wildcards.sample}.mito.bam -o {output.mito} 
        samtools index {output.mito}
        rm out/human/{wildcards.sample}/{wildcards.sample}.mito.sai out/human/{wildcards.sample}/{wildcards.sample}.mito.bam
        """


rule human_hg19:
    input:
        ref = config["hg19"],
        reads = "out/qc_reads/{sample}/{sample}.rmdup.merged.qc.fastq.gz"
    output:
        human_reads = "out/human/{sample}/{sample}.auto.sorted.bam"
    threads: lambda wildcards: threads_num[wildcards.sample]
    envmodules:
        "bwa", "samtools", "gcc"
    resources:
        mem_mb=lambda wildcards, input, attempt: max(2*input.size_mb, 10000) ** 1 * attempt,
        disk_mb=lambda wildcards, input, attempt: max(2*input.size_mb, 10000) ** 1 * attempt
    shell:
        """
        bwa aln -t {threads} {input.ref} {input.reads} > out/human/{wildcards.sample}/{wildcards.sample}.auto.sai
        bwa samse {input.ref} out/human/{wildcards.sample}/{wildcards.sample}.auto.sai {input.reads} | samtools view -q 30 -bSh > out/human/{wildcards.sample}/{wildcards.sample}.auto.bam
        samtools sort out/human/{wildcards.sample}/{wildcards.sample}.auto.bam -o {output.human_reads}
        samtools index {output.human_reads}
        rm out/human/{wildcards.sample}/{wildcards.sample}.auto.sai out/human/{wildcards.sample}/{wildcards.sample}.auto.bam
        """

rule human_QualiMap:
    input:
        human_auto = "out/human/{sample}/{sample}.auto.sorted.bam",
        human_mito = "out/human/{sample}/{sample}.mito.sorted.bam"
    output:
        mito_qualimap = directory("out/human/{sample}/{sample}.mito.qualimap"),
        auto_qualimap = directory("out/human/{sample}/{sample}.auto.qualimap")
    envmodules:
        "qualimap", "gcc", "java"
    threads:
        4
    shell:
        """
		unset DISPLAY
        qualimap bamqc -bam {input.human_mito} -outdir {output.mito_qualimap}
        qualimap bamqc -bam {input.human_auto} -outdir {output.auto_qualimap}
        """

rule human_mapDamage:
    input:
        human_auto = "out/human/{sample}/{sample}.auto.sorted.bam",
        human_mito = "out/human/{sample}/{sample}.mito.sorted.bam",
        rCRS = config["rCRS"],
        hg19 = config["hg19"]
    output:
        mito_damage = directory("out/human/{sample}/{sample}.mito.mapDamage"),
        mito_plot = report("out/human/{sample}/{sample}.mito.mapDamage/Fragmisincorporation_plot.pdf", category= "MapDamage"),
        rescaled_mito = "out/human/{sample}/{sample}.mito.sorted.rescaled.bam",
        rescaled_auto = "out/human/{sample}/{sample}.auto.sorted.rescaled.bam",
        auto_damage = directory("out/human/{sample}/{sample}.auto.mapDamage"),
        auto_plot = report("out/human/{sample}/{sample}.auto.mapDamage/Fragmisincorporation_plot.pdf", category= "MapDamage"),
    envmodules:
        "mapdamage", "samtools", "python3", "R/3.5.1", "gcc"
    threads:
        8
    resources:
        mem_mb=lambda wildcards, input, attempt: max(2*input.size_mb, 10000) ** 1 * attempt,
        disk_mb = lambda wildcards, input, attempt: max(2*input.size_mb, 10000) ** 1 * attempt
    shell:
        """
        module load mapdamage
        module load R/3.5.1
        mapDamage -i {input.human_mito} -r {input.rCRS} -d {output.mito_damage} --rescale --rescale-out {output.rescaled_mito}
        samtools index {output.rescaled_mito}
        mapDamage -i {input.human_auto} -r {input.hg19} -d {output.auto_damage} --rescale --rescale-out {output.rescaled_auto}
        samtools index {output.rescaled_auto}
        """


rule sex_assignment:
    input:
        bam = "out/human/{sample}/{sample}.auto.sorted.bam"
    output:
        skglnd_sex = "out/human/{sample}/{sample}.skglnd.sex",
        mtnk_sex = "out/human/{sample}/{sample}.mtnk.sex"
    envmodules:
        "python2", "samtools", "R/3.5.1", "gcc" ## Note for Myriam #"miniconda/4.8.3"
    shell:
        """
        samtools index {input.bam}
        samtools view {input.bam} | python2 src/scripts/skoglund_xy.py > {output.skglnd_sex}
        samtools idxstats {input.bam} > out/human/{wildcards.sample}/{wildcards.sample}.auto.idxstats
        Rscript src/scripts/sex_determination_mittnik.R out/human/{wildcards.sample}/{wildcards.sample}.auto > {output.mtnk_sex}
        """

rule schmutzi:
    input:
        human_mito = "out/human/{sample}/{sample}.mito.sorted.bam",
        rCRS = config["rCRS"]
    output:
        md_bam = "out/human/{sample}/{sample}.mito.sorted.MD.bam",
        schmutzi = "out/human/{sample}/{sample}.schmutzi",
        fasta = "out/human/{sample}/{sample}.mito.fasta",
    envmodules:
        "schmutzi", "gcc", "samtools"
    shell:
        """
        samtools calmd -b {input.human_mito} {input.rCRS} > {output.md_bam}
        samtools index {output.md_bam}
        contDeam.pl --library double --length 10 --out {output.schmutzi} {output.md_bam}
        schmutzi.pl --notusepredC --uselength --ref {input.rCRS} {output.schmutzi}\
        /apps/schmutzi/20171024/alleleFreqMT/eurasian/freqs/ {output.md_bam}
        log2fasta -q 30 out/human/{wildcards.sample}/{wildcards.sample}.mito.schmutzi_final_endo.log > {output.fasta}
        """


rule Y_haplogroup:
    input:
        bam = "out/human/{sample}/{sample}.auto.sorted.bam",
        hg19 = config["hg19"]
    output:
        #out_dir = directory("out/human/{sample}/Y_haplogroup"),
        #y_hg = "out/human/{sample}/Y_haplogroup/{sample}.Y_chr.hg",
        y_vcf = "out/human/{sample}/Y_haplogroup/{sample}.y.vcf",
        y_haplo = directory("out/human/{sample}/Y_haplogroup/{sample}.yhaplo")
    envmodules:
        "samtools", "gcc", "yhaplo", "bcftools"
    shell:
        """
        bcftools mpileup -Q 30 -q 30 -f {input.hg19} {input.bam} -r chrY | bcftools call -c -o {output.y_vcf}
        callHaplogroups.py -i {output.y_vcf} -c -hp -ds -dsd -as -asd -o {output.y_haplo}
        """

        #python src/scripts/Yleaf/Yleaf.py -bam {input.bam} -pos src/scripts/Yleaf/Position_files/WGS_hg19.txt -out {output.out_dir} -r 1 -q 20 -b 90 
        #python src/scripts/Yleaf/predict_haplogroup.py -input {output.out_dir} -out {output.y_hg}



###Contamination checks 
#run angsd
rule chrX_contam:
    input:
        angsd_path = "/home/apps/angsd/0.918",
        bam = "out/human/{sample}/{sample}.auto.sorted.bam"

    output:
        angsd_R = "out/human/{sample}/angsd.chrX.R.contam",
        angsd_conta = "out/human/{sample}/angsd.chrX.contam"
    envmodules:
        "angsd", "R", "gcc"
    shell:
        """
        angsd -i {input.bam} -r chrX:5000000-154900000 -doCounts 1  -iCounts 1 -minMapQ 30 -minQ 30 -out out/human/{wildcards.sample}/angsdput

        Rscript {input.angsd_path}/R/contamination.R mapFile={input.angsd_path}/RES/chrX.unique.gz hapFile={input.angsd_path}/RES/HapMapChrX.gz countFile=out/human/{wildcards.sample}/angsdput.icnts.gz mc.cores=4 > {output.angsd_R}
        {input.angsd_path}/misc/contamination -a out/human/{wildcards.sample}/angsdput.icnts.gz -h /apps/angsd/0.918/RES/HapMapChrX.gz -d 2 -e 100 2 > {output.angsd_conta}
        """
#./angsd -i my.bam -r X:5000000-154900000 -doCounts 1  -iCounts 1 -minMapQ 30 -minQ 20
#do jackKnife in R
rule install_haplogrep:
    output:
        "src/scripts/haplogrep"
    shell:
        "cd src/scripts/; curl -sL haplogrep.now.sh | bash"
rule mito_haplogroup:
    input:
        rescaled_mito = "out/human/{sample}/{sample}.mito.sorted.rescaled.bam",
        rCRS = config["rCRS"]
    output:
        mito_vcf = "out/human/{sample}/{sample}.mito.sorted.rescaled.vcf",
        mito_hg = "out/human/{sample}/{sample}.mito.hg"
    envmodules:
        "samtools", "bcftools", "gcc"
    shell:
       """
       bcftools mpileup -B -Ou -d 1000 -q 30 -f {input.rCRS} {input.rescaled_mito} | bcftools call -mv --ploidy 1 -Ou -o {output.mito_vcf}

       ./src/scripts/haplogrep classify --in {output.mito_vcf} --format vcf --out {output.mito_hg} --extend-report
        """
##Phenotypic traits bed file
##rescaled auto hg19 bam for the vcf generation
#Then compare to the bed file
##//scratch/eurac/maixner/human/SNPs
#####################################################################
#---------------------------------------------------------#
#-------------Taxonomic Classifications------------------#
rule diamond_nr:
    input:
        merged = "out/qc_reads/{sample}/{sample}.rmdup.merged.qc.fastq.gz",
        reads1 = "out/qc_reads/{sample}/{sample}.R1.rmdup.unmerged.qc.fastq.gz",
        reads2 = "out/qc_reads/{sample}/{sample}.R2.rmdup.unmerged.qc.fastq.gz",
        diamond_db = config["diamond_db"]
    output:
        diamond_tab = "out/taxonomic_classifications/diamond/{sample}/{sample}.nr.tab",
        discarded_tab = "out/taxonomic_classifications/diamond/{sample}/{sample}.discarded.nr.tab"
    threads: lambda wildcards: threads_num[wildcards.sample]
    envmodules:
        "diamond", "gcc"
    resources:
        mem_mb= lambda wildcards, input, attempt: max(threads_num[wildcards.sample], 4) *10000 * attempt,
        disk_mb= lambda wildcards, input, attempt: max(threads_num[wildcards.sample], 4) *10000 * attempt
        #mem_mb=lambda wildcards, input, attempt: max(2*input.size_mb, 100000) * 0.25 * attempt
    shell:
        """
        diamond blastx -p {threads} -d {input.diamond_db} \
        -q {input.merged} -o {output.diamond_tab} -b 12 -c 1

        zcat {input.reads1} {input.reads2} | diamond blastx -p {threads} -d {input.diamond_db} \
        -o {output.discarded_tab} -b 12 -c 1
        """
rule MEGAN:
    input:
        diamond_tab = "out/taxonomic_classifications/diamond/{sample}/{sample}.nr.tab",
        discarded_tab = "out/taxonomic_classifications/diamond/{sample}/{sample}.discarded.nr.tab",
        megan_prot_db = config["megan_prot"]
    output:
        diamond_rma6 = "out/taxonomic_classifications/diamond/{sample}/{sample}.nr.rma6",
        discarded_rma6 = "out/taxonomic_classifications/diamond/{sample}/{sample}.discarded.nr.rma6",
    threads: lambda wildcards: threads_num[wildcards.sample]
    envmodules:
        "megan", "gcc"
    resources:
        mem_mb= lambda wildcards, input, attempt: max(threads_num[wildcards.sample], 4) *10000 * attempt,
        disk_mb= lambda wildcards, input, attempt: max(threads_num[wildcards.sample], 4) *10000 * attempt
    shell:
        """
        #unset DISPLAY
        blast2rma \
            --in {input.diamond_tab} \
            --format "BlastTab" \
            --blastMode "BlastX" \
            --out {output.diamond_rma6} \
            --minPercentIdentity 80 \
            --minSupportPercent 0.0001 \
            --threads {threads} \
            --acc2taxa {input.megan_prot_db}
        
		blast2rma \
            --in {input.discarded_tab} \
            --format "BlastTab" \
            --blastMode "BlastX" \
            --out {output.discarded_rma6} \
            --minPercentIdentity 80 \
            --minSupportPercent 0.0001 \
            --threads {threads} \
            --acc2taxa {input.megan_prot_db}
        """
rule metaphlan:
    input:
        merged = "out/qc_reads/{sample}/{sample}.rmdup.merged.qc.fastq.gz",
        reads1 = "out/qc_reads/{sample}/{sample}.R1.rmdup.unmerged.qc.fastq.gz",
        reads2 = "out/qc_reads/{sample}/{sample}.R2.rmdup.unmerged.qc.fastq.gz",
        db = config["metaphlan_db"]
    output:
        merged = "out/taxonomic_classifications/metaphlan/{sample}/{sample}.def.metaphlan",
        reads1 = "out/taxonomic_classifications/metaphlan/{sample}/{sample}.F.def.metaphlan",
        reads2 = "out/taxonomic_classifications/metaphlan/{sample}/{sample}.R.def.metaphlan",
        sp_table = "out/taxonomic_classifications/metaphlan/{sample}/merged_abundance_table_species.txt",
        heatmap = report("out/taxonomic_classifications/metaphlan/{sample}/{sample}.heatmap.pdf", category="MetaPhlAn")
    threads: lambda wildcards: threads_num[wildcards.sample]
    envmodules:
        "metaphlan", "python3/3.8.5", "gcc"
    resources:
        mem_mb= lambda wildcards, input, attempt: max(threads_num[wildcards.sample], 4) * 5000 * attempt,
        disk_mb= lambda wildcards, input, attempt: max(threads_num[wildcards.sample], 4) * 5000 * attempt
    shell:
        """
        metaphlan --bowtie2db {input.db} --input_type fastq --nproc {threads} --min_mapq_val 25 --no_map --add_viruses --read_min_len 25 {input.merged}  > {output.merged}

        metaphlan --bowtie2db {input.db} --input_type fastq --nproc {threads} --min_mapq_val 25 --no_map --add_viruses --read_min_len 25 {input.reads1}  > {output.reads1}

        metaphlan --bowtie2db {input.db} --input_type fastq --nproc {threads} --min_mapq_val 25 --no_map --add_viruses --read_min_len 25 {input.reads2}  > {output.reads2}

        merge_metaphlan_tables.py out/taxonomic_classifications/metaphlan/{wildcards.sample}/*def.metaphlan > out/taxonomic_classifications/metaphlan/{wildcards.sample}/merged_abundance_table.txt

        grep -E "s__|clade" out/taxonomic_classifications/metaphlan/{wildcards.sample}/merged_abundance_table.txt | sed 's/^.*s__//g' | cut -f1,3- | sed -e 's/clade_name/body_site/g' > {output.sp_table}
        
        python3 /apps/metaphlan/3.0.1/lib/python3.8/site-packages/metaphlan/utils/hclust2/hclust2.py \
        -i {output.sp_table} \
        -o {output.heatmap} \
        --f_dist_f correlation \
        --s_dist_f euclidean \
        --cell_aspect_ratio 0.5 \
        -l \
        --flabel_size 4 \
        --slabel_size 4 \
        --max_flabel_len 100 \
        --max_slabel_len 100 \
        --minv 0.1 \
        --dpi 300 --no_fclustering --no_sclustering
        """
rule kraken2:
    input:
        merged = "out/qc_reads/{sample}/{sample}.rmdup.merged.qc.fastq.gz",
        reads1 = "out/qc_reads/{sample}/{sample}.R1.rmdup.unmerged.qc.fastq.gz",
        reads2 = "out/qc_reads/{sample}/{sample}.R2.rmdup.unmerged.qc.fastq.gz",
        db = config["kraken_db"]
    output:
        out = "out/taxonomic_classifications/kraken/{sample}/{sample}.kraken.out",
        krona = "out/taxonomic_classifications/kraken/{sample}/{sample}.krona.html",
        report = "out/taxonomic_classifications/kraken/{sample}/{sample}.kraken.report",
        bracken = "out/taxonomic_classifications/kraken/{sample}/{sample}.kraken.bracken",
        out_dis = "out/taxonomic_classifications/kraken/{sample}/{sample}.discarded.kraken.out",
        krona_dis = "out/taxonomic_classifications/kraken/{sample}/{sample}.discarded.krona.html",
        report_dis = "out/taxonomic_classifications/kraken/{sample}/{sample}.discarded.kraken.report",
        bracken_dis = report("out/taxonomic_classifications/kraken/{sample}/{sample}.discarded.kraken.bracken", category="Kraken")
    threads: 
        8
    envmodules:
        "krona", "bracken", "kraken", "perl/5.26.1"
    resources:
        mem_mb=lambda wildcards, input, attempt: max(2*input.size_mb, 10000) ** 1 * attempt
    shell:
        """
        kraken2 --db {input.db} \
        --report {output.report} \
        --threads {threads} \
        --output {output.out} {input.merged}

        kraken2 --db {input.db} \
        --report {output.report_dis} \
        --threads {threads} \
        --output {output.out_dis} \
        --paired {input.reads1} {input.reads2}

        ktImportTaxonomy -q 2 -t 3 {output.out} -o {output.krona}
        rm -r out/taxonomic_classifications/kraken/{wildcards.sample}/{wildcards.sample}.krona.html.files

        ktImportTaxonomy -q 2 -t 3 {output.out_dis} -o {output.krona_dis}
        rm -r out/taxonomic_classifications/kraken/{wildcards.sample}/{wildcards.sample}.discarded.krona.html.files


        bracken -d {input.db} -i {output.report} -o {output.bracken} -l S

        bracken -d {input.db} -i {output.report_dis} -o {output.bracken_dis} -l S
        """
#‘mem_mb=max(2*input.size_mb, 1000)’ 
rule summary_stats:
    input:
        reads1 = "data/{sample}.R1.fastq.gz",
        reads_merged = "out/qc_reads/{sample}/{sample}.merged.qc.fastq.gz",
        reads1_unmerged = "out/qc_reads/{sample}/{sample}.R1.unmerged.qc.fastq.gz",
        dedup = "out/qc_reads/{sample}/{sample}.rmdup.merged.qc.fastq.gz",
        reads1_rmdup = "out/qc_reads/{sample}/{sample}.R1.rmdup.unmerged.qc.fastq.gz",
        reads2_rmdup = "out/qc_reads/{sample}/{sample}.R2.rmdup.unmerged.qc.fastq.gz",
        mito = "out/human/{sample}/{sample}.mito.sorted.bam",
        auto = "out/human/{sample}/{sample}.auto.sorted.bam",
        skglnd_sex = "out/human/{sample}/{sample}.skglnd.sex",
        mtnk_sex = "out/human/{sample}/{sample}.mtnk.sex",
        auto_cov = "out/human/{sample}/{sample}.auto.qualimap",
        mito_cov = "out/human/{sample}/{sample}.mito.qualimap"
    output:
        stats = "out/qc_reads/{sample}/{sample}.seq_stats.txt",
        plt = report("out/qc_reads/{sample}/{sample}.seq_stats.txt.pdf")
    envmodules:
        "samtools", "R"
    threads: 
        2
    shell:
        """
        echo "---------------------General_stats----------------------------" >> {output.stats}
        count=$(zgrep -c "^+$" {input.reads1})
        echo "Total number of raw reads "\t" $count" >> {output.stats}
        merge=$(zgrep -c "^+$" {input.reads_merged})
        echo "Total number of merged reads "\t" $merge" >> {output.stats}
        echo "Percentage of merged reads "\t" $(echo "$merge*100/$count" | bc -l) %" >> {output.stats}
        echo "" >>  {output.stats}
        echo "---------------------Deduplication_stats----------------------" >> {output.stats}
        dedup=$(zgrep -c "^+$" {input.dedup})
        echo "Total number of deduplicated reads "\t" $dedup" >> {output.stats}
        echo "Percentage of duplication "\t" $(echo "100-$(echo "$dedup*100/$merge" | bc -l)") %">> {output.stats}
        echo "Total number of deduplicated discarded forward reads "\t" " $(zgrep -c "^+$" {input.reads1_rmdup}) >> {output.stats}
        echo "Total number of deduplicated discarded reverse reads "\t" " $(zgrep -c "^+$" {input.reads2_rmdup}) >> {output.stats}
        echo "" >>  {output.stats}
        echo "---------------------Human_mitochondrial_DNA_stats------------" >> {output.stats}
        mito=$(samtools view -c -F 260 {input.mito} )
        echo "Number of mapped reads "\t" $mito" >> {output.stats}
        echo "Percentage of mitochondrial reads "\t" $(echo "100*$mito/$dedup" | bc -l) %" >> {output.stats}
        echo "Mitochondrial genome coverage "\t" $(grep "mean coverageData" {input.mito_cov}/genome_results.txt | cut -d"=" -f2)" >> {output.stats}
        echo $(grep "std coverageData" {input.mito_cov}/genome_results.txt)>> {output.stats}
        echo "" >>  {output.stats}
        echo "---------------------Human_autosomal_DNA_stats----------------" >> {output.stats}
        auto=$(samtools view -c -F 260 {input.auto} )
        echo "Number of mapped reads "\t" $auto" >> {output.stats}
        echo "Human endogenous DNA content "\t" $(echo "100*$auto/$dedup" | bc -l) %" >> {output.stats}
        echo "Human genome coverage "\t" $(grep "mean coverageData" {input.auto_cov}/genome_results.txt| cut -d"=" -f2)" >> {output.stats}
        echo $(grep "std coverageData" {input.auto_cov}/genome_results.txt)>> {output.stats}
        echo "" >>  {output.stats}
        echo "---------------------Human_sex_assignment----------------------" >> {output.stats}
        echo "Skoglund sex assignment: " >> {output.stats}
        head -2 {input.skglnd_sex} >> {output.stats}
        echo "" >>  {output.stats}
        echo "Mittnik sex assignment: " >> {output.stats}
        echo $(grep "Sex assignment" {input.mtnk_sex} | cut -c 6- | rev | cut -c 2- | rev) >> {output.stats}
        echo "" >>  {output.stats}
        Rscript src/scripts/PlotStats.R {output.stats}
        """
