rule assign_paired_reads:
    input:
        kaiju_db_index = Path(config["input"]["kaiju_db_index"]),
        kaiju_db_nodes = Path(config["input"]["kaiju_db_nodes"]),
        R1 = scratch_dict["trimmed_reads"] / "{sample}_1_trimmed.fastq.gz",
        R2 = scratch_dict["trimmed_reads"] / "{sample}_2_trimmed.fastq.gz",
    output:
        scratch_dict["kaiju"]["kaiju_assignment"] / "{sample}_paired_kaiju.out"
    conda:
        "../envs/kaiju.yaml"
    log:
        "logs/kaiju/assign_paired_reads/{sample}.log"
    benchmark:
        "benchmark/kaiju/assign_paired_reads/{sample}.benchmark"
    threads: 10
    resources:
        mem_mb = 100000
    shell:
        """
        kaiju -z {threads} \
            -t {input.kaiju_db_nodes} \
            -f {input.kaiju_db_index} \
            -i {input.R1} \
            -j {input.R2} \
            -o {output} 2> {log}
        """


rule assign_unpaired_reads:
    input:
        kaiju_db_index = Path(config["input"]["kaiju_db_index"]),
        kaiju_db_nodes = Path(config["input"]["kaiju_db_nodes"]),
        R = scratch_dict["trimmed_reads"] / "{sample}_unpaired_trimmed.fastq.gz",
    output:
        scratch_dict["kaiju"]["kaiju_assignment"] / "{sample}_unpaired_kaiju.out"
    conda:
        "../envs/kaiju.yaml"
    log:
        "logs/kaiju/assign_unpaired_reads/{sample}.log"
    benchmark:
        "benchmark/kaiju/assign_unpaired_reads/{sample}.benchmark"
    threads: 10
    resources:
        mem_mb = 100000
    shell:
        """
        kaiju -z {threads} \
            -t {input.kaiju_db_nodes} \
            -f {input.kaiju_db_index} \
            -i {input.R} \
            -o {output} 2> {log}
        """


rule add_taxons_paired_reads:
    input:
        kaiju_db_nodes = Path(config["input"]["kaiju_db_nodes"]),
        kaiju_db_names = Path(config["input"]["kaiju_db_names"]),
        kaiju_assignment = scratch_dict["kaiju"]["kaiju_assignment"] / "{sample}_paired_kaiju.out"
    output:
        scratch_dict["kaiju"]["taxon_assignment"] / "{sample}_paired_kaiju.out"
    conda:
        "../envs/kaiju.yaml"
    log:
        "logs/kaiju/add_taxons_paired_reads/{sample}.log"
    benchmark:
        "benchmark/kaiju/add_taxons_paired_reads/{sample}.benchmark"
    shell:
        """
        kaiju-addTaxonNames -u -p \
            -t {input.kaiju_db_nodes} \
            -n {input.kaiju_db_names} \
            -i {input.kaiju_assignment} \
            -o {output} 2> {log}
        """


rule add_taxons_unpaired_reads:
    input:
        kaiju_db_nodes = Path(config["input"]["kaiju_db_nodes"]),
        kaiju_db_names = Path(config["input"]["kaiju_db_names"]),
        kaiju_assignment = scratch_dict["kaiju"]["kaiju_assignment"] / "{sample}_unpaired_kaiju.out"
    output:
        scratch_dict["kaiju"]["taxon_assignment"] / "{sample}_unpaired_kaiju.out"
    conda:
        "../envs/kaiju.yaml"
    log:
        "logs/kaiju/add_taxons_unpaired_reads/{sample}.log"
    benchmark:
        "benchmark/kaiju/add_taxons_unpaired_reads/{sample}.benchmark"
    shell:
        """
        kaiju-addTaxonNames -u -p \
            -t {input.kaiju_db_nodes} \
            -n {input.kaiju_db_names} \
            -i {input.kaiju_assignment} \
            -o {output} 2> {log}
        """

rule filter_taxons:
    input:
        scratch_dict["kaiju"]["taxon_assignment"] / "{sample}_{readpairing}_kaiju.out"
    output:
        scratch_dict["kaiju"]["filtering_assignments"] / "{taxon}" / "{sample}_{readpairing}_kaiju.out"
    log:
        "logs/kaiju/filter_taxons/{taxon}.{sample}.{readpairing}.log"
    benchmark:
        "benchmark/kaiju/filter_taxons/{taxon}.{sample}.{readpairing}.benchmark"
    params:
        filtr=lambda wildcards, output: config["taxons_of_interst"][wildcards.taxon]
    shell:
        """
        echo "{params.filtr}" > {log}
        cat {input} | grep "{params.filtr}" > {output}
        """

rule list_taxon_assigned:
    input:
        scratch_dict["kaiju"]["filtering_assignments"] / "{taxon}" / "{sample}_{readpairing}_kaiju.out"
    output:
        scratch_dict["kaiju"]["filtering_assignments"] / "{taxon}" / "{sample}_{readpairing}_kaiju.lst"
    log:
        "logs/kaiju/list_taxon_assigned/{taxon}.{sample}.{readpairing}.log"
    benchmark:
        "benchmark/kaiju/list_taxon_assigned/{taxon}.{sample}.{readpairing}.benchmark"
    shell:
        "cat {input} | awk '{{print $2}}' > {output} 2> {log}"

rule extract_fwd_reads:
    input:
        assignment = scratch_dict["kaiju"]["filtering_assignments"] / "{taxon}" / "{sample}_paired_kaiju.lst",
        read_group = scratch_dict["trimmed_reads"] / "{sample}_1_trimmed.fastq.gz",
    output:
        scratch_dict["kaiju"]["taxon_assigned_reads"] / "{taxon}" / "{sample}_1_trimmed.fastq",
    log:
        "logs/kaiju/extract_reads/{taxon}.{sample}.fwd.log"
    benchmark:
        "benchmark/kaiju/extract_reads/{taxon}.{sample}.fwd.benchmark"
    conda:
        "../envs/seqtk.yaml"
    shell:
        "seqtk subseq {input.read_group} {input.assignment} > {output} 2> {log}"

rule extract_rev_reads:
    input:
        assignment = scratch_dict["kaiju"]["filtering_assignments"] / "{taxon}" / "{sample}_paired_kaiju.lst",
        read_group = scratch_dict["trimmed_reads"] / "{sample}_2_trimmed.fastq.gz",
    output:
        scratch_dict["kaiju"]["taxon_assigned_reads"] / "{taxon}" / "{sample}_2_trimmed.fastq",
    log:
        "logs/kaiju/extract_reads/{taxon}.{sample}.rev.log"
    benchmark:
        "benchmark/kaiju/extract_reads/{taxon}.{sample}.rev.benchmark"
    conda:
        "../envs/seqtk.yaml"
    shell:
        "seqtk subseq {input.read_group} {input.assignment} > {output} 2> {log}"

rule extract_unpaired_reads:
    input:
        assignment = scratch_dict["kaiju"]["filtering_assignments"] / "{taxon}" / "{sample}_paired_kaiju.lst",
        read_group = scratch_dict["trimmed_reads"] / "{sample}_unpaired_trimmed.fastq.gz",
    output:
        scratch_dict["kaiju"]["taxon_assigned_reads"] / "{taxon}" / "{sample}_unpaired_trimmed.fastq",
    log:
        "logs/kaiju/extract_reads/{taxon}.{sample}.unpaired.log"
    log:
        "logs/kaiju/extract_reads/{taxon}.{sample}.unpaired.log"
    conda:
        "../envs/seqtk.yaml"
    shell:
        "seqtk subseq {input.read_group} {input.assignment} > {output} 2> {log}"

rule zip:
    input:
        scratch_dict["kaiju"]["taxon_assigned_reads"] / "{taxon}" / "{anything}.fastq",
    output:
        scratch_dict["kaiju"]["taxon_assigned_reads"] / "{taxon}" / "{anything}.fastq.gz",
    resources:
        mem_mb=100000,
    conda:
        "../envs/gzip.yaml"
    shell:
        "gzip {input}"
