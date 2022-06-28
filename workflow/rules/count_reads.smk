rule get_fastq_statistics:
    input:
        merged = expand(scratch_dict["kaiju"]["merged_reads"] / "{{taxon}}" / "{sample}_merged.fastq.gz", sample=SAMPLES),
        unmerged_R1 = expand(scratch_dict["kaiju"]["merged_reads"] / "{{taxon}}" / "{sample}_1_unmerged.fastq.gz", sample=SAMPLES),
        unmerged_R2 = expand(scratch_dict["kaiju"]["merged_reads"] / "{{taxon}}" / "{sample}_2_unmerged.fastq.gz", sample=SAMPLES),
        unpaired = expand(scratch_dict["kaiju"]["merged_reads"] / "{{taxon}}" / "{sample}_unpaired.fastq.gz", sample=SAMPLES),
    output:
        scratch_dict["read_stats"] / "{taxon}_read_stats.tsv"
    conda:
        "../envs/seqkit.yaml"
    params:
        ext = ".fastq.gz"
    shell:
        "seqkit stats $(dirname {input.merged[0]})/*{params.ext} -T > {output}"

rule make_table:
    input:
        taxon_read_stat_files=expand(scratch_dict["read_stats"] / "{taxon}_read_stats.tsv", taxon=list(config["taxons_of_interst"].keys()))
    output:
        read_count = results_dict['taxon_read_counts'],
        base_count = results_dict['taxon_base_counts']
    conda:
        "../envs/deep_variant_calling.yaml"
    params:
        ext = ".fastq.gz"
    script:
        "../scripts/taxon_read_count_table.py"
