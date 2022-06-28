
rule convert_sam2bam:
    input:
        scratch_dict["read_mapping"] / "{filter}" / "{sample}_mapped.sam",
    output:
        temp(scratch_dict["read_mapping"] / "{filter}" / "{sample}_mapped.bam"),
    conda:
        "../envs/samtools.yaml"
    log:
        "logs/samtools/convert_sam2bam/{sample}.{filter}.log"
    shell:
        "samtools view -S -b {input} > {output} 2> {log}"

rule sort_bam:
    input:
        scratch_dict["read_mapping"] / "{filter}" / "{sample}_mapped.bam",
    output:
        scratch_dict["read_mapping"] / "{filter}" / "{sample}_mapped_sorted.bam"
    conda:
        "../envs/samtools.yaml"
    log:
        "logs/samtools/sort_bam/{sample}.{filter}.log"
    shell:
        "samtools sort {input} -o {output} > {log}"

rule bam_coverage:
    input:
        scratch_dict["read_mapping"] / "{filter}" / "{sample}_mapped_sorted.bam"
    output:
        scratch_dict["read_mapping"] / "{filter}" / "{sample}_depth.tsv"
    conda:
        "../envs/samtools.yaml"
    log:
        "logs/remove_PCR_duplicates/samtools_depth/{sample}.{filter}.log"
    shell:
        "samtools depth -o {output} {input}"


rule plot_bam_depths:
    input:
        scratch_dict["read_mapping"] / "{filter}" / "{sample}_depth.tsv"
    output:
        depth_histogram = results_dict["depth_histograms"] / "{filter}" / "{sample}_depth_histogram.png",
        depth_genome = results_dict["depth_genome"] / "{filter}" / "{sample}_depth_genome.png",
    conda:
        "../envs/deep_variant_calling.yaml"
    script:
        "../scripts/depth_histogram.py"

# rule index_bam:
#     input:
#         scratch_dict["read_mapping"] / "{filter}" / "{sample}_mapped_sorted.bam",
#     output:
#         temp(scratch_dict["read_mapping"] / "{filter}" / "{sample}_mapped_sorted.bam.bai"),
#     resources:
#         mem_mb=100000,
#     conda:
#         "../envs/samtools.yaml"
#     log:
#         "logs/samtools/index_bam/{sample}.log"
#     shell:
#         "samtools index -b {input} > {log}"


# rule index_fasta:
#     input:
#         reference_genome_file,
#     output:
#         reference_genome_index_file,
#     resources:
#         mem_mb=100000,
#     conda:
#         "../envs/samtools.yaml"
#     shell:
#         "samtools faidx {input}"
    