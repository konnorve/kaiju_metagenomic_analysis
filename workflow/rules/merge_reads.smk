rule run_merge:
    input:
        r1 = scratch_dict["kaiju"]["taxon_assigned_reads"] / "{taxon}" / "{sample}_1_trimmed.fastq.gz",
        r2 = scratch_dict["kaiju"]["taxon_assigned_reads"] / "{taxon}" / "{sample}_2_trimmed.fastq.gz"
    output:
        merged = scratch_dict["kaiju"]["merged_reads"] / "{taxon}" / "{sample}_merged.fastq.gz",
        o1 = scratch_dict["kaiju"]["merged_reads"] / "{taxon}" / "{sample}_1_unmerged.fastq.gz",
        o2 = scratch_dict["kaiju"]["merged_reads"] / "{taxon}" / "{sample}_2_unmerged.fastq.gz",
        ihist= scratch_dict["kaiju"]["merged_reads_histogram"] / "{sample}_{taxon}_bbmerge_histogram_strict.txt",
    conda:
        "../envs/bbtools.yaml"
    log:
        "logs/merge_reads/run_merge/{sample}.{taxon}.log"
    benchmark:
        "benchmark/merge_reads/run_merge/{sample}.{taxon}.benchmark"
    shell:
        "bbmerge.sh "
        "in1={input.r1} "
        "in2={input.r2} "
        "out={output.merged} "
        "outu1={output.o1} "
        "outu2={output.o2} "
        "ihist={output.ihist} "
        "interleaved=f "
        "strict=t"

rule link_unpaired_to_merged_dir:
    input:
        scratch_dict["kaiju"]["taxon_assigned_reads"] / "{taxon}" / "{sample}_unpaired_trimmed.fastq.gz"
    output:
        scratch_dict["kaiju"]["merged_reads"] / "{taxon}" / "{sample}_unpaired.fastq.gz",
    shell:
        "ln {input} {output}; sleep 1"