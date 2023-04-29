import pandas as pd

configfile: "config.yaml"

sampleTable = pd.read_csv("samples.tsv", index_col="sample", sep="\t")
SAMPLES = sampleTable.index.values

def get_genome_nucl(wildcards):
    return sampleTable.at[wildcards.sample, 'genome_file']

def get_genome(wildcards):
    if config.get("translate_cds", 1):
        return 'genomes_prot/'+wildcards.sample+'.faa.gz'
    else:
        return get_genome_nucl(wildcards)

def get_taxonomy(wildcards):
    return sampleTable.at[wildcards.sample, 'taxonomy']

def get_biomass(wildcards):
    return sampleTable.at[wildcards.sample, 'biomass']

def get_medium(wildcards):
    return sampleTable.at[wildcards.sample, 'medium']


localrules:
    gapseq_medium,
    all


rule all:
    input:
        expand("models/{sample}/{sample}.RDS", sample=SAMPLES),
    output:
        touch("finished_recon")


rule pyrodigal:
    input:
        genome=get_genome_nucl
    output:
        prot="genomes_prot/{sample}.faa.gz"
    threads: 1
    resources:
        mem=config.get("translate_mem", 1),
        time=config.get("translate_time", 1)
    conda:
        "../envs/pyrodigal.yaml"
    script:
        "../scripts/pyrodigal_predict.py"


rule gapseq_find:
    input:
        get_genome
    params:
        b=config.get("find_b", 1),
        taxonomy=get_taxonomy
    output:
        rxn="models/{sample}/{sample}-all-Reactions.tbl.gz",
        pwy="models/{sample}/{sample}-all-Pathways.tbl.gz"
    threads: config.get("find_threads", 1)
    resources:
        mem=config.get("find_mem", 1),
        time=config.get("find_time", 1)
    log:
        "logs/find/{sample}.log"
    shell:
        """
        gapseq find -p all -b {params.b} -t {params.taxonomy} -m {params.taxonomy} -K {threads} -v 0 -O -f models/{wildcards.sample} {input} > {log}
        gzip -f models/{wildcards.sample}/{wildcards.sample}-all-Reactions.tbl
        gzip -f models/{wildcards.sample}/{wildcards.sample}-all-Pathways.tbl
        """
        
rule gapseq_find_transport:
    input:
        get_genome
    output:
        "models/{sample}/{sample}-Transporter.tbl.gz"
    threads: config.get("transport_threads", 1)
    resources:
        mem=config.get("transport_mem", 1),
        time=config.get("transport_time", 1)
    log:
        "logs/transport/{sample}.log"
    shell:
        """
        gapseq find-transport -b 200 -K {threads} -v 0 -f models/{wildcards.sample} {input} > {log}
        gzip -f models/{wildcards.sample}/{wildcards.sample}-Transporter.tbl
        """
        
rule gapseq_draft:
    input:
        rxn="models/{sample}/{sample}-all-Reactions.tbl.gz",
        pwy="models/{sample}/{sample}-all-Pathways.tbl.gz",
        trsp="models/{sample}/{sample}-Transporter.tbl.gz",
        genome=get_genome,
    params:
        biomass=get_biomass,
        u=config.get("draft_u", 1),
        l=config.get("draft_l", 1)
    output:
        draft="models/{sample}/{sample}-draft.RDS",
        rxnWeights="models/{sample}/{sample}-rxnWeights.RDS",
        rxnXgenes="models/{sample}/{sample}-rxnXgenes.RDS",
        xml="models/{sample}/{sample}-draft.xml.gz"
    threads: config.get("draft_threads", 1)
    resources:
        mem=config.get("draft_mem", 1),
        time=config.get("draft_time", 1)
    log:
        "logs/draft/{sample}.log"
    shell:
        """
        gapseq draft -r {input.rxn} -t {input.trsp} -b {params.biomass} -c {input.genome} -p {input.pwy} -u {params.u} -l {params.l} -f models/{wildcards.sample} > {log}
        gzip -f models/{wildcards.sample}/{wildcards.sample}-draft.xml
        """
        
rule gapseq_medium:
    input:
        model="models/{sample}/{sample}-draft.RDS",
        pwy="models/{sample}/{sample}-all-Pathways.tbl.gz"
    params:
        c=config.get("medium_c", 1)
    output:
        "models/{sample}/{sample}-medium.csv"
    log:
        "logs/medium/{sample}.log"
    shell:
        """
        par_c={params.c}
        if [ -n "$par_c" ]; then
            gapseq medium -m {input.model} -p {input.pwy} -c {params.c} -f models/{wildcards.sample} > {log}
        else
            gapseq medium -m {input.model} -p {input.pwy} -f models/{wildcards.sample} > {log}
        fi
        """
        
rule gapseq_fill:
    input:
        draft="models/{sample}/{sample}-draft.RDS",
        medium=get_medium,
        rxnWeights="models/{sample}/{sample}-rxnWeights.RDS",
        rxnXgenes="models/{sample}/{sample}-rxnXgenes.RDS"
    params:
        b=config.get("fill_b", 1)
    output:
        model="models/{sample}/{sample}.RDS",
        xml="models/{sample}/{sample}.xml.gz"
    threads: config.get("fill_threads", 1)
    resources:
        mem=config.get("fill_mem", 1),
        time=config.get("fill_time", 1)
    log:
        "logs/fill/{sample}.log"
    shell:
        """
        if grep -q cpd11640 "{input.medium}"; then
            gapseq fill -m {input.draft} -n {input.medium} -c {input.rxnWeights} -g {input.rxnXgenes} -b {params.b} -e highH2 -f models/{wildcards.sample} > {log}
        else
             gapseq fill -m {input.draft} -n {input.medium} -c {input.rxnWeights} -g {input.rxnXgenes} -b {params.b} -f models/{wildcards.sample} > {log}
        fi
        
        gzip -f models/{wildcards.sample}/{wildcards.sample}.xml
        """

