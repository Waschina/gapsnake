import pandas as pd
import numpy as np

configfile: "config.yaml"

sampleTable = pd.read_csv("samples.tsv", index_col="sample", sep="\t")
SAMPLES = sampleTable.index.values

# define genome after optional translation
if config.get("translate_cds", 1):
    sampleTable['genome_use'] = "genomes_prot/" + SAMPLES + ".faa.gz"
else:
    sampleTable['genome_use'] = sampleTable['genome_file']

# assign groups for bundle processing
sampleTable['Group'] = list(pd.Series(range(len(sampleTable))) // int(config.get("find_batch_size", 1)))
ngroups=int(np.ceil(len(sampleTable) / int(config.get("find_batch_size", 1))))
GROUPS = np.unique(list(sampleTable['Group']))


#print(sampleTable)
print(ngroups)

def get_genome_orig(wildcards):
    return sampleTable.at[wildcards.sample, 'genome_file']

def get_genome(wildcards):
    return sampleTable.at[wildcards.sample, 'genome_use']

def get_taxonomy(wildcards):
    return sampleTable.at[wildcards.sample, 'taxonomy']

def get_taxonomy_bygroup(wildcards):
    return " ".join(sampleTable[sampleTable['Group'] == int(wildcards.group)]['taxonomy'])

def get_sample_bygroup(wildcards):
    return " ".join(sampleTable[sampleTable['Group'] == int(wildcards.group)].index)

def get_biomass(wildcards):
    return sampleTable.at[wildcards.sample, 'biomass']

def get_medium(wildcards):
    return sampleTable.at[wildcards.sample, 'medium']
    
def get_translation_table(wildcards):
    return str(sampleTable.at[wildcards.sample, 'translation_table'])

def get_group_file(wildcards):
    return "group_" + str(sampleTable.at[wildcards.sample, 'Group']) + ".txt"


# Rules

rule all:
    input:
        expand("models/{sample}/{sample}.txt", sample=SAMPLES),
        expand("group_{group}.txt", group=GROUPS),
    output:
        touch("finished_recon")

rule pyrodigal:
    input:
        genome=get_genome_orig
    params:
        translation_table=get_translation_table
    output:
        prot="genomes_prot/{sample}.faa.gz"
    threads: 1
    resources:
        mem_mb=config.get("translate_mem", 1) * 1000,
        time=config.get("translate_time", 1)
    conda:
        "../envs/pyrodigal.yaml"
    script:
        "../scripts/pyrodigal_predict.py"


rule gapseq_find:
    input:
        genomes=lambda wildcards: sampleTable[sampleTable['Group'] == int(wildcards.group)]['genome_use'].tolist()
    output: temp("group_{group}.txt")
    params:
        b=config.get("find_b", 1),
        splids=get_sample_bygroup,
        taxonomy=get_taxonomy_bygroup,
        gapseqdir=config.get("gapseq_dir",1),
        copygapseqdir=config.get("copy_gapseq_dir", 1),
        searchterm=config.get("search_term", 1)
    threads: config.get("find_threads", 1)
    resources:
        mem_mb=config.get("find_mem", 1) * 1000,
        time=config.get("find_time", 1)
    log: "logs/find/group_{group}.log"
    benchmark: "benchmark/find/group_{group}.tsv"
    shell:
        """
        if [ {params.copygapseqdir} == "True" ]; then
            rsync -av -q --exclude=".*" {params.gapseqdir} {resources.tmpdir}/gapsnake_{wildcards.group} > {log}
            echo "{resources.tmpdir}/{wildcards.group}" >> {log}
        fi
        
        gsfind() {{
            # redirect gapseq if gapseq dir is copied to tmp-dir
            if [ {params.copygapseqdir} == "True" ]; then
                gapseq="{resources.tmpdir}/gapsnake_{wildcards.group}/gapseq/./gapseq"
            else
                gapseq="{params.gapseqdir}/./gapseq"
            fi
            
            idx=$1
            genomes=({input.genomes})
            geno=${{genomes[$idx]}}
            taxonomy=({params.taxonomy})
            tax=${{taxonomy[$idx]}}
            splids=({params.splids})
            spl=${{splids[$idx]}}
            
            (( ncores=({threads}+3-1)/3 ));
            
            # Reactions and Pathways
            date > logs/find/$spl.log
            echo "Starting reaction and pathway prediction\n" >> logs/find/$spl.log
            $gapseq find -p {params.searchterm} -b {params.b} -t $tax -m $tax -K $ncores -v 0 -O -f models/$spl $geno >> logs/find/$spl.log 
            
            # Transporters
            date >> logs/find/$spl.log
            echo "Starting transporter prediction\n" >> logs/find/$spl.log
            $gapseq find-transport -b {params.b} -K $ncores -v 0 -f models/$spl $geno >> logs/find/$spl.log
            
            # gzip all ".tbl" files
            date >> logs/find/$spl.log
            echo "Compressing files\n" >> logs/find/$spl.log
            gzip -f models/$spl/$spl-*.tbl >> logs/find/$spl.log
            
        }}
        export -f gsfind
        
        genomes=({input.genomes})
        
        parallel --jobs {threads} gsfind ::: $(seq 0 $((${{#genomes[@]}} - 1))) > {log}
        
        # check if everything is there
        splids=({params.splids})
        output_check=true
        
        for sample in "${{splids[@]}}"; do
            # Define the required file paths
            files=(
                "models/${{sample}}/${{sample}}-{params.searchterm}-Reactions.tbl.gz"
                "models/${{sample}}/${{sample}}-{params.searchterm}-Pathways.tbl.gz"
                "models/${{sample}}/${{sample}}-Transporter.tbl.gz"
            )
    
            # Check if each file exists
            for file in "${{files[@]}}"; do
                if [[ ! -f $file ]]; then
                    echo "Missing file: $file"
                    all_files_exist=false  # Set flag to false if a file is missing
                fi
            done
        done
        
        # cleaning up
        if [ {params.copygapseqdir} == "True" ]; then
            rm -r {resources.tmpdir}/gapsnake_{wildcards.group}
        fi
        
        # finalizing
        if [[ $output_check ]]; then
            echo "Group  reconstructions succeeded." > {output}
        else
            exit 1
        fi
        
        """

rule gapseq_draft:
    input:
        grpstamp=get_group_file,
        genome=get_genome
    output: touch("models/{sample}/{sample}.txt")


