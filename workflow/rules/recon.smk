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

# define final gapseq files
sampleTable['rxn_file'] = "models/" + SAMPLES + "/" + SAMPLES + "-" + config.get("search_term", 1) + "-Reactions.tbl.gz"
sampleTable['pwy_file'] = "models/" + SAMPLES + "/" + SAMPLES + "-" + config.get("search_term", 1) + "-Pathways.tbl.gz"
sampleTable['trs_file'] = "models/" + SAMPLES + "/" + SAMPLES + "-Transporter.tbl.gz"

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

def get_rxn_file_bygroup(wildcards):
    return " ".join(sampleTable[sampleTable['Group'] == int(wildcards.group)]['rxn_file'])

def get_pwy_file_bygroup(wildcards):
    return " ".join(sampleTable[sampleTable['Group'] == int(wildcards.group)]['pwy_file'])

def get_trs_file_bygroup(wildcards):
    return " ".join(sampleTable[sampleTable['Group'] == int(wildcards.group)]['trs_file'])



# Rules

localrules:
    gapseq_medium,
    all

rule all:
    input:
        expand("models/{sample}/{sample}.RDS", sample=SAMPLES),
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


rule gapseq_tar:
    output: temp("gapseq.tar.gz")
    params:
        gapseqdir=config.get("gapseq_dir",1),
        copygapseqdir=config.get("copy_gapseq_dir",1)
    resources:
        mem_mb=4000,
        time=1
    threads: 1
    shell:
        """
        if [ {params.copygapseqdir} == "True" ]; then
            currdir=`pwd`
            cd {params.gapseqdir}/../
            tar --exclude='.*' -czf $currdir/{output} gapseq/
        else
            touch {output}
        fi
        """

rule gapseq_find:
    input:
        genomes=lambda wildcards: sampleTable[sampleTable['Group'] == int(wildcards.group)]['genome_use'].tolist(),
        gapseqarchive=rules.gapseq_tar.output
    output: temp("group_{group}.txt")
    params:
        b=config.get("find_b", 1),
        splids=get_sample_bygroup,
        taxonomy=get_taxonomy_bygroup,
        rxnfiles=get_rxn_file_bygroup,
        pwyfiles=get_pwy_file_bygroup,
        trsfiles=get_trs_file_bygroup,
        gapseqdir=config.get("gapseq_dir",1),
        copygapseqdir=config.get("copy_gapseq_dir", 1),
        searchterm=config.get("search_term", 1)
    threads: config.get("find_threads", 1)
    resources:
        mem_mb=config.get("find_mem", 1) * 1000,
        time=config.get("find_time", 1),
        attempt=lambda wildcards, attempt: '_{}.log'.format(attempt)
    log: "logs/find/find_group_{group}"
    benchmark: "benchmark/find/group_{group}.tsv"
    shell:
        """
        mkdir -p models
        echo "Starting find group {wildcards.group}" > {log}{resources.attempt}
        if [ {params.copygapseqdir} == "True" ]; then
            #rsync -av -q --exclude=".*" {params.gapseqdir} {resources.tmpdir}/gapsnake_{wildcards.group} >> {log}{resources.attempt}
            mkdir -p {resources.tmpdir}/gapsnake_{wildcards.group}
            cp {input.gapseqarchive} {resources.tmpdir}/gapsnake_{wildcards.group}/{input.gapseqarchive}
            tar -xzf {resources.tmpdir}/gapsnake_{wildcards.group}/{input.gapseqarchive} -C {resources.tmpdir}/gapsnake_{wildcards.group}/
            rm {resources.tmpdir}/gapsnake_{wildcards.group}/{input.gapseqarchive}
            echo "{resources.tmpdir}/{wildcards.group}" >> {log}{resources.attempt}
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
            
            ncores=`echo "{threads}*3/$njobs"| bc`
            mkdir -p models/$spl
            ncores=$(( {threads} < ncores ? {threads} : ncores ))
            
            echo "ncores for find: $ncores" > logs/find/${{spl}}{resources.attempt}
            echo "gapseq: $gapseq" >> logs/find/${{spl}}{resources.attempt}
            
            # Reactions and Pathways
            date >> logs/find/${{spl}}{resources.attempt}
            echo "Starting reaction and pathway prediction" >> logs/find/${{spl}}{resources.attempt}
            $gapseq find -p {params.searchterm} -b {params.b} -t $tax -m $tax -K $ncores -v 0 -O -f models/$spl $geno >> logs/find/${{spl}}{resources.attempt}
            
            # Transporters
            date >> logs/find/${{spl}}{resources.attempt}
            echo "Starting transporter prediction" >> logs/find/${{spl}}{resources.attempt}
            $gapseq find-transport -b {params.b} -K $ncores -v 0 -f models/$spl $geno >> logs/find/${{spl}}{resources.attempt}
            
            # gzip all ".tbl" files
            date >> logs/find/${{spl}}{resources.attempt}
            echo "Compressing files" >> logs/find/${{spl}}{resources.attempt}
            gzip -f models/$spl/$spl-*.tbl >> logs/find/${{spl}}{resources.attempt}
            
        }}
        export -f gsfind
        
        genomes=({input.genomes})
        
        # check if really all samples still need to be processed, if not, give only indices that still need to be done to GNU parallel
        rxnf=({params.rxnfiles})
        pwyf=({params.pwyfiles})
        trsf=({params.trsfiles})
        spls=({params.splids})
        idsall=(`seq 0 $((${{#genomes[@]}} -1))`)
        idstodo=()
        
        for i in "${{!idsall[@]}}"; do
            if [[ ! -d "models/${{spls[i]}}" || ! -e "${{rxnf[i]}}" || ! -e "${{pwyf[i]}}" || ! -e "${{trsf[i]}}" ]]; then
                idstodo+=("${{idsall[i]}}")
            fi
        done
        
        k=${{#idstodo[@]}}
        n=${{#idsall[@]}}
        echo "Remaining genomes in group({wildcards.group}): $k / $n" >> {log}{resources.attempt}
        
        # finally parallel processing of remaining samples
        nthreads={threads}
        njobs=$(( nthreads < k ? nthreads : k ))
        export njobs
        parallel --env njobs --jobs $njobs gsfind ::: ${{idstodo[@]}}  >> {log}{resources.attempt}
        
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
    params:
        biomass=get_biomass,
        u=config.get("draft_u", 1),
        l=config.get("draft_l", 1),
        searchterm=config.get("search_term", 1)
    output:
        draft="models/{sample}/{sample}-draft.RDS",
        rxnWeights="models/{sample}/{sample}-rxnWeights.RDS",
        rxnXgenes="models/{sample}/{sample}-rxnXgenes.RDS",
        xml="models/{sample}/{sample}-draft.xml.gz"
    threads: config.get("draft_threads", 1)
    resources:
        mem_mb=config.get("draft_mem", 1) * 1000,
        time=config.get("draft_time", 1)
    log:
        "logs/draft/{sample}.log"
    shell:
        """
        gapseq draft -r models/{wildcards.sample}/{wildcards.sample}-{params.searchterm}-Reactions.tbl.gz -t models/{wildcards.sample}/{wildcards.sample}-Transporter.tbl.gz -b {params.biomass} -c {input.genome} -p models/{wildcards.sample}/{wildcards.sample}-{params.searchterm}-Pathways.tbl.gz -u {params.u} -l {params.l} -f models/{wildcards.sample} > {log}
        gzip -f models/{wildcards.sample}/{wildcards.sample}-draft.xml
        """

rule gapseq_medium:
    input:
        model="models/{sample}/{sample}-draft.RDS",
        grpstamp=get_group_file,
    params:
        c=config.get("medium_c", 1),
        searchterm=config.get("search_term", 1)
    output:
        "models/{sample}/{sample}-medium.csv"
    log:
        "logs/medium/{sample}.log"
    shell:
        """
        par_c={params.c}
        if [ -n "$par_c" ]; then
            gapseq medium -m {input.model} -p models/{wildcards.sample}/{wildcards.sample}-{params.searchterm}-Pathways.tbl.gz -c {params.c} -f models/{wildcards.sample} > {log}
        else
            gapseq medium -m {input.model} -p models/{wildcards.sample}/{wildcards.sample}-{params.searchterm}-Pathways.tbl.gz -f models/{wildcards.sample} > {log}
        fi
        """
        
rule gapseq_fill:
    input:
        draft="models/{sample}/{sample}-draft.RDS",
        medium=get_medium,
        rxnWeights="models/{sample}/{sample}-rxnWeights.RDS",
        rxnXgenes="models/{sample}/{sample}-rxnXgenes.RDS"
    params:
        b=config.get("fill_b", 1),
        mingr=config.get("fill_mingr", 1)
    output:
        model="models/{sample}/{sample}.RDS",
        xml="models/{sample}/{sample}.xml.gz"
    threads: config.get("fill_threads", 1)
    resources:
        mem_mb=config.get("fill_mem", 1) * 1000,
        time=config.get("fill_time", 1)
    log:
        "logs/fill/{sample}.log"
    shell:
        """
        if grep -q cpd11640 "{input.medium}"; then
            gapseq fill -m {input.draft} -n {input.medium} -c {input.rxnWeights} -g {input.rxnXgenes} -b {params.b} -e highH2 -f models/{wildcards.sample} -k {params.mingr} > {log}
        else
             gapseq fill -m {input.draft} -n {input.medium} -c {input.rxnWeights} -g {input.rxnXgenes} -b {params.b} -f models/{wildcards.sample} -k {params.mingr} > {log}
        fi
        
        gzip -f models/{wildcards.sample}/{wildcards.sample}.xml
        """
