##############################################
# Snakefile – Bac_pipeline_AC/input_from_prokka
# Flow: Prokka → geNomad → Abricate → PADLOC → CRISPRidentify → Defense-finder
# Activate existing conda environments SIN --use-conda
##############################################

##############################################
# Instructions
# Install all tools in separate environments tailored to each one, including snakemake
# If the environment name differs from the one listed in this file, make the necessary adjustments in each case
# Environment locations are at your discretion; adjust them as your needs
# Output file locations for each tool are at your discretion; adjust them as your needs
# Output folders for each tool are stored in separate directories named after the tool followed by the identifier of the strain being analyzed
# All these folders are contained within a mother directory named after the strain identifier
# The environment and output locations in this document are examples and should be modified to suit your needs
##############################################

# To run the code, go to the location of the Snakefile and run by
# snakemake --cores 16
# Adjust the number of cores as the capacity of your computer

import os, glob

# ====== Base routes (examples)======
BASE     = "/root/proyectos/analisis_genomico/Bac_pipeline_AC/input_from_prokka"
DATA     = f"{BASE}/Data_inputs/Simon"
ANALYSIS = f"{BASE}/Analysis"
DBROOT   = f"{BASE}/DB"
ENVS     = "/root/proyectos/analisis_genomico/Bac_pipeline_AC/input_from_prokka/envs"

# ====== Conda environments (existing paths) ======
PROKKA_ENV         = f"{ENVS}/prokka_env"
GENOMAD_ENV        = f"{ENVS}/genomad_env"
ABRICATE_ENV       = f"{ENVS}/abricate_env"
PADLOC_ENV         = f"{ENVS}/padloc_env"
DEFENSEF_ENV       = f"{ENVS}/defensefinder_env_new"  # ajusta si tu env se llama distinto
CRISPRIDENTIFY_ENV = "/root/miniconda3/envs/crispridentify_py38"
CCTYPER_ENV        = f"{ENVS}/cctyper"      # AJUSTA al nombre real si difiere
ECTYPER_ENV        = f"{ENVS}/ectyper_env"      # AJUSTA
RESFINDER_ENV      = f"{ENVS}/resfinder_env"    # AJUSTA
MLST_ENV           = f"{ENVS}/mlst_env"

# ====== DBs / markers ======
GENOMAD_DB     = f"{DBROOT}/genomad_db"
PADLOC_DB_DIR  = f"{DBROOT}/padloc_db"
PADLOC_DB_MARK = f"{PADLOC_DB_DIR}/version.txt"
ABRICATE_READY = f"{DBROOT}/.abricate_db_ready"
DEFENSEF_DB_DIR  = f"{DBROOT}/defensefinder_db"
DEFENSEF_DB_MARK = f"{DEFENSEF_DB_DIR}/version.txt"
RESFINDER_DB_DIR   = "/root/resfinder_dbs/resfinder_db"
POINTFINDER_DB_DIR = "/root/resfinder_dbs/pointfinder_db"
RESFINDER_DB_MARK  = "/root/resfinder_dbs/.resfinder_pointfinder_db_ready"
ECTYPER_MASH_SKETCH = "/root/resfinder_dbs/ectyper_db/EnteroRef_GTDBSketch_20231003_V2.msh"
ECTYPER_MASH_MARK   = "/root/resfinder_dbs/ectyper_db/.ectyper_mash_ready"


# ====== Discovering entries ======
def collect_inputs(root):
    hits = {}
    for ext in (".fna", ".fa", ".fasta", ".fas"):
        for f in glob.glob(f"{root}/LCUS_*/**/*{ext}", recursive=True):
            # Take the LCUS_XXXX segment (folder) from the path. 
            # LCUS = identifier for the bacterial collection of Louis-Charles Fortier's laboratory at the University of Sherbrooke; XXXX = strain identification number.
            parts = os.path.normpath(f).split(os.sep)
            sample = next((p for p in parts if p.startswith("LCUS_")), None)
            if sample:
                hits.setdefault(sample, os.path.relpath(f, root))
    return hits

dict_data = collect_inputs(DATA)
samples   = sorted(dict_data)

if not samples:
    raise ValueError(f"No FASTA was found in {DATA}/LCUS_*/**/*.fna|fa|fasta|fas")

# ===================== Objective =====================
rule all:
    input:
        # per sample (example)
        expand(f"{ANALYSIS}" + "/{sample}/genome_analysis/prokka_{sample}/{sample}.fna",           sample=samples),
        expand(f"{ANALYSIS}" + "/{sample}/genome_analysis/geNomad_{sample}",                       sample=samples),
        expand(f"{ANALYSIS}" + "/{sample}/genome_analysis/abricate_{sample}/{sample}_summary.tab", sample=samples),
        expand(f"{ANALYSIS}" + "/{sample}/genome_analysis/padloc_{sample}",                        sample=samples),
        expand(f"{ANALYSIS}" + "/{sample}/genome_analysis/crispridentify_{sample}",                sample=samples),
        expand(f"{ANALYSIS}" + "/{sample}/genome_analysis/cctyper_{sample}",                       sample=samples),
        expand(f"{ANALYSIS}" + "/{sample}/genome_analysis/ECTyper_{sample}",                       sample=samples),
        expand(f"{ANALYSIS}" + "/{sample}/genome_analysis/resfinder_pointfinder_{sample}",         sample=samples),
        expand(f"{ANALYSIS}" + "/{sample}/genome_analysis/mlst_{sample}",                          sample=samples),

        # DBs
        f"{GENOMAD_DB}/version.txt",
        ABRICATE_READY,
        PADLOC_DB_MARK,
        DEFENSEF_DB_MARK,
        RESFINDER_DB_MARK,
        ECTYPER_MASH_MARK,

# ===================== PROKKA =====================
rule prokka:
    input:
        lambda w: f"{DATA}/{dict_data[w.sample]}"
    output:
        all_prokka = directory(f"{ANALYSIS}" + "/{sample}/genome_analysis/prokka_{sample}"),
        fna        = f"{ANALYSIS}" + "/{sample}/genome_analysis/prokka_{sample}/{sample}.fna",
        faa        = f"{ANALYSIS}" + "/{sample}/genome_analysis/prokka_{sample}/{sample}.faa"
    log:
        f"{ANALYSIS}" + "/{sample}/log/annotation/prokka_{sample}.log"
    message: " Prokka: {wildcards.sample}"
    shell:
        r"""
        source /root/miniconda3/etc/profile.d/conda.sh
        conda activate {PROKKA_ENV}
        prokka --outdir {output.all_prokka} --force --prefix {wildcards.sample} \
               --addgenes --locustag {wildcards.sample} {input} > {log} 2>&1
        """

# ===================== geNomad DB =====================
rule genomad_db:
    output:
        f"{GENOMAD_DB}/version.txt"
    message: "Verifying/downloading geNomad DB"
    shell:
        r"""
        if [ ! -f "{output}" ]; then
            mkdir -p "{GENOMAD_DB}"
            source /root/miniconda3/etc/profile.d/conda.sh
            conda activate {GENOMAD_ENV}
            genomad download-database "{GENOMAD_DB}"
            [ -f "{output}" ] || echo "v1.x (mark)" > "{output}"
        else
            echo "✔ geNomad DB already exist"
        fi
        """

# ===================== geNomad =====================
rule genomad:
    input:
        file = f"{ANALYSIS}" + "/{sample}/genome_analysis/prokka_{sample}/{sample}.fna",
        DB   = GENOMAD_DB
    output:
        directory(f"{ANALYSIS}" + "/{sample}/genome_analysis/geNomad_{sample}")
    log:
        f"{ANALYSIS}" + "/{sample}/log/geNomad/{sample}.log"
    threads: 8
    message: "geNomad: {wildcards.sample}"
    shell:
        r"""
        source /root/miniconda3/etc/profile.d/conda.sh
        conda activate {GENOMAD_ENV}
        genomad end-to-end --cleanup --splits {threads} {input.file} {output} {input.DB} > {log} 2>&1
        """

# ===================== Abricate DBs (Mark) =====================
rule abricate_db_ready:
    output:
        ABRICATE_READY
    message: "Preparing Abricate DBs"
    shell:
        r"""
        if [ ! -f "{output}" ]; then
            source /root/miniconda3/etc/profile.d/conda.sh
            conda activate {ABRICATE_ENV}
            abricate --setupdb
            touch "{output}"
        else
            echo "✔ Abricate DBs ready"
        fi
        """

# ===================== Abricate (All DBs instaled) =====================
rule abricate:
    input:
        fna   = f"{ANALYSIS}" + "/{sample}/genome_analysis/prokka_{sample}/{sample}.fna",
        ready = ABRICATE_READY
    output:
        summary = f"{ANALYSIS}" + "/{sample}/genome_analysis/abricate_{sample}/{sample}_summary.tab"
    message:
        "Abricate (All DBs): {wildcards.sample}"
    shell:
        r"""
        source /root/miniconda3/etc/profile.d/conda.sh
        conda activate {ABRICATE_ENV}

        outdir="{ANALYSIS}/{wildcards.sample}/genome_analysis/abricate_{wildcards.sample}"
        mkdir -p "$outdir"

        DBS="card resfinder argannot megares plasmidfinder ecoli_vf ncbi vfdb"
        
        for db in $DBS; do
            abricate --db "$db" {input.fna} > "$outdir/{wildcards.sample}_${{db}}.tab"
        done

        # Create the summary only if there are files .tab
        if compgen -G "$outdir/*.tab" > /dev/null; then
            abricate --summary "$outdir"/*.tab > "{output.summary}"
        else
            : > "{output.summary}"
        fi
        """

# ===================== PADLOC DB (Mark) =====================
rule padloc_db:
    output:
        PADLOC_DB_MARK
    message:
        "PADLOC-DB Verifying/updating"
    shell:
        r"""
        if [ ! -f "{output}" ]; then
            source /root/miniconda3/etc/profile.d/conda.sh
            conda activate {PADLOC_ENV}
            padloc --db-update || true
            mkdir -p "$(dirname "{output}")"
            echo "PADLOC DB updated" > "{output}"
        else
            echo "✔ PADLOC-DB is already updated (local mark)."
        fi
        """

# ===================== PADLOC =====================
rule padloc:
    input:
        fna    = f"{ANALYSIS}" + "/{sample}/genome_analysis/prokka_{sample}/{sample}.fna",
        dbmark = PADLOC_DB_MARK   # <— aquí
    output:
        directory(f"{ANALYSIS}" + "/{sample}/genome_analysis/padloc_{sample}")
    log:
        f"{ANALYSIS}" + "/{sample}/log/PADLOC/{sample}.log"
    threads: 8
    message: "PADLOC: {wildcards.sample}"
    shell:
        r"""
        source /root/miniconda3/etc/profile.d/conda.sh
        conda activate {PADLOC_ENV}
        mkdir -p {output}
        padloc --fna {input.fna} --outdir {output} --cpu {threads} > {log} 2>&1
        """

# ===================== CRISPRidentify =====================
CRISPRIDENTIFY_PY = "/root/proyectos/analisis_genomico/CRISPRidentify/CRISPRidentify.py"

rule crispridentify:
    input:
        fna = f"{ANALYSIS}" + "/{sample}/genome_analysis/prokka_{sample}/{sample}.fna"
    output:
        directory(f"{ANALYSIS}" + "/{sample}/genome_analysis/crispridentify_{sample}")
    log:
        f"{ANALYSIS}" + "/{sample}/log/CRISPRidentify/{sample}.log"
    message: "CRISPRidentify: {wildcards.sample}"
    shell:
        r"""
        source /root/miniconda3/etc/profile.d/conda.sh
        conda activate {CRISPRIDENTIFY_ENV}
        mkdir -p {output}
        python "{CRISPRIDENTIFY_PY}" \
            --file {input.fna} \
            --result_folder {output} \
            --cas True \
            --fasta_report True > {log} 2>&1
        """
# ===================== Restfinder_Pointfinder-finder =====================
rule resfinder_pointfinder_db:
    output:
        RESFINDER_DB_MARK
    message: "Verifying DBs ResFinder/PointFinder (marker)"
    shell:
        r"""
        set -euo pipefail
        mkdir -p "{RESFINDER_DB_DIR}" "{POINTFINDER_DB_DIR}"
        [ -f "{RESFINDER_DB_DIR}/INSTALL.py" ] || (echo "Missing resfinder_db" && exit 1)
        [ -f "{POINTFINDER_DB_DIR}/INSTALL.py" ] || (echo "Missing pointfinder_db" && exit 1)
        echo "OK" > "{output}"
        """

rule resfinder_pointfinder:
    input:
        fna   = f"{ANALYSIS}" + "/{sample}/genome_analysis/prokka_{sample}/{sample}.fna",
        db_ok = RESFINDER_DB_MARK
    output:
        directory(f"{ANALYSIS}" + "/{sample}/genome_analysis/resfinder_pointfinder_{sample}")
    log:
        f"{ANALYSIS}" + "/{sample}/log/resfinder_pointfinder/{sample}.log"
    threads: 8
    message: "ResFinder+PointFinder: {wildcards.sample}"
    shell:
        r"""
        set -euo pipefail
        source /root/miniconda3/etc/profile.d/conda.sh
        conda activate {RESFINDER_ENV}

        mkdir -p {output}

        # Species: ResFinder accepts scientific names; for E. coli, it is "Escherichia coli".
        # -l y -t are coverage and identity (0-1). Adjust according to your judgment.
        python -m resfinder \
            -ifa {input.fna} \
            -o {output} \
            -s "Escherichia coli" \
            --acquired --point \
            -l 0.60 -t 0.80 \
            -db_res "{RESFINDER_DB_DIR}" \
            -db_point "{POINTFINDER_DB_DIR}" \
            > {log} 2>&1
        """


# ===================== defensefinder =====================
rule defensefinder_db:
    output:
        DEFENSEF_DB_MARK
    message:
        "📥 Defense-finder: update models + set CasFinder 3.1.0"
    shell:
        r"""
        set -euo pipefail
        source /root/miniconda3/etc/profile.d/conda.sh
        conda activate {DEFENSEF_ENV}

        mkdir -p "$(dirname "{output}")"

        # 1) Update models (this installs CasFinder 3.1.1 by defect)
        defense-finder update -f --force-reinstall

        # 2) Replazed CasFinder by 3.1.0 (avoids errors from versions 2.0 vs 2.1)
        rm -rf /root/.macsyfinder/models/CasFinder
        mkdir -p /tmp/casfinder_fix && cd /tmp/casfinder_fix
        wget -q -O CasFinder-3.1.0.tar.gz https://github.com/macsy-models/CasFinder/archive/refs/tags/3.1.0.tar.gz
        tar -xzf CasFinder-3.1.0.tar.gz
        mkdir -p /root/.macsyfinder/models/CasFinder
        rsync -a CasFinder-3.1.0/ /root/.macsyfinder/models/CasFinder/
        rm -rf /tmp/casfinder_fix

        echo "DefenseFinder models ok; CasFinder pinned to 3.1.0" > "{output}"
        """


rule defense_finder:
    input:
        faa    = f"{ANALYSIS}" + "/{sample}/genome_analysis/prokka_{sample}/{sample}.faa",
        dbmark = DEFENSEF_DB_MARK
    output:
        directory(f"{ANALYSIS}" + "/{sample}/genome_analysis/defensefinder_{sample}")
    log:
        f"{ANALYSIS}" + "/{sample}/log/defensefinder/{sample}.log"
    threads: 8
    message: "Defense-finder: {wildcards.sample}"
    shell:
        r"""
        source /root/miniconda3/etc/profile.d/conda.sh
        conda activate {DEFENSEF_ENV}

        mkdir -p {output}

        defense-finder run {input.faa} -o {output} > {log} 2>&1
        """
# ===================== CCTyper =====================
rule cctyper:
    input:
        fna = f"{ANALYSIS}" + "/{sample}/genome_analysis/prokka_{sample}/{sample}.fna"
    output:
        outdir = directory(f"{ANALYSIS}" + "/{sample}/genome_analysis/cctyper_{sample}")
    log:
        f"{ANALYSIS}" + "/{sample}/log/cctyper/{sample}.log"
    threads: 8
    message: "CCTyper: {wildcards.sample}"
    shell:
        r"""
        set -euo pipefail
        source /root/miniconda3/etc/profile.d/conda.sh
        conda activate {CCTYPER_ENV}

        # Output father
        mkdir -p "$(dirname "{output.outdir}")"

        # Create a temporary path THAT DOES NOT EXIST (cctyper will create it)
        tmpbase="$(mktemp -u)"
        tmpdir="${{tmpbase}}_cctyper"

        # Cleaning just in case
        rm -rf "$tmpdir"

        # Execute: the output must NOT exist
        cctyper -t {threads} {input.fna} "$tmpdir" > {log} 2>&1

        # Minimal validation
        if [ ! -d "$tmpdir" ] || [ -z "$(ls -A "$tmpdir" 2>/dev/null)" ]; then
            echo "ERROR: CCTyper did not generate files (tmpdir empty or not created)." >> {log}
            exit 1
        fi

        # Replace final output
        rm -rf "{output.outdir}"
        mv "$tmpdir" "{output.outdir}"
        """

# ===================== ECTyper =====================        
rule ectyper_mash_db:
    output:
        ECTYPER_MASH_MARK
    shell:
        r"""
        set -euo pipefail
        [ -s "{ECTYPER_MASH_SKETCH}" ] || (echo "Missing ECTyper MASH sketch: {ECTYPER_MASH_SKETCH}" && exit 1)
        echo "OK" > "{output}"
        """

rule ECTyper:
    input:
        fna    = f"{ANALYSIS}" + "/{sample}/genome_analysis/prokka_{sample}/{sample}.fna",
        mashok = ECTYPER_MASH_MARK
    output:
        directory(f"{ANALYSIS}" + "/{sample}/genome_analysis/ECTyper_{sample}")
    log:
        f"{ANALYSIS}" + "/{sample}/log/ECTyper/{sample}.log"
    threads: 8
    message: "ECTyper: {wildcards.sample}"
    shell:
        r"""
        set -euo pipefail
        source /root/miniconda3/etc/profile.d/conda.sh
        conda activate {ECTYPER_ENV}

        # Avoid blocks if a stale lock remains:
        rm -f /root/proyectos/analisis_genomico/Bac_pipeline_AC/input_from_prokka/envs/ectyper_env/lib/python*/site-packages/ectyper/Data/.lock || true

        mkdir -p {output}
        ectyper -i {input.fna} -o {output} -c {threads} --verify --reference {ECTYPER_MASH_SKETCH} --pathotype > {log} 2>&1
        """

# ===================== mlst =====================    
rule mlst:
    input:
        fna = f"{ANALYSIS}" + "/{sample}/genome_analysis/prokka_{sample}/{sample}.fna"
    output:
        outdir = directory(f"{ANALYSIS}" + "/{sample}/genome_analysis/mlst_{sample}"),
        tsv    = f"{ANALYSIS}" + "/{sample}/genome_analysis/mlst_{sample}/{sample}_mlst.tsv"
    log:
        f"{ANALYSIS}" + "/{sample}/log/mlst/{sample}.log"
    threads: 1
    message: "MLST: {wildcards.sample}"
    shell:
        r"""
        set -euo pipefail
        source /root/miniconda3/etc/profile.d/conda.sh
        conda activate {MLST_ENV}

        mkdir -p {output.outdir}

        # mlst prints one TSV line per genome; we save the output to a file
        mlst --scheme ecoli_achtman_4 {input.fna} > {output.tsv} 2> {log}
        """




