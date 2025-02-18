# Resource Allocation ------------------
def allocate_mem_RCS(wildcards, attempt):
    return 4500*attempt+2500

def allocate_time_RCS(wildcards, attempt):
    return 5*attempt+5

# --------------------------------------


rule run_calico_solo:
    input:
        f"{config['hashsolo_demux_pipeline']['h5ad_bustools_dir']}{config['fold_struct_demux']}{config['hashsolo_demux_pipeline']['bustools_h5ad']}",
        starsolo_out=f"{config['STARsolo_pipeline']['bams_dir']}{config['fold_struct']}{config['STARsolo_pipeline']['genefull_matrix']}" #get_STARsolo_mat

    # priority: 8
   
    output:
        f"{config['hashsolo_demux_pipeline']['calico_solo_dir']}{config['fold_struct_demux']}{config['hashsolo_demux_pipeline']['calico_solo_h5ad']}"
  
    params:
        mito=config['max_mito_percentage'],  # Max mitochodrial genes content per cell
        min_genes=config['min_genes_per_cell'], # Min #genes per cell
        min_cells=config['min_cells_per_gene'],  # Min #cells expressing a gene for it to pass the filter
        genes_info=config['gene_info_file'], # File containing gene names and gene ids for annotations
        mito_prefix=config['mito_prefix'] # Mitochondrial genes' (names') prefix

    # For snakemake < v8
    # threads: 2

    conda: "../envs/hashsolo.yaml"

    resources:
        cpus_per_task=2, # For snakemake > v8
        mem_mb=allocate_mem_RCS,
        time_min=allocate_time_RCS

    shell: 
        """
        python3 helper_py_scripts/create_h5ad_from_calico_solo.py {input[0]} {input[1]} {output} \
            {params.genes_info} -m {params.mito} -g {params.min_genes} \
            -c {params.min_cells} --mito_prefix {params.mito_prefix}
        sleep 100
        """
