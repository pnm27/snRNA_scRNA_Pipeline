from lib.aligner import count_matrix as aligner_count_matrix



rule run_calico_solo:
    input:
        bustools_h5ad=f"{config['modules']['hashsolo']['h5ad_bustools_dir']}{config['modules']['hashsolo']['fold_struct']}{config['modules']['hashsolo']['bustools_h5ad']}",
        starsolo_out=aligner_count_matrix(config)

    output:
        f"{config['modules']['hashsolo']['calico_solo_dir']}{config['modules']['hashsolo']['fold_struct']}{config['modules']['hashsolo']['calico_solo_h5ad']}"
  
    params:
        mito=config['experiment']['qc']['max_mito_percentage'],  # Max mitochodrial genes content per cell
        min_genes=config['experiment']['qc']['min_genes_per_cell'], # Min #genes per cell
        min_cells=config['experiment']['qc']['min_cells_per_gene'],  # Min #cells expressing a gene for it to pass the filter
        genes_info=config['modules']['STARsolo']['genome_pick']['gene_info_file'], # File containing gene names and gene ids for annotations
        mito_prefix=config['references']['mito_prefix'], # Mitochondrial genes' (names') prefix
        script="helper_py_scripts/create_h5ad_from_calico_solo.py"

    conda: "../envs/hashsolo.yaml"

    resources:
        **get_resources("run_calico_solo")

    shell: 
        """
        python3 {params.script} {input.bustools_h5ad} \
            {input.starsolo_out} {output} {params.genes_info} \
            -m {params.mito} -g {params.min_genes} \
            -c {params.min_cells} --mito_prefix {params.mito_prefix}
        sleep 100
        """
