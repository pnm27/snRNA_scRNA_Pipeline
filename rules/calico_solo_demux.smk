from lib.aligner import count_matrix as aligner_count_matrix



rule run_calico_solo:
    input:
        bustools_h5ad=f"{config['hashsolo_demux_pipeline']['h5ad_bustools_dir']}{config['fold_struct_demux']}{config['hashsolo_demux_pipeline']['bustools_h5ad']}",
        starsolo_out=aligner_count_matrix(config)

    output:
        f"{config['hashsolo_demux_pipeline']['calico_solo_dir']}{config['fold_struct_demux']}{config['hashsolo_demux_pipeline']['calico_solo_h5ad']}"
  
    params:
        mito=config['max_mito_percentage'],  # Max mitochodrial genes content per cell
        min_genes=config['min_genes_per_cell'], # Min #genes per cell
        min_cells=config['min_cells_per_gene'],  # Min #cells expressing a gene for it to pass the filter
        genes_info=config['STARsolo_pipeline']['genome_pick']['gene_info_file'], # File containing gene names and gene ids for annotations
        mito_prefix=config['mito_prefix'], # Mitochondrial genes' (names') prefix
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
