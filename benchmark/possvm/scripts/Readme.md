# Possvm benchmark

## Default parameters
1. Transform leaves of original trees to ncbi taxid format  
    `python benchmark/possvm/scripts/convert_code2taxid.py`  

2. Annote possvm trees with:  
    - emmaper results  
    - Pfam fams  
    - OGs at euk and metazoa level  
    - Group assigned with possvm  
    - Refog from homeoDB  

    `python benchmark/possvm/scripts/annotate_possvm_trees.py`  

3. Run OG Delination command line  
    -   `python og_delineation.py --tree benchmark/possvm/annotated_trees/annotated_transform_ANTP.genes.iqtree.treefile --raw_alg benchmark/possvm/transform_alg/tranform_ANTP.genes.l.fasta  --user_taxonomy taxonomy/e6.taxa.sqlite --output_path benchmark/possvm/results_ogs/antp &>benchmark/possvm/results_ogs/antp.log`  

    -   `python og_delineation.py --tree benchmark/possvm/annotated_trees/annotated_transform_PRD.genes.iqtree.treefile --raw_alg benchmark/possvm/transform_alg/tranform_PRD.genes.l.fasta --user_taxonomy taxonomy/e6.taxa.sqlite --output_path benchmark/possvm/results_ogs/prd/ &>benchmark/possvm/results_ogs/prd.log`  
    
    -   `python og_delineation.py --tree benchmark/possvm/annotated_trees/annotated_transform_TALE.genes.iqtree.treefile --raw_alg benchmark/possvm/transform_alg/tranform_TALE.genes.l.fasta --user_taxonomy taxonomy/e6.taxa.sqlite --output_path benchmark/possvm/results_ogs/tale/ &>benchmark/possvm/results_ogs/tale.log`  


4. Get stats for benchmark  
    -   python stats_posvvm.py /data/projects/og_delineation/benchmark/possvm/original_data/reference.antp_all_ite_mcl.csv /data/projects/og_delineation/benchmark/possvm/results_ogs/antp/recovery_ogs_info.tsv /data/projects/og_delineation/benchmark/possvm/original_data/ANTP.possom.ortholog_groups.csv /data/projects/og_delineation/benchmark/possvm/results_ogs/antp

    -   python stats_posvvm.py /data/projects/og_delineation/benchmark/possvm/original_data/reference.tale_all_ite_mcl.csv /data/projects/og_delineation/benchmark/possvm/results_ogs/tale/recovery_ogs_info.tsv /data/projects/og_delineation/benchmark/possvm/original_data/TALE.possom.ortholog_groups.csv /data/projects/og_delineation/benchmark/possvm/results_ogs/tale

    -   python stats_posvvm.py /data/projects/og_delineation/benchmark/possvm/original_data/reference.prds_all_ite_mcl.csv /data/projects/og_delineation/benchmark/possvm/results_ogs/prd/recovery_ogs_info.tsv /data/projects/og_delineation/benchmark/possvm/original_data/PRD.possom.ortholog_groups.csv /data/projects/og_delineation/benchmark/possvm/results_ogs/prd


*********************************************
## Species overlap for Eukaryota 0.0


