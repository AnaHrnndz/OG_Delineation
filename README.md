# Orthologous Group Delineation

Orthologous Groups Delineation (OGD) algorithm attempt to identify orthologous groups in trees based on duplication events.

MAIN FUNCTION, run all steps of the analysis:

        1.  Load files needed for the analysis.  
                Tree, taxonomy, reference tree (species tree), taxonomy counter  
        2.  Tree setup (Pre-analysis). Run some basic analysis:   
                - Resolve polytomies  
                - Rooting;Midpoint or MinVar  
                - NCBI annotation  
                - Save original species and original leaf names  
        3. Outliers and Dups score:  
                - Detect long branches  
                - Detect taxonomical outliers  
                - Caldulate: Species overlap, Score1, Score2, Inpalalogs rate, Duplication score and Linage lost  
        4. Detect HQ-Duplications:  
                Select high quality duplication that create Basal Orthologs Groups (OGs)  
        5. Get OGs
        6. Optionally skip get all orthologs pairs
        7. Optionally modify OGs by recovering sequences (see RECOVERY PIPELINE)
        8. Optionally add annotations from emapper  (see EMAPPER ANNOTATION)
                8.1 Run emapper and annotate with treeprofiler
                8.2 Only annotate with treeprofiler
        9. Annotate root 
        10. Flag seqs out OGs
        11. Write output files


RECOVERY PIPELINE:

    If user provide a alignment file AND some sequences (members) are not included in at least one Core-OGs
            1. Write a fasta file with all raw seqs left out of the Core-OGs and for each Core-OGs write a fasta file
                1.1 Regular mode: Re-aling sequences
                1.2 Fast mode: use aligned sequnces
            2. Create HMM file for each Core-OGs fasta file and Build HMM DB
            3. Run Hmmscan to assign seqs out Core-OGs to Core-OGs
            4. For each seq out Core-OGs, select best Core-OGs
            5. Expand_hmm: Create a dict with all the OGs that have recover some seq
            6. update_og_info: Update og_info with all the seqs taht have been recovered
            7. update_tree: Update tree with all the seqs taht have been recovered
            8. write_ogs_info for recovery
            9. write_best_match
            10. Tar the folder with all the fasta file, HMMs, etc

EMAPPER ANNOTATION:

    Annotate tree sequences with eggnog-mapper, user needs to provide alignment file