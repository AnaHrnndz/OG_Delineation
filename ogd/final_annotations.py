from collections import defaultdict
import ogd.utils as utils



##  8.Annotate root ##
def annotate_root(ogs_info, t, name_tree, total_mems_in_tree, sp_set, seqs_in_ogs, recover_seqs, taxonomy_db, args):

    """
    Add properties in root
    Add info about the root in ogs_info dict
    Need for web
    """

    # Save parameters as property

    parametrs2save = ["lineage_thr@"+str(args.lineage_thr), "best_tax_thr@"+str(args.best_tax_thr), 
                    "sp_loss_perc@"+str(args.sp_loss_perc), "inherit_out@"+str(args.inherit_out), "sp_ovlap_all@"+str(args.so_all)]
    
    if args.so_bact != None:
        parametrs2save.append("sp_ovlap_bact@"+str(args.so_bact))
    if args.so_euk != None:
        parametrs2save.append("sp_ovlap_euk@"+str(args.so_euk))
    if args.so_arq != None:
        parametrs2save.append("sp_ovlap_arq@"+str(args.so_arq))
    
    parameters_str = '|'.join(parametrs2save)

    t.add_prop('parameters', parameters_str)

    # Save general results as properties, needed for web

    seqs_out_ogs = total_mems_in_tree.difference(seqs_in_ogs)
    
    results2save = ["Tree_name@"+name_tree, "Total_seqs@"+str(len(total_mems_in_tree)), "Total_species@"+str(len(sp_set)),
                    "Seqs_in_OGs@"+str(len(seqs_in_ogs)), "Recovery_seqs@"+str(len(recover_seqs)), "Seqs_out_OGs@"+str(len(seqs_out_ogs)), 
                    "Num_OGs@"+str(len(ogs_info))]
    
    general_results_str = '|'.join(results2save)
    t.add_prop('general_result', general_results_str)


    # Save taxlevel OGs info, needed for web

    lca_all_dups = set()
    taxlev2ogs = defaultdict(set)
    taxlev2mems = defaultdict(set)

    for dup_node in t.search_nodes(node_create_og='True'):
        lca_dup = dup_node.props.get('lca_node')
        
        mems_dup = set(dup_node.props.get('leaves_in'))
        lca_all_dups.add(lca_dup)
        taxlev2ogs[lca_dup].add(dup_node.name)
        taxlev2mems[lca_dup].update(mems_dup)

   
    taxlev_list = []
    #for taxlev, og in taxlev2ogs.items():
        
        #if isinstance(taxlev, str):
        #taxlev = int(taxlev)

        # if (str(taxonomyogd_env_db).split('.')[1]) == 'ncbi_taxonomy':
            # sci_name = taxonomy_db.get_taxid_translator([taxlev])[taxlev]
        
        # elif (str(taxonomy_db).split('.')[1]) == 'gtdb_taxonomy':
            # sci_name = taxlev
          
        # sci_name_taxid = sci_name+'_'+str(taxlev)
        # ogs_str = '_'.join(list(og))
        # num_mems = len(taxlev2mems[taxlev])

        # tax_str = '|'.join([sci_name_taxid, ogs_str,str(num_mems)])

        # taxlev_list.append(tax_str)

    # taxlev_str = '@'.join(taxlev_list)
    # t.add_prop('taxlev2ogs', taxlev_str)

    t.add_prop("OGD_annot", True)


##  9. Flag seqs out OGs ##
def flag_seqs_out_og(t, seqs_in_ogs, total_mems_in_tree):

    """
        Flags seqs that do no belong to any OG
        Could be for taxonomic outlier or for branch lenght
    """
    seqs_out_og = total_mems_in_tree.difference(seqs_in_ogs)
    for leaf in t:
        if leaf.name in seqs_out_og:
            leaf.add_prop('seq_out_og', "true")

    t, props = utils.run_clean_properties(t)

    return t

