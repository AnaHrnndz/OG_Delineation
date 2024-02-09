from collections import defaultdict
import utils


def get_messy_groups(t, taxonomy_db):

    messy_ogs = defaultdict(dict)
    c=0

    for n in t.traverse():

        lca_target = n.props.get('lca_node')

        if n.is_root:
            if n.children[0].props.get('node_is_og') == 'True' and n.children[1].props.get('node_is_og') == 'True':
                pass
            elif n.props.get('node_is_og') == 'True':
                pass

            else:
                for ch in n.children:
                    all_mems = ch.props.get('_leaves_in', set())
                    mem2remove = set()
                    nodes2ckech = set(ch.search_nodes(node_is_og='True', lca_dup=lca_target))

                    ogs_down = set()
                    for n_ in nodes2ckech:
                        ogs_down.add(n_.name)

                    nodes2ckech.update(ch.search_nodes(evoltype_2='D', lca_node=lca_target))
                    for n_ in nodes2ckech:
                        mem2remove.update(n_.props.get('_leaves_in'))

                    diff = all_mems.difference(mem2remove)
                    if  len(diff)> 0:
                        name = 'mOG_'+str(c)
                        c+=1

                        #messy_ogs[str(lca_target)][name] = diff
                        messy_ogs[name] = defaultdict()
                        messy_ogs[name]['TaxoLevel'] = str(lca_target)
                        messy_ogs[name]['SciName_TaxoLevel'] = taxonomy_db.get_taxid_translator([lca_target])[lca_target]
                        messy_ogs[name]['AssocNode'] = ch.name
                        messy_ogs[name]['OG_down'] = '|'.join(list(ogs_down))
                        messy_ogs[name]['Mems'] = set()
                        sp_set = set()
                        for lname in diff:
                            if 'taxo_outlier' not in ch[lname].props:
                                old_mOG = ch[lname].props.get('mOG', str())
                                new_mOG = old_mOG+'@'+name
                                ch[lname].add_prop('mOG', new_mOG)
                                messy_ogs[name]['Mems'].add(lname)
                                sp_set.add(ch[lname].props.get('taxid'))
                        messy_ogs[name]['NumSP'] = len(sp_set)
                        messy_ogs[name]['NumMems'] = len(messy_ogs[name]['Mems'])


        else:
            if n.props.get('evoltype_2')=='D':

                ogs_up, dups_up = utils.check_nodes_up(n)
                if n.children[0].props.get('node_is_og') == 'True' and n.children[1].props.get('node_is_og') == 'True':
                    pass
                elif len(ogs_up)>0:
                    pass

                else:

                    for ch in n.children:
                        all_mems = ch.props.get('_leaves_in', set())
                        mem2remove = set()
                        nodes2ckech = set(ch.search_nodes(node_is_og='True', lca_dup=lca_target))
                        nodes2ckech.update(ch.search_nodes(evoltype_2='D', lca_node=lca_target))
                        for n_ in nodes2ckech:
                            mem2remove.update(n_.props.get('_leaves_in'))

                        diff = all_mems.difference(mem2remove)

                        if  len(diff)> 1:
                            name = 'mOG_'+str(c)
                            c+=1

                            #messy_ogs[str(lca_target)][name] = diff
                            messy_ogs[name] = defaultdict()
                            messy_ogs[name]['TaxoLevel'] = str(lca_target)
                            messy_ogs[name]['SciName_TaxoLevel'] = taxonomy_db.get_taxid_translator([lca_target])[lca_target]
                            messy_ogs[name]['AssocNode'] = ch.name
                            messy_ogs[name]['OG_down'] = '|'.join(list(ogs_down))
                            messy_ogs[name]['Mems'] = set()
                            sp_set = set()
                            for lname in diff:
                                if 'taxo_outlier' not in ch[lname].props:
                                    old_mOG = ch[lname].props.get('mOG', str())
                                    new_mOG = old_mOG+'@'+name
                                    ch[lname].add_prop('mOG', new_mOG)
                                    messy_ogs[name]['Mems'].add(lname)
                                    sp_set.add(ch[lname].props.get('taxid'))
                            messy_ogs[name]['NumSP'] = len(sp_set)
                            messy_ogs[name]['NumMems'] = len(messy_ogs[name]['Mems'])



    return t,  messy_ogs


# def annotate_messy_ogs(messy_ogs):
    # #OG_name        TaxoLevel       SciName_TaxoLevel       AssocNode       NumSP   OG_down OG_up   NumSeqs
    # for taxid, mog in messy_ogs.items():
        # for mog_name, mog_mems in mog.items()
        # sp_set = set()




# def write_mogs(messy_ogs, path_out):
    # fout_mogs = open(path_out+'/mOGs.tsv', 'w')

    # for taxid, mog in messy_ogs.items():
        # for mog_name, mog_mems in mog.items():
            # fout_mogs.write(mog_name+'\t'+str(taxid)+'\t'+str(len(mog_mems))+'\t'+','.join(list(mog_mems))+'\n')
    # fout_mogs.close()
