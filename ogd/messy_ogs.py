from collections import defaultdict
import utils


def get_messy_groups(t, taxonomy_db):

    messy_ogs = defaultdict(dict)
    seqs_in_messy_ogs = set()
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
                    all_mems = ch.props.get('leaves_in', set())
                    mem2remove = set()
                    nodes2ckech = set(ch.search_nodes(node_is_og='True', lca_dup=lca_target))
                    nodes2ckech.update(ch.search_nodes(evoltype_2='D', lca_node=lca_target))
                    for n_ in nodes2ckech:
                        mem2remove.update(n_.props.get('leaves_in'))

                    diff = all_mems.difference(mem2remove)
                    if  len(diff)>1:

                        c, messy_ogs, seqs_in_messy_ogs = new_mog(c, messy_ogs, ch, taxonomy_db, seqs_in_messy_ogs, lca_target, diff)

        else:
            if n.props.get('evoltype_2')=='D':

                if n.children[0].props.get('node_is_og') == 'True' and n.children[1].props.get('node_is_og') == 'True':
                    pass

                else:

                    for ch in n.children:
                        all_mems = ch.props.get('leaves_in', set())
                        mem2remove = set()
                        nodes2ckech = set(ch.search_nodes(node_is_og='True', lca_dup=lca_target))
                        nodes2ckech.update(ch.search_nodes(evoltype_2='D', lca_node=lca_target))
                        for n_ in nodes2ckech:
                            mem2remove.update(n_.props.get('leaves_in'))

                        diff = all_mems.difference(mem2remove)
                        if  len(diff)>1:

                            c, messy_ogs, seqs_in_messy_ogs = new_mog(c, messy_ogs, ch, taxonomy_db, seqs_in_messy_ogs, lca_target, diff)



            # elif n.props.get('node_is_og') == 'True':

                # node_lin = n.props.get('lineage')

                # # Get all lineages in OG
                # lin2mems = defaultdict(set)

                # all_mems = n.props.get('_leaves_in', set())

                # for mem in all_mems:

                    # for lin in n[mem].props.get('lineage').split('|'):

                        # if int(lin) not in node_lin:
                            # lin2mems[int(lin)].add(mem)


                # for lin, mems in lin2mems.items():

                    # nodes2check = set(ch.search_nodes(node_create_og='True',lca_node=lin))
                    # mems2remove = set()
                    # for n_ in nodes2check:
                        # mems2remove.update(n_.props.get('_leaves_in'))

                    # diff = mems.difference(mems2remove)
                    # if  len(diff)>1:
                        # c, messy_ogs, seqs_in_messy_ogs = new_mog(c, messy_ogs, n, taxonomy_db, seqs_in_messy_ogs, lin, diff)



    messy_ogs = add_mogs_up_down(messy_ogs, t)

    return t,  messy_ogs, seqs_in_messy_ogs



def new_mog(c, messy_ogs, ch, taxonomy_db, seqs_in_messy_ogs, lca_target, diff):

    name = 'mOG_'+str(c)
    c+=1

    messy_ogs[name] = defaultdict()
    messy_ogs[name]['TaxoLevel'] = str(lca_target)
    if lca_target == 'r_root':
        messy_ogs[name]['SciName_TaxoLevel'] = 'r_root'
    else:
        if (str(taxonomy_db).split('.')[1]) == 'ncbi_taxonomy':
            messy_ogs[name]['SciName_TaxoLevel'] = taxonomy_db.get_taxid_translator([lca_target])[lca_target]
        elif (str(taxonomy_db).split('.')[1]) == 'gtdb_taxonomy':
            messy_ogs[name]['SciName_TaxoLevel'] = lca_target
        
    messy_ogs[name]['AssocNode'] = ch.name
    messy_ogs[name]['OG_down'] = set()
    messy_ogs[name]['OG_up'] = set()
    messy_ogs[name]['Mems'] = set()

    for lname in diff:

        if lname not in ch.props.get('leaves_out'):
            old_mOG = ch[lname].props.get('mOG', str())
            new_mOG = old_mOG+'@'+name
            ch[lname].add_prop('mOG', new_mOG)
            messy_ogs[name]['Mems'].add(lname)
            seqs_in_messy_ogs.add(lname)
            ogs_up, dups_up = utils.check_nodes_up(ch[lname])
            messy_ogs[name]['OG_down'].update(set(ogs_up))

    messy_ogs[name]['NumMems'] = len(messy_ogs[name]['Mems'])

    return c ,messy_ogs, seqs_in_messy_ogs


def add_mogs_up_down(messy_ogs, t):
    for mog_name, mog_info in messy_ogs.items():
        all_mogs = set()
        mog_level = int(mog_name.split('_')[1])

        for s in mog_info['Mems']:
            all_mogs.update(set(t[s].props.get('mOG').split('@')))


        for mog in all_mogs:
            if mog != '':
                lev = int(mog.split('_')[1])
                if lev == mog_level:
                    pass
                elif lev > mog_level:
                    messy_ogs[mog_name]['OG_down'].add(mog)
                elif lev < mog_level:
                    messy_ogs[mog_name]['OG_up'].add(mog)

        if len(messy_ogs[mog_name]['OG_down']) == 0:
            messy_ogs[mog_name]['OG_down'].add('-')

        if len(messy_ogs[mog_name]['OG_up']) == 0:
            messy_ogs[mog_name]['OG_up'].add('-')

    return messy_ogs


def annotate_messy_og(t, messy_ogs):


    for mog_name, info in messy_ogs.items():
        sp_set = set()
        recseq_down = set()
        ogs_up = messy_ogs[mog_name]['OG_down']

        for m in messy_ogs[mog_name]['Mems']:
            sp_set.add(t[m].props.get('taxid'))

        for og_ in  ogs_up:
            if og_ != '-':
                if not og_.startswith('mOG_') :
                    if t[og_].props.get('recovery_seqs')!= None:
                        for s in t[og_].props.get('recovery_seqs'):
                            recseq_down.add(s)
                            sp_set.add(t[s].props.get('taxid'))


        if recseq_down != None:
            messy_ogs[mog_name]['RecoverySeqs'] = recseq_down
            messy_ogs[mog_name]['NumRecoverySeqs'] = len(recseq_down)
        else:
            messy_ogs[mog_name]['RecoverySeqs'] = ['-']
            messy_ogs[mog_name]['NumRecoverySeqs'] = str(0)


        messy_ogs[mog_name]['NumSP'] = len(sp_set)


    return messy_ogs

