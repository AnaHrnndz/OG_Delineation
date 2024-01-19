import os
import shutil
import glob
import subprocess
import tarfile
import re
from collections import defaultdict
from ete4 import SeqGroup



####    FUNCTIONS FOR RECOVERING SEQUENCES PIPELINE ####

def recover_sequences(tree, t, alg, total_mems_in_tree, total_mems_in_ogs, name_tree, path_out, ogs_info, mode):

    name_fam = os.path.basename(tree).split('.')[0]
    pathout = path_out+'/'+name_fam+'_aln_hmm'
    if not os.path.exists(pathout):
        os.mkdir(pathout)
    fasta = alg

    # 1. Write a fasta file
    run_write_fastas(t, fasta, name_tree, pathout, ogs_info, total_mems_in_ogs, total_mems_in_tree, mode)

    if mode == "regular":  # regular mode: Re-aling sequences
        run_alignments(pathout)
    else:
        pass  # we just go with what we had

    # 2. Create HMM file for each Core-OGs fasta file and Build HMM DB
    run_create_hmm_og(pathout)

    # 3. Run Hmmscan to assign seqs out Core-OGs to Core-OGs
    tblfile = run_hmmscan(pathout)

    # 4. For each seq out Core-OGs, select best Core-OGs
    best_match = get_best_match(tblfile)
    recover_seqs = set(best_match.keys())

    t = expand_hmm(t, best_match)

    total_mems_in_ogs.update(recover_seqs)


    # 9. write_best_match
    write_best_match(best_match, path_out, name_tree)

    output_filename = path_out+'/'+name_fam+'.tar.gz'

    # 10. Tar the folder with all the fasta file, HMMs, etc
    make_tarfile(output_filename, pathout)
    shutil.rmtree(pathout)

    return recover_seqs, best_match


def run_write_fastas(t, fasta, name_tree, path_out, ogs_info, total_mems_in_ogs, total_mems_in_tree, mode):

    """
        Function calls  write_outog_seqs() AND write_og_seqs_regular_mode() OR write_og_seqs_fast_mode()
    """

    print('-Fastas')

    fasta = SeqGroup(fasta)

    diff = (total_mems_in_tree.difference(total_mems_in_ogs))

    name_fam = name_tree.split('.')[0]
    not_og_fasta = path_out+'/'+name_fam+'.notOG.faa'

    write_outog_seqs(fasta, diff, not_og_fasta)

    if mode == "regular":
        write_og_seqs_regular_mode(t, fasta, ogs_info, path_out)
    elif mode == "fast":
        write_og_seqs_fast_mode(t, fasta, ogs_info, path_out)

    return

def write_outog_seqs(fasta, diff, not_og_fasta):

    """
        Write fasta file with seqs that do not belong to any OG
    """

    with open(not_og_fasta, 'w') as  f_out:
        for name_seq in diff:
            aa = fasta.get_seq(name_seq)
            clean_aa = aa.replace('-','')
            f_out.write('>'+name_seq+'\n'+clean_aa+'\n')

def write_og_seqs_regular_mode(t, fasta, ogs_info, path_out):

    """
        Write one fasta file per each OG
        In regular mode, sequences will be realing, so gaps are removed
    """

    set_trees  = set(t.search_nodes(node_is_og='True'))


    for subtree in set_trees:
        if 'is_root' not in subtree.props:
            lca = str(subtree.props.get('lca_node'))
            node_name = subtree.name
            with open(path_out+'/'+node_name+'_'+lca+'.raw_fasta.faa', 'w') as f_out:
                list_mems = list(subtree.props.get('_leaves_in'))
                if len(list_mems) >0:
                    for m in list_mems:
                        aa = fasta.get_seq(m)
                        clean_aa = aa.replace('-','')
                        f_out.write('>'+m+'\n'+clean_aa+'\n')

def write_og_seqs_fast_mode(t, fasta, ogs_info, path_out):

    """
        Write one fasta file per each OG
        In fast mode, sequences wont be realing, so gaps are keept
    """

    set_trees  = set(t.search_nodes(node_is_og='True'))

    for subtree in set_trees:
        if 'is_root' not in subtree.props:
            lca = str(subtree.props.get('lca_node'))
            node_name = subtree.name
            with open(path_out+'/'+node_name+'_'+lca+'.aln', 'w') as f_out:
                list_mems = list(subtree.props.get('_leaves_in'))
                if len(list_mems) >0:
                    for m in list_mems:
                        aa = fasta.get_seq(m)
                        f_out.write('>'+m+'\n'+aa+'\n')

def run_alignments(path_out):

    """
        Function call run_mafft() or run_famsa() depend of the size of OG
    """

    total_aln = glob.glob(path_out+'/*.raw_fasta.faa')
    for aln in total_aln:
        aln_seq = SeqGroup(aln)
        name_og = os.path.basename(aln).split('.')[0]

        if len(aln_seq) <= 100:
            run_mafft(aln,path_out, name_og)
        else:
            run_famsa(aln, path_out, name_og)

def run_mafft(raw_fasta, path_out, name_og):

    """
        Run mafft software for OGs < 100 seqs
    """

    aln_fasta = path_out+'/'+name_og+'.aln'
    subprocess.run(("mafft --auto --anysymbol %s > %s" %(raw_fasta, aln_fasta)), shell = True)

def run_famsa(raw_fasta,path_out, name_og):

    """
        Run famsa software for OGs > 100 seqs
    """

    aln_fasta = path_out+'/'+name_og+'.aln'
    subprocess.run(("famsa %s %s" %(raw_fasta, aln_fasta)), shell = True)

def run_create_hmm_og(path):
    """
        Build HHMM file for each OG-aln fasta
        Create HMM DB with all the HMM files

        TODO: split in 2¿?
    """

    aln_list = glob.glob(path+'/*.aln')
    for aln in aln_list:
        out_hmm = aln+'.hmm'
        subprocess.run(("hmmbuild %s %s" %(out_hmm, aln)), shell = True)

    hmm_list = glob.glob(path+'/*.hmm')
    hmm_db = path+'/hmm_db.hmm'
    for hmm in hmm_list:
        subprocess.run(("cat %s >>%s" %(hmm, hmm_db)), shell = True)


    subprocess.run(("hmmpress %s" %(hmm_db)), shell = True)

def run_hmmscan(path):

    """
        Run Hmmscan
    """

    hmm_db = path+'/hmm_db.hmm'
    outfile = path+'/result_hmmscan.tsv'
    tblfile = path+'/tblout.tsv'
    domtblfile = path+'/domfile.tsv'
    seqs_not_og = glob.glob(path+'/*.notOG.faa')[0]
    print(("hmmscan -o %s %s %s" %(outfile, hmm_db, seqs_not_og)))
    subprocess.run(("hmmscan -o %s --tblout %s --domtblout %s  %s %s" %(outfile, tblfile, domtblfile, hmm_db, seqs_not_og)), shell = True)
    return tblfile

def get_best_match(tblout):

    """
        For each seq,  get best hit from Hmmscan
    """

    best_match = defaultdict(dict)

    with open(tblout) as f:
        for line in f:
            if not line.startswith('#'):
                line = re.sub(' +','\t',line)
                info = line.split('\t')
                name_og = info[0]
                name_seq = info[2]
                score = float(info[5])

                k_ = str()
                score_ = float()
                if name_seq in best_match.keys():
                    for k, val in best_match[name_seq].items():
                        k_ = k
                        score_ = val

                    if score > score_:
                        best_match[name_seq].pop(k_)
                        best_match[name_seq][name_og] = score
                else:
                    best_match[name_seq][name_og] = score


    return(best_match)

def write_best_match(best_match, path, name_tree):

    """
        Write a table with the best match result from hmmscan for each seq
    """
    name_fam = name_tree.split('.',1)[0]
    with open(path+'/'+name_fam+'.best_match.tsv', 'w') as fout:
        for seq, best_og in best_match.items():
            for best_og_name in best_og.keys():
                fout.write('\t'.join([seq, best_og_name])+'\n')

def expand_hmm(t, best_match):

    """
        Create a dict with all the OGs that have recover some seq

        TODO: change name
    """

    for seq_name, best_og in best_match.items():
        for k in best_og.keys():
            best_og_name = k.split('_')[0]

        if 'recovery_seqs' in t[best_og_name].props:
            t[best_og_name].props.get('recovery_seqs').update(seq_name)

        else:
            recovery_set = set()
            recovery_set.add(seq_name)
            t[best_og_name].add_prop('recovery_seqs',  recovery_set)

        if t[best_og_name].props.get('ogs_up') != '-':
            for nup in t[best_og_name].props.get('ogs_up'):
                if 'recovery_seqs' in t[nup].props:
                    t[nup].props.get('recovery_seqs').update(seq_name)

                else:
                    recovery_set = set()
                    recovery_set.add(seq_name)
                    t[nup].add_prop('recovery_seqs',  recovery_set)
    return t

def update_taxlevel2ogs(glob_taxlev2ogs, og_info_recovery, glob_og_info) :

    """
        update taxlevel2og after recovered step
        Needed for web
    """

    for tax, info in glob_taxlev2ogs.items():
        for og in info['ogs_names']:
            if og in og_info_recovery.keys():
                glob_taxlev2ogs[tax]['mems'].update(set(og_info_recovery[og]))

                ogs_up = glob_og_info[og]['ogs_up'].split(',')
                for og_up in ogs_up:

                    if og_up != '-':
                        lca_tax = glob_og_info[og_up]['lca_dup']
                        glob_taxlev2ogs[lca_tax]['mems'].update(set(glob_og_info[og_up]['ogs_mems']))

    return(glob_taxlev2ogs)

def make_tarfile(output_filename, source_dir):

    """
        Create tar.gz file
    """

    with tarfile.open(output_filename, "w:gz") as tar:
        tar.add(source_dir, arcname=os.path.basename(source_dir))
