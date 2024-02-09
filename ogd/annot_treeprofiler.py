
import subprocess

def run_annot_treeprofiler(t, path2emmaper, path2alg, path_out):

    t.write(outfile=path_out+'tmp2.nw')

    # Run treerprofiler
    subprocess.run(("treeprofiler annotate --tree %s'/tmp2.nw' --input-type newick  --internal-parser name \
        --emapper-annotations %s/result_emapper.emapper.annotations --emapper-pfam %s/result_emapper.emapper.hmm_hits \
        -o %s --alignment %s" %(path_out, path2emmaper, path2emmaper, path_out, path2alg)), shell = True)

    # # Load new treeprofiler annotated tree
    # t_treeprof = PhyloTree(path_out)
    input()
    return t_treeprof
