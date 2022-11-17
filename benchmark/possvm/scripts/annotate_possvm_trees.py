from ete4 import NCBITaxa, Tree
from collections import defaultdict, Counter
import glob, os

ncbi = NCBITaxa()

tax_code2sci_name = defaultdict()
with open('/data/projects/og_delineation/benchmark/possvm/original_data/taxon_sampling.csv') as f:
    for line in f:
        info = line.split('\t')
        tax_code = info[0]
        sci_name = info[1]
        tax_code2sci_name[tax_code] = sci_name
    


sci_name2taxid = defaultdict()
for tax_code, sci_name in tax_code2sci_name.items():
    taxid_dict = ncbi.get_name_translator([sci_name])
    
    if taxid_dict:
        sci_name2taxid[sci_name] = str(taxid_dict[sci_name][0])
    else:
        if sci_name == 'Fungia sp.':
            taxid = '2635559'
            sci_name2taxid[sci_name] = taxid
        elif sci_name == 'Hoilungia hongkongensis H13':
            taxid = '1964320'
            sci_name2taxid[sci_name] = taxid
        elif sci_name == 'Trichoplax adhaerens H1':
            taxid = '10228'
            sci_name2taxid[sci_name] = taxid




#RefOGs HomeoDB
seq2refog = defaultdict(dict)
for file in glob.glob("/data/projects/og_delineation/benchmark/possvm/original_data/*_all_ite_mcl.csv"):
    class_name = file.split('.')[1].split('_')[0].upper()
    with open(file) as f:
        for line in f:
            if not line.startswith('sseqid'):
                info = line.split('\t')
                ori_seq_name = info[0]
                tax_code = ori_seq_name.split('_',1)[0]
                sci_name = tax_code2sci_name[tax_code]
                taxid = sci_name2taxid[sci_name]
                new_seq_name = taxid+'.'+ori_seq_name.split('_',1)[1]
                seq2refog[class_name][new_seq_name] = info[1].strip()


print(seq2refog)
#Possvm OGs
seq2possvm_og = defaultdict(dict)
for file in glob.glob("/data/projects/og_delineation/benchmark/possvm/original_data/*.ortholog_groups.csv"):
    class_name = os.path.basename(file).split('.')[0]
    with open(file) as f:
        for line in f:
            if not line.startswith('gene'):
                info = line.split('\t')
                ori_seq_name = info[0]
                tax_code = ori_seq_name.split('_',1)[0]
                sci_name = tax_code2sci_name[tax_code]
                taxid = sci_name2taxid[sci_name]
                new_seq_name = taxid+'.'+ori_seq_name.split('_',1)[1]
                seq2possvm_og[class_name][new_seq_name] = info[1].strip()

print(seq2possvm_og.keys())
# #Emapper
seq2Pname = defaultdict()
seq2OG_euk = defaultdict()
seq2OG_metazoa = defaultdict()
for file in glob.glob("/data/projects/og_delineation/benchmark/possvm/emapper_results/*.annotations"):
    with open(file) as f:
        for line in f:
            if not line.startswith('#'):
                info = line.split('\t')
                ori_seq_name = info[0]
                tax_code = ori_seq_name.split('_',1)[0]
                sci_name = tax_code2sci_name[tax_code]
                taxid = sci_name2taxid[sci_name]
                new_seq_name = taxid+'.'+ori_seq_name.split('_',1)[1]

                pname = info[8]
                seq2Pname[new_seq_name] = pname
                ogs = info[4].split(',')
                for og in ogs:
                    if og.split('|')[0].split('@')[1] == '2759':
                        seq2OG_euk[new_seq_name] = og.split('|')[0].split('@')[0]
                    elif og.split('|')[0].split('@')[1] == '33208':
                        seq2OG_metazoa[new_seq_name] = og.split('|')[0].split('@')[0]



# #Pfam
seq2pfam = defaultdict(dict)
for file in glob.glob("/data/projects/og_delineation/benchmark/possvm/emapper_results/*.pfam"):
    with open(file) as f:
        for line in f:
            if not line.startswith('#'):
                info = line.split('\t')
                ori_seq_name = info[0]
                tax_code = ori_seq_name.split('_',1)[0]
                sci_name = tax_code2sci_name[tax_code]
                taxid = sci_name2taxid[sci_name]
                new_seq_name = taxid+'.'+ori_seq_name.split('_',1)[1]

                dom_name = info[1]
                dom_start = info[7]
                dom_end = info[8]
                seq2pfam[new_seq_name][dom_name] = [dom_start, dom_end]





# #Load Trees

for file in glob.glob("/data/projects/og_delineation/benchmark/possvm/transform_trees/*.treefile"):
    class_name = (os.path.basename(file).split('.')[0].split('_')[1])


    file_name = os.path.basename(file)
    t = Tree(file)

    all_props = {'RefOG','Possvm_OG','Pname', 'OG_Euk', 'OG_Metazoa', 'Pfam_Doms', 'support'}

    #print(t.props)
    for l in t:
       

        if l.name in seq2refog[class_name].keys():
            refog = seq2refog[class_name][l.name]
        else:
            refog = '-'
        
        if l.name in seq2possvm_og[class_name].keys():
            possvm_og = seq2possvm_og[class_name][l.name]
        else:
            possvm_og = '-'
        
        if l.name in seq2Pname.keys():
            pname = seq2Pname[l.name]
        else:
            pname = '-'

        if l.name in seq2OG_euk.keys():
            og_euk = seq2OG_euk[l.name]
        else:
            og_euk = '-'
        
        if l.name in seq2OG_metazoa.keys():
            og_met = seq2OG_metazoa[l.name]
        else:
            og_met = '-'

        doms =  []
        if l.name in seq2pfam.keys():
            doms_arq = seq2pfam[l.name]
           
            for d_name, info in doms_arq.items():
                start = str(info[0])
                end = str(info[1])
                
                dom = d_name+'@'+start+'@'+end
                doms.append(dom)
        

    
        l.add_prop('RefOG', refog)
        l.add_prop('Possvm_OG', possvm_og)
        l.add_prop('Pname', pname)
        l.add_prop('OG_Euk', og_euk)
        l.add_prop('OG_Metazoa', og_met)
        l.add_prop('Pfam_Doms', doms)
        

    for n in t.traverse():
       
        if not n.is_leaf():
            refog_list = list()
            possvm_list = list()
            pname_list = list()
            og_euk_list = list()
            og_met_list = list()

            list_leaves = n.get_leaves()
            for l in list_leaves:
                
                refog_list.append(l.props.get('RefOG'))
                possvm_list.append(l.props.get('Possvm_OG'))
                pname_list.append(l.props.get('Pname'))
                og_euk_list.append(l.props.get('OG_Euk'))
                og_met_list.append(l.props.get('OG_Metazoa'))


            counter_refog = Counter(refog_list)
            counter_possvm = Counter(possvm_list)
            counter_pname = Counter(pname_list)
            counter_og_euk = Counter(og_euk_list)
            counter_og_met = Counter(og_met_list)

            common_refog = counter_refog.most_common(1)[0][0]
            n_most = counter_refog.most_common(1)[0][1]
            perc = int(n_most) / len(refog_list) 
            r_perc = round(perc, 2)
            if common_refog == '-' and len(counter_refog) >1:
                common_refog = counter_refog.most_common(2)[1][0]
                n_most = counter_refog.most_common(2)[1][1]
                perc = int(n_most) / len(refog_list) 
                r_perc = round(perc, 2)

            n.add_prop('RefOG', common_refog+'_'+str(r_perc))
            
            common_possvm = counter_possvm.most_common(1)[0][0]
            n.add_prop('Possvm_OG', common_possvm)

            common_pname = counter_pname.most_common(1)[0][0]
            n.add_prop('Pname', common_pname)

            common_og_euk = counter_og_euk.most_common(1)[0][0]
            n.add_prop('OG_Euk', common_og_euk)

            common_og_met = counter_og_met.most_common(1)[0][0]
            n.add_prop('OG_Metazoa', common_og_met)


            doms =  []
            first_leaf = next(n.iter_leaves())
            if first_leaf.name in seq2pfam.keys():
                doms_arq = seq2pfam[first_leaf.name]
                doms =  []
                for d_name, info in doms_arq.items():
                
                    start = info[0]
                    end = info[1]
                    dom =  d_name+'@'+start+'@'+end
                    doms.append(dom)

            n.add_prop('Pfam_Doms', doms)
            
    print(t.props)
    outfile = '/data/projects/og_delineation/benchmark/possvm/annotated_trees/annotated_'+file_name
    (t.write(format=1, outfile = outfile,properties=all_props, format_root_node = True))
