import re
from ete4 import NCBITaxa
from collections import defaultdict
import sys
import json

ncbi = NCBITaxa()
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


refog_table = sys.argv[1]
refog2members = defaultdict(set)
mems2refogs = defaultdict()
with open(refog_table) as f:
    for line in f:
        if not line.startswith('sseqid'):
            info = line.split('\t')
            ori_seq_name = info[0]
            tax_code = ori_seq_name.split('_',1)[0]
            sci_name = tax_code2sci_name[tax_code]
            taxid = sci_name2taxid[sci_name]
            new_seq_name = taxid+'.'+ori_seq_name.split('_',1)[1]

            refog = info[1].strip()
            refog2members[refog].add(new_seq_name)

            
            mems2refogs[new_seq_name] = refog

#print(refog2members)
refog_1mem = list()
for refog, mems in refog2members.items():
    if len(mems) == 1:
        refog_1mem.append(refog)

print(len(refog2members))
for r in refog_1mem:
    del refog2members[r]
print(len(refog2members))

myogs_table = sys.argv[2]
myogs2members = defaultdict(set)
myogs2refog_included = defaultdict(set)
with open(myogs_table) as f:
    for line in f:
        if not line.startswith('#'):
            info = line.split('\t')
            og_name = info[0]
            mems = set(info[11].split('|'))
            clean_mems = set()
            for m in mems:
                clean_mems.add(m.strip())
            
            myogs2members[og_name] = clean_mems
            
            for m in clean_mems:
                if m in mems2refogs.keys():
                    myogs2refog_included[og_name].add(mems2refogs[m])

            
#print(len(myogs2members))
#print(myogs2members)


possvm_og_table = sys.argv[3]
possvm2mems = defaultdict(set)
with open(possvm_og_table) as f:
    for line in f:
        if not line.startswith('gene'):
            info = line.split('\t')
            ori_seq_name = info[0]
            tax_code = ori_seq_name.split('_',1)[0]
            sci_name = tax_code2sci_name[tax_code]
            taxid = sci_name2taxid[sci_name]
            new_seq_name = taxid+'.'+ori_seq_name.split('_',1)[1]
            
            possvm_og = info[1].strip()
            possvm2mems[possvm_og].add(new_seq_name)



# Refog vs myogs 
# Refog vs possvm

refog_vs_myogs = defaultdict(dict)
refog_vs_possvm = defaultdict(dict)

fscore_refog_vs_myog = defaultdict(dict)
fscore_refog_vs_possvm = defaultdict(dict)

recall_refog_vs_myog = defaultdict(dict)
recall_refog_vs_possvm = defaultdict(dict)

prec_refog_vs_myog = defaultdict(dict)
prec_refog_vs_possvm = defaultdict(dict)

for refog, refog_mems in refog2members.items():
    if len(refog_mems) > 0:

# TP: Genes from refog that were assigned to the same possvm/myog
# FP: Genes from different refog that were assigned to the same possvm/myog
# FN: Genes from refog that were assigned to different possvm/myog

# precision: TP / TP +FP
# recall: TP / TP+FN
# F-score :  F1 = 2 * (precision * recall) / (precision + recall)


        for myog, myog_mems in myogs2members.items():

            refog_vs_myogs[refog][myog] = dict()
            refog_vs_myogs[refog][myog]['precision'] = int
            refog_vs_myogs[refog][myog]['recall'] = int
            refog_vs_myogs[refog][myog]['Fscore'] = int
            
            tp = refog_mems&myog_mems

            miss_refog = myog_mems.difference(refog_mems)
            fp_set = set()
            for s in miss_refog:
                if s in mems2refogs.keys():
                    fp_set.add(s)

            fn = refog_mems.difference(myog_mems)

            if len(tp) == 0:
                prec = 0
                recall = 0
                f = 0
            elif (len(tp)+ len(fp_set)) == 0:
                prec = 0
                recall = len(tp) / (len(tp) + len(fn))
                f = 0
            elif (len(tp) + len(fn)) == 0:    
                prec = len(tp) / (len(tp) + len(fp_set))
                recall = 0
                f = 0
            else:
            
                prec = len(tp) / (len(tp) + len(fp_set))
                recall = len(tp) / (len(tp) + len(fn))
                
                f = 2 * (prec * recall) / (prec + recall)

            
            refog_vs_myogs[refog][myog]['precision'] = prec
            refog_vs_myogs[refog][myog]['recall'] = recall
            refog_vs_myogs[refog][myog]['Fscore'] = f

            fscore_refog_vs_myog[refog][myog] = f
            recall_refog_vs_myog[refog][myog] = recall
            prec_refog_vs_myog[refog][myog] = prec
                


        
        for possvm_og, possvm_mems in possvm2mems.items():

            refog_vs_possvm[refog][possvm_og] = dict()
            refog_vs_possvm[refog][possvm_og]['precision'] = int
            refog_vs_possvm[refog][possvm_og]['recall'] = int
            refog_vs_possvm[refog][possvm_og]['Fscore'] = int

            tp = refog_mems&possvm_mems

            miss_refog = possvm_mems.difference(refog_mems)
            fp_set = set()
            for s in miss_refog:
                if s in mems2refogs.keys():
                    fp_set.add(s)

            fn = refog_mems.difference(possvm_mems)

            if len(tp) == 0:
                prec = 0
                recall = 0
                f = 0
            elif (len(tp)+ len(fp_set)) == 0:
                prec = 0
                recall = len(tp) / (len(tp) + len(fn))
                f = 0
            elif (len(tp) + len(fn)) == 0:    
                prec = len(tp) / (len(tp) + len(fp_set))
                recall = 0
                f = 0
            else:
            
                prec = len(tp) / (len(tp) + len(fp_set))
                recall = len(tp) / (len(tp) + len(fn))
                
                f = 2 * (prec * recall) / (prec + recall)
            
            refog_vs_possvm[refog][possvm_og]['precision'] = prec
            refog_vs_possvm[refog][possvm_og]['recall'] = recall
            refog_vs_possvm[refog][possvm_og]['Fscore'] = f

            fscore_refog_vs_possvm[refog][possvm_og] = f
            recall_refog_vs_possvm[refog][possvm_og] = recall
            prec_refog_vs_possvm[refog][possvm_og] = prec

outdir = sys.argv[4]            
        
with open(outdir+'/ogd_fscore.json', 'w') as f:
    json.dump(fscore_refog_vs_myog, f)

with open(outdir+'/ogd_recall.json', 'w') as f:
    json.dump(recall_refog_vs_myog, f)

with open(outdir+'/ogd_prec.json', 'w') as f:
    json.dump(prec_refog_vs_myog, f)



best_fscore_my_og = defaultdict(list)
for ref_og, results in fscore_refog_vs_myog.items():
    if len(refog2members[ref_og]) >0:   
        best_fval = 0
        best_name = str()
        for my_og, fval in results.items():
            if fval > best_fval:
                best_fval = fval
                best_name = my_og
        best_fscore_my_og[ref_og] =[ best_name, str(best_fval).replace('.', ',')]
        #print(ref_og, best_fscore_my_og[ref_og])

with open(outdir+'/best_ogd_fscore.tsv', 'w') as f:
    for rname, best_info in best_fscore_my_og.items():
        f.write('\t'.join(map(str,[rname,'\t'.join(best_info), '\n'])))



not_recover_refog = set()
best_recall_my_og = defaultdict(list)
for ref_og, results in recall_refog_vs_myog.items():
    if len(refog2members[ref_og]) >0:   
        best_fscore_name = best_fscore_my_og[ref_og][0]
        if best_fscore_name == '':
            not_recover_refog.add(ref_og)
        else:
            best_recall_value = str(results[best_fscore_name])
    
    
    # best_recall = 0
    # best_name = str()
    # for my_og, recall in results.items():
    #     if recall > best_recall:
    #         best_recall = recall
    #         best_name = my_og
            best_recall_my_og[ref_og] =[ best_fscore_name, best_recall_value.replace('.', ',')]

with open(outdir+'/best_ogd_recall.tsv', 'w') as f:
    for rname, best_info in best_recall_my_og.items():
        f.write('\t'.join(map(str,[rname,'\t'.join(best_info), '\n'])))





best_prec_my_og = defaultdict(list)
for ref_og, results in prec_refog_vs_myog.items():
    best_fscore_name = best_fscore_my_og[ref_og][0]
    if best_fscore_name == '':
            not_recover_refog.add(ref_og)
    else:
        best_prec_value = str(results[best_fscore_name])
    # best_prec = 0
    # best_name = str()
    # for my_og, prec in results.items():
    #     if prec > best_prec:
    #         best_prec = prec
    #         best_name = my_og
        best_prec_my_og[ref_og] =[ best_fscore_name, best_prec_value.replace('.', ',')]

with open(outdir+'/best_ogd_prec.tsv', 'w') as f:
    for rname, best_info in best_prec_my_og.items():
        f.write('\t'.join(map(str,[rname,'\t'.join(best_info), '\n'])))

    






with open(outdir+'/posvvm_fscore.json', 'w') as f:
    json.dump(fscore_refog_vs_possvm, f)


with open(outdir+'/posvvm_recall.json', 'w') as f:
    json.dump(recall_refog_vs_possvm, f)

with open(outdir+'/posvvm_prec.json', 'w') as f:
    json.dump(prec_refog_vs_possvm, f)




fscore_best_possvm = defaultdict(list)
for ref_og, results in fscore_refog_vs_possvm.items():
    best_fval = 0
    best_name = str()
    for possvm_og, fval in results.items():
        if fval > best_fval:
            best_fval = fval
            best_name = possvm_og
    fscore_best_possvm[ref_og] =[best_name, str(best_fval).replace('.', ',')]
    #print(ref_og, fscore_best_possvm[ref_og])


with open(outdir+'/best_possvm_fscore.tsv', 'w') as f:
    for rname, best_info in fscore_best_possvm.items():
        f.write('\t'.join(map(str,[rname,'\t'.join(best_info), '\n'])))




recall_best_possvm = defaultdict(list)
for ref_og, results in recall_refog_vs_possvm.items():
    
    best_name = fscore_best_possvm[ref_og][0]
    best_recall = results[best_name]
    # for possvm_og, recall in results.items():
        # if recall > best_recall:
            # best_recall = recall
            # best_name = possvm_og
    recall_best_possvm[ref_og] =[best_name, str(best_recall).replace('.', ',')]

with open(outdir+'/best_possvm_recall.tsv', 'w') as f:
    for rname, best_info in recall_best_possvm.items():
        f.write('\t'.join(map(str,[rname,'\t'.join(best_info), '\n'])))




prec_best_possvm = defaultdict(list)
for ref_og, results in prec_refog_vs_possvm.items():
    
    best_name = fscore_best_possvm[ref_og][0]
    best_prec = results[best_name]
    # for possvm_og, prec in results.items():
        # if prec > best_prec:
            # best_prec = prec
            # best_name = possvm_og
    prec_best_possvm[ref_og] =[best_name, str(best_prec).replace('.', ',')]

with open(outdir+'/best_possvm_prec.tsv', 'w') as f:
    for rname, best_info in prec_best_possvm.items():
        f.write('\t'.join(map(str,[rname,'\t'.join(best_info), '\n'])))

#
# print(not_recover_refog)

