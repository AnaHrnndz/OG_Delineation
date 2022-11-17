from ete4 import NCBITaxa
from collections import defaultdict
import sys

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

myogs_table = sys.argv[2]
myogs2members = defaultdict(set)
myogs2refog_included = defaultdict(set)
with open(myogs_table) as f:
    for line in f:
        if not line.startswith('#'):
            info = line.split('\t')
            og_name = info[0]
            mems = set(info[10].split('|'))
            clean_mems = set()
            for m in mems:
                clean_mems.add(m.strip())
            
            myogs2members[og_name] = clean_mems
            
            for m in clean_mems:
                if m in mems2refogs.keys():
                    myogs2refog_included[og_name].add(mems2refogs[m])

            

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


best_myOGs = set()
#refog vs myogs and vs possvm
print('#', '\t'.join(['RefOG_name', 'Len_RefOG', 'Best_myOG_name', 'Best_myOG_score', 'Diff_myOG', 'Best_possvm_name', 'Best_possvm_score', 'Diff_possvm']))
for refog, refog_mems in refog2members.items():
    best_match_sens = 0.0
    best_match_spec = 0.0
    best_match_name = str()
    total_score = len(refog_mems)
    best_diff_mems = set()
    best_included = 0
    best_included_refog = 0

    for myog, myog_mems in myogs2members.items():
        #shared_mems = (len(set(refog_mems&myog_mems))/ len(refog_mems)) / len(myog_mems)
        sens = round((len(set(refog_mems&myog_mems))/ len(refog_mems)),3) 
        spec = round((len(set(refog_mems&myog_mems))/ len(myog_mems)), 3)
        included_refog = myogs2refog_included[myog]
        
        
        
        diff_mems = set(refog_mems.difference(myog_mems))
        
        if spec >= best_match_spec:
            if sens >= best_match_sens:
                best_match_sens= sens
                best_match_spec = spec
                best_match_name = myog
                best_diff_mems = diff_mems
                best_included_refog = len(included_refog)
        elif len(included_refog) > 0 and len(included_refog) < best_included_refog:
            best_match_sens= sens
            best_match_spec = spec
            best_match_name = myog
            best_diff_mems = diff_mems
            best_included_refog = len(included_refog)
                   
    best_myOGs.add(best_match_name)

    
    spec_2 = len(myogs2refog_included[best_match_name])
    
    best_possvm_score = 0
    best_possvm_name = str()
    diff_possvm = set

    for possvm_og, possvm_mems in possvm2mems.items():
        shared_mems = set(refog_mems&possvm_mems)
        diff_mems = set(refog_mems.difference(possvm_mems))
        if len(shared_mems) > best_possvm_score:
            best_possvm_score = len(shared_mems)
            best_possvm_name = possvm_og
            diff_possvm = diff_mems

           
    
    print('\t'.join((map(str,[refog, len(refog_mems), best_match_name, best_match_sens, best_match_spec, spec_2,len(myogs2members[best_match_name]), len(best_diff_mems), best_possvm_name, best_possvm_score, len(diff_possvm)]))))

print('#########################')
print(len(best_myOGs))
print('#########################')
#myogs vs refog and vs possm
print('#', '\t'.join(['MyOG_name', 'Len_MyOG', 'Best_RefOG_name', 'Best_RefOG_score', 'RefOG_included', 'Diff_RefOG', 'Best_possvm_name', 'Best_possvm_score', 'Possvm_included','Diff_possvm']))
for myog, myog_mems in myogs2members.items():
    
    if myog in best_myOGs:
        best_refog_score = 0
        best_refog_name = '-'
        refog_included = set()
        diff_refog = set()

        for refog, refog_mems in refog2members.items():
            shared_mems = set(myog_mems&refog_mems)
            diff_mems = set(myog_mems.difference(refog_mems))
            if len(shared_mems) >0:
                refog_included.add(refog)
            if len(shared_mems) > best_refog_score:
                best_refog_score = len(shared_mems)
                best_refog_name = refog
                diff_refog = diff_mems


        best_possvm_score = 0
        best_possvm_name = '-'
        possvm_included = set()
        diff_possvm = set()

        for possvm_og, possvm_mems in possvm2mems.items():
            shared_mems = set(myog_mems&possvm_mems)
            diff_mems = set(myog_mems.difference(possvm_mems))
            if len(shared_mems) >0:
                possvm_included.add(possvm_og)
            if len(shared_mems) > best_refog_score:
                best_possvm_score = len(shared_mems)
                best_possvm_name = possvm_og
                diff_possvm = diff_mems


        print('\t'.join((map(str,[myog, len(myog_mems),  best_refog_name, best_refog_score, refog_included, len(diff_refog), best_possvm_name, best_possvm_score, possvm_included,len(diff_possvm)]))))
