from collections import defaultdict
import json


tale_possvm2seqs = defaultdict(set)
with open('/data/projects/og_delineation/benchmark/possvm/original_data/TALE.possom.ortholog_groups.csv') as f:
    for num, (line) in enumerate(f): 
        if num != 0:
            info = line.strip().split('\t')
            tale_possvm2seqs[info[1]].add(info[0])

tale_refog2seqs = defaultdict(set)
with open('/data/projects/og_delineation/benchmark/possvm/original_data/reference.tale_all_ite_mcl.csv') as f:
    for num, (line) in enumerate(f): 
        if num != 0:
            info = line.strip().split('\t')
            tale_refog2seqs[info[1]].add(info[0])

tale_fpI = defaultdict()
for pog, poss_seqs in tale_possvm2seqs.items():
    hit = set()
    for refog, ref_seqs in tale_refog2seqs.items():
        if len(poss_seqs&ref_seqs) >0:
            hit.add(refog)
    if len(hit) == 0:
        tale_fpI[pog] = len(poss_seqs)
        print (pog+'\t'+str(len(poss_seqs)))
    # else:
        # print(pog, hit)

json.dump(tale_fpI, open('tale_fpI.json', 'w'))
print('**********************************************')
# ###############################################3

antp_possvm2seqs = defaultdict(set)
with open('/data/projects/og_delineation/benchmark/possvm/original_data/ANTP.possom.ortholog_groups.csv') as f:
    for num, (line) in enumerate(f): 
        if num != 0:
            info = line.strip().split('\t')
            antp_possvm2seqs[info[1]].add(info[0])

antp_refog2seqs = defaultdict(set)
with open('/data/projects/og_delineation/benchmark/possvm/original_data/reference.antp_all_ite_mcl.csv') as f:
    for num, (line) in enumerate(f): 
        if num != 0:
            info = line.strip().split('\t')
            antp_refog2seqs[info[1]].add(info[0])

antp_fpI = defaultdict()
for pog, poss_seqs in antp_possvm2seqs.items():
    hit = set()
    
    for refog, ref_seqs in antp_refog2seqs.items():
        if len(poss_seqs&ref_seqs) >0:
            hit.add(refog)
            #refog2poss[refog].add(pog)
    
    
    if len(hit) == 0:
        antp_fpI[pog] = len(poss_seqs)
        print (pog+'\t'+str(len(poss_seqs)))
    # elif len(hit) > 1:
        # print(pog, hit)

json.dump(antp_fpI, open('antp_fpI.json', 'w'))

# for refog, possog in refog2poss.items():
    # if len(possog) >1:
        # print(refog, possog)
print('**********************************************')

#######################################################


prd_possvm2seqs = defaultdict(set)
with open('/data/projects/og_delineation/benchmark/possvm/original_data/PRD.possom.ortholog_groups.csv') as f:
    for num, (line) in enumerate(f): 
        if num != 0:
            info = line.strip().split('\t')
            prd_possvm2seqs[info[1]].add(info[0])

prd_refog2seqs = defaultdict(set)
with open('/data/projects/og_delineation/benchmark/possvm/original_data/reference.prds_all_ite_mcl.csv') as f:
    for num, (line) in enumerate(f): 
        if num != 0:
            info = line.strip().split('\t')
            prd_refog2seqs[info[1]].add(info[0])
prd_fpI = defaultdict(set)
for pog, poss_seqs in prd_possvm2seqs.items():
    hit = set()
    
    for refog, ref_seqs in prd_refog2seqs.items():
        if len(poss_seqs&ref_seqs) >0:
            hit.add(refog)
           # refog2poss[refog].add(pog)
    
    
    if len(hit) == 0:
       print (pog+'\t'+str(len(poss_seqs)))
       prd_fpI[pog] = len(poss_seqs)
    # elif len(hit) > 1:
        # print(pog, hit)
json.dump(prd_fpI, open('prds_fpI.json', 'w'))
# for refog, possog in refog2poss.items():
    # if len(possog) >1:
        # print(refog, possog)
