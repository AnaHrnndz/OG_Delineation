import sys
from collections import defaultdict
import math


outfile = open(sys.argv[2], 'w')


results_hmmscan = defaultdict(dict)
with open(sys.argv[1]) as infile:
    
    for line in infile:
        line = line+'\n'
        line = "\s".join(line.split()).replace('\s', '\t')
        info = line.split('\t')
        target = info[2]
        hit = info[0]
        score = (float(info[5]))
        
        if results_hmmscan[target]:
            if results_hmmscan[target]['eval'] < score:
                #print(target, eval_)
                results_hmmscan[target]['hit'] = hit
                results_hmmscan[target]['eval'] = score
        else:
            results_hmmscan[target]['hit'] = hit
            results_hmmscan[target]['eval'] = score


for target, info in results_hmmscan.items():
    
    outfile.write('\t'.join([target,'\t'.join(map(str, info.values()))])+'\n')


outfile.close()
        