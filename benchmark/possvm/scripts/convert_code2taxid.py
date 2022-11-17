from ete4 import NCBITaxa, Tree, SeqGroup
from collections import defaultdict
import glob, os
import re

ncbi = NCBITaxa()

tax_code2sci_name = defaultdict()
with open('/data/projects/og_delineation/benchmark/possvm/original_data/taxon_sampling.csv') as f:
    for line in f:
        info = line.split('\t')
        tax_code = info[0]
        sci_name = info[1]
        tax_code2sci_name[tax_code] = sci_name
    
print(tax_code2sci_name)

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

print(sci_name2taxid)



#os.chdir("/data/projects/og_delineation/benchmark/possvm/original_data")
for file in glob.glob("/data/projects/og_delineation/benchmark/possvm/original_data/*.treefile"):
    file_name = os.path.basename(file)
    t = Tree(file)

    for l in t:
        tax_code = l.name.split('_',1)[0]
        seq_name = l.name.split('_',1)[1]
        sci_name = tax_code2sci_name[tax_code]
        taxid = sci_name2taxid[sci_name]
        #print(l.name, tax_code, sci_name, taxid)
        seq_name_clean = re.sub('-', '_', seq_name)

        l.name = taxid+'.'+seq_name_clean

    #print(t)
    outfile = '/data/projects/og_delineation/benchmark/possvm/transform_trees/transform_'+ file_name
    t.write(outfile = outfile)

for file in glob.glob("/data/projects/og_delineation/benchmark/possvm/original_data/*.fasta"):
    file_name = os.path.basename(file)
    aln = SeqGroup(file)
    with open("/data/projects/og_delineation/benchmark/possvm/transform_alg/tranform_"+file_name, 'w') as f_out:
        for num, (s_name, seq, _) in enumerate(aln):
            tax_code = s_name.split('_',1)[0]
            seq_name = s_name.split('_',1)[1]
            seq_name_clean = re.sub('-', '_', seq_name)
            sci_name = tax_code2sci_name[tax_code]
            taxid = sci_name2taxid[sci_name]
            #print(l.name, tax_code, sci_name, taxid)

            new_name = taxid+'.'+seq_name_clean
            f_out.write('>'+new_name+'\n'+seq+'\n')


   