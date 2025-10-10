from ete4 import  PhyloTree, NCBITaxa, GTDBTaxa
from collections import defaultdict
import os
import json




def load_tree_local(tree=None, taxonomy = None, sp_delimitator = None):

    """
        Load Tree from file
    """
    
    mssg1 = f"""        -Load tree: {os.path.basename(tree)}"""
    print(mssg1)
  
    
    #t = PhyloTree(open(tree), parser = 0)
    t = PhyloTree(open(tree))
    t.resolve_polytomy()
   
    t.set_species_naming_function(lambda node: node.name.split(sp_delimitator)[0])

    return t


def load_taxonomy(taxonomy=None, user_taxonomy=None):

    """
        Load taxonomy from local server or from user file
        Local server:
            -Eggnog 5
            -Eggnog 6
            -GTDB v207
    """

    if taxonomy == 'NCBI':
        if user_taxonomy != None:
            taxonomy_db = NCBITaxa(user_taxonomy, memory = True)
        else:
            taxonomy_db = NCBITaxa(memory = True)

    elif taxonomy == 'GTDB':
        if user_taxonomy != None:
            taxonomy_db = GTDBTaxa(user_taxonomy, memory = True)
        else:
            taxonomy_db = GTDBTaxa()

    mssg2 = f"""        -Load taxonomy: {taxonomy_db}"""
    print(mssg2)
    

    return taxonomy_db


def load_reftree(rtree=None, t=None, taxonomy_db=None):

    """
        Get reference tree (species tree) from input tree or user provide it
    """
    
    if rtree != None:
        mssg3 = f"""        -Load reftree: from user"""
        print(mssg3)
        
        reftree = PhyloTree(open(rtree), parser = 0)
    else:
        mssg3 = f"""        -Load reftree: from gene tree"""
        print(mssg3)
        
        reftree = get_reftree(t, taxonomy_db)


    taxonomy_db.annotate_tree(reftree,  taxid_attr="name")
    
    return reftree


def get_reftree(t, taxonomy_db):

    """
        Create reference tree (species tree) if user do not provide it
    """


    taxid_list = t.get_species()

    reftree = taxonomy_db.get_topology(taxid_list)
    
    return reftree


def load_taxonomy_counter(reftree=None, user_taxonomy_counter=None):

    """
        Get number of species per each taxonomic level
        Get it from reference tree or user provide it
    """
    
    if user_taxonomy_counter:
        mssg5 = f"""        -Load taxonomy counter: from user"""
        print(mssg5)

        if isinstance(user_taxonomy_counter, dict):
            level2sp_mem = user_taxonomy_counter
        else:
            with open(user_taxonomy_counter) as levels:
                level2sp_mem = json.load(levels)
    else:
        mssg5 = f"""        -Load taxonomy counter: from gene tree"""
        print(mssg5)

        level2sp_mem = get_taxonomy_counter(reftree)

    return level2sp_mem


def get_taxonomy_counter(reftree, taxonomy = None):

    """
        Create Taxonomy counter if user de not provide it
    """

    level2sp_mem = defaultdict(set)
    for l in reftree:  
        
        if l.name.isdigit():
            lin = l.props.get('lineage')
        else:
            lin = l.props.get('named_lineage')
        
        for tax in lin:
            level2sp_mem[str(tax)].add(l.name)
    
    return level2sp_mem


