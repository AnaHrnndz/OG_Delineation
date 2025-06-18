
from ete4.smartview import Layout, TextFace, RectFace, BoxFace, BoxFace
from collections import  OrderedDict, defaultdict
import json


with open("./ogd/egg7_color_taxid.json", "r") as f:
    colors_taxid = json.load(f)

with open("./ogd/egg7_color_sciname.json", "r") as f:
    colors_sciname = json.load(f)



def get_level(node, level=1):
    if node.is_root:
        return level+1
    else:
        return get_level(node.up, level +1)


def draw_tree_eggnog(tree):
    yield  {
    'collapsed': {'shape': 'outline'}
    }


def draw_node_leafname(node, collapsed):
    
    if node.is_leaf:

        sci_name = node.props.get('sci_name')
        name_seq = node.name.split('.',1)[1]

        return [TextFace(sci_name, style={'fill': 'black'},
                         column=0, position='right'),
                TextFace(name_seq, style={'fill': 'grey'},
                         column=1, position='right')]
    if collapsed:
        text = node.props.get('lca_node_name')
        return [TextFace(text, style={'fill': 'black'},position="right", column=1)]


def draw_node_evoltype(node):

        if node.props.get('monophyletic_og'):

            lca = node.props.get('lca_dup')
            color = colors_taxid.get(lca,'orange')
            
            return {'dot': {'shape': 'square', 'radius': 4, 'fill': color } }
                

        if node.props.get('evoltype_2') == 'S':
            return {'dot': {'radius': 4, 'fill': 'blue' } }
        elif node.props.get('evoltype_2') == 'D':
            return {'dot': {'radius': 4, 'fill': 'red' } }
        elif node.props.get('evoltype_2') == 'FD':
            return {'dot': {'radius': 4, 'fill': 'Coral' } }
    

def draw_node_species_overlap(node):
    
    if node.props.get('so_score', '0.0'):
        so = str(round(float(node.props.get('so_score', '0.0')),3))
        return [TextFace( so, style={'fill': 'green'}, position = "top", column = 0, fs_min=8, fs_max=10)]
       

def draw_node_branch_lenght(node):
    
    dist = str(round(float(node.props.get('dist', '0.0')),3))
    return [TextFace( dist, style={'fill': 'grey'}, position = "bottom", column = 0, fs_min=8, fs_max=10)]


def draw_node_support(node):
    
    support = str(round(float(node.props.get('support', '0.0')),3))
    return [TextFace( support, style={'fill': 'red'}, position = "bottom", column = 0, fs_min=8, fs_max=10)]


def draw_node_background_og(node):

    if node.props.get('monophyletic_og'):
        lca = node.props.get('lca_node_name')
        
        color = colors_sciname.get(lca, 'orange')
        return {'box': {'fill':  color } }
            
    

def draw_node_lca_rects(node, collapsed):
    if node.props.get('lca_node_name'):
        lca = node.props.get('lca_node_name')
        color = colors_sciname.get(lca, 'grey')
        lca_face = TextFace(lca, rotation=90, style={'fill': 'black'})
        level = get_level(node)+7
        return [ RectFace(wmax= 30, style={'fill': color, 'stroke': 'grey'}, column=level, text=lca_face, position = 'aligned') ]