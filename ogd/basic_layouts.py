
from ete4.smartview import Layout, TextFace, RectFace, BoxFace
from collections import  OrderedDict, defaultdict
import json
from pathlib import Path


json_path_taxid = Path(__file__).parent / "data/smartview/egg7_color_taxid.json"
json_path_sciname = Path(__file__).parent / "data/smartview/egg7_color_sciname.json"

with open(json_path_taxid, "r") as f:

    colors_taxid = json.load(f)

with open(json_path_sciname, "r") as f:
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
            color = colors_taxid.get(lca, 'orange')
            # OG root: el más grande + borde oscuro para que resalte sobre el fondo del mismo color
            return {'dot': {'shape': 'hexagon', 'radius': 5, 'fill': color, 'opacity': 0.75}}

        if node.props.get('evoltype_2') == 'S':
            # especiación: círculo pequeño (son mayoría, conviene no recargar)
            return {'dot': {'radius': 3, 'fill': 'blue'}}
        elif node.props.get('evoltype_2') == 'D':
            # duplicación: círculo mediano con borde
            return {'dot': {'radius': 3, 'fill': 'red'}}
        elif node.props.get('evoltype_2') == 'FD':
            return {'dot': {'radius': 3, 'fill': 'Coral'}}
    

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



def _soften(color, amount=0.6):
    """Blend a hex color toward white (0 = original, 1 = white) for a softer tone."""
    if not (isinstance(color, str) and color.startswith('#') and len(color) == 7):
        return color  # named colors (e.g. 'LightGrey') are already soft: leave as-is
    r, g, b = (int(color[i:i+2], 16) for i in (1, 3, 5))
    r = round(r + (255 - r) * amount)
    g = round(g + (255 - g) * amount)
    b = round(b + (255 - b) * amount)
    return f'#{r:02x}{g:02x}{b:02x}'



def draw_node_background_og(node):

    if node.props.get('monophyletic_og'):
        lca = node.props.get('lca_node_name')
        
        color = colors_sciname.get(lca, 'orange')

        return {'box': {'fill': color, 'fill-opacity': 0.35,
                        'stroke': color, 'stroke-width': 1.5, 'stroke-opacity': 0.9,
                        'rx': 8, 'ry': 8}}
       
            
    

def draw_node_lca_rects(node, collapsed):
    if node.props.get('lca_node_name'):
        lca = node.props.get('lca_node_name')
        color = colors_sciname.get(lca, 'grey')
        lca_face = TextFace(lca, rotation=90, style={'fill': 'black'})
        level = get_level(node)+10
        return [ RectFace(wmax= 30, style={'fill': color, 'stroke': 'grey'}, column=level, text=lca_face, position = 'aligned') ]



