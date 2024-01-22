import re

def parse_taxid(node):

    #TODO: add argument for split gen name
    return node.name.split('.')[0]


def run_clean_properties(t):

    """
        Clean problematic characters from some properties
        Call clean_string()
    """

    # clean strings
    all_props = set()

    #if tree came from web server is  str format,
    if isinstance(t, str):
        t = PhyloTree(t)

    for n in t.traverse():
        for string in ("sci_name", "lca_node_name", "common_name"):
            prop = n.props.get(string)
            if prop:
                n.props[string] = clean_string(prop)

        all_props.update(set(n.props.keys()))


    return t, all_props


def clean_string(string):


    """
        Remove problematic characters for newick format
    """

    clean_string = re.sub(r"'|\[|\]|\=|\.|\-|\:", "", string)
    return clean_string


def check_nodes_up(node):

    """
        Find OGs in upper nodes
    """

    ogs_up = set()
    dups_up = list()
    while node.up:
        if node.up.props.get('node_is_og'):
            if not node.up.props.get('is_root'):
                ogs_up.add(node.up.props.get('name'))
                dups_up.append(node.up.up.props.get('name'))
            else:
                ogs_up.add(node.up.props.get('name'))
        node = node.up

    return ogs_up, dups_up