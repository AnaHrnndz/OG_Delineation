from ete4.smartview import Layout
import ogd.basic_layouts as bl
import ogd.emapper_layouts as el

evol_event_layout = Layout(name='Evolutionary events', draw_node = bl.draw_node_evoltype)
scinames_layout = Layout(name='Scientific names', draw_node = bl.draw_node_leafname)
species_overlap_layout = Layout(name='Species overlap', draw_node = bl.draw_node_species_overlap)
branch_legth_layout = Layout(name='Branch lenght', draw_node = bl.draw_node_branch_lenght)
support_layout = Layout(name='Support', draw_node = bl.draw_node_support)
og_box_layout = Layout(name='OGs brackground', draw_node = bl.draw_node_background_og)
tree_style_layout = Layout('Tree style', draw_tree=bl.draw_tree_eggnog)
lca_layout = Layout('LCA', draw_node=bl.draw_node_lca_rects)

ko_layout = Layout('kegg_KO', draw_node=el.draw_kegg_ko)
kpath_layout = Layout('kegg_pathway', draw_node=el.draw_kegg_path)
pname_layout = Layout('Pref_name', draw_node=el.draw_pref_name)
Basal_OG_layout = Layout('Basal_OG', draw_node=el.draw_basal_og)



all_layouts = [evol_event_layout, scinames_layout,species_overlap_layout, branch_legth_layout, support_layout, og_box_layout, tree_style_layout, lca_layout,   
    ko_layout, kpath_layout, pname_layout, Basal_OG_layout]