from ete4.smartview import Layout, TextFace, RectFace, SeqFace, BASIC_LAYOUT
from ete4.smartview.faces import Face
import ete4.smartview.graphics as gr
from ete4.smartview.coordinates import Size, Box, make_box

import json
from pathlib import Path
from io import StringIO
import re


DOMAIN2COLOR = Path(__file__).parent / "./pfam2color.json"

def get_colormap():
    with open(Path(__file__).parent / DOMAIN2COLOR) as handle:
        _pfam2color = json.load(handle)
    return _pfam2color

colormap = get_colormap()

class SeqMotifFaceNew(Face):
    """A face for visualizing sequence motifs as rounded rectangles with connecting lines."""

    def __init__(self, rects, len_alg=None, wmax=400, hmax=30, 
                seq_format='()', box_corner_radius=10, 
                box_opacity=0.7, box_storke_width=0.1,
                fgcolor='black', bgcolor='#bcc3d0', 
                gap_color='grey', gap_linewidth=0.5,
                font_color='black', max_fsize=12, ftype='sans-serif',
                position='aligned', column=0, anchor=None):
        """
        :param rects: List of motif regions, each defined as (start, end, color, label).
        :param wmax: Total width of the sequence visualization.
        :param hmax: Height of the motif boxes.
        :param corner_radius: Radius for rounded corners.
        :param gap_color: Color of the connecting lines.
        :param gap_linewidth: Width of the connecting lines.
        :param position: Position of the face in the layout.
        :param column: Column index for alignment.
        :param anchor: Alignment anchor.
        """
        super().__init__(position, column, anchor)
        self.rects = rects
        self.len_alg = len_alg
        self.wmax = wmax
        self.hmax = hmax
        self.seq_format = seq_format
        self.box_corner_radius = box_corner_radius
        self.box_opacity = 0.7
        self.box_storke_width = 0.1
        self.fgcolor = fgcolor
        self.bgcolor = bgcolor
        self.gap_color = gap_color
        self.gap_linewidth = gap_linewidth
        self.font_color = font_color
        self.max_fsize = max_fsize
        self.ftype = ftype
        
        #self.triangles = {'^': 'top', '>': 'right', 'v': 'bottom', '<': 'left'}

        #self.gaps = self._compute_gaps(rects, wmax)

    def _compute_gaps(self, rects, seq_start, seq_end):
        """Compute the gaps (empty spaces between domain)."""
        gaps = []
        prev_end = seq_start

        for start_x, end_x, _, _, _, _, _, _ in rects:
            if start_x > prev_end:
                gaps.append((prev_end, start_x))
            prev_end = end_x

        if prev_end < seq_end:
            gaps.append((prev_end, seq_end))

        return gaps

    def draw(self, nodes, size, collapsed, zoom=(1, 1), ax_ay=(0, 0), r=1):
        dx, dy = size
        zx, zy = zoom
        
        # Compute width and height with respect to zoom 
        if dx <= 0:
            w = self.wmax
        
        h = min(zy * r * dy, self.hmax) if dy > 0 else self.hmax
        
        graphics = []

        # Normalize per-node sequence range
        # Sort rects by start
        sorted_rects = sorted(self.rects, key=lambda x: x[0])
        if not sorted_rects:
            return graphics, Size(w, h / (r * zy))

        
        #seq_start = min(start for start, _, _, _, _, _, _, _ in sorted_rects)
        seq_start = 0
        max_domain_end = max(end for _, end, *_ in sorted_rects)
        seq_end = max(max_domain_end, self.len_alg or w)
        
        seq_length = seq_end - seq_start
        
        if seq_length == 0:
            seq_length = 1  # avoid division by zero

        scale = w / seq_length

        # Compute gaps per draw call
        gaps = self._compute_gaps(sorted_rects, seq_start, seq_end)
        
        # Draw connecting lines in the gaps
        for gap_start, gap_end in gaps:
            
            gap_width = (gap_end - gap_start) * scale

            if gap_width > 0:
                size_obj = Size(gap_width, h / (r * zy))
                line_height = h / (r * zy) / 2
                box = make_box(((gap_start - seq_start) * scale, line_height), size_obj)
                style = {
                    'stroke': self.gap_color, 
                    'stroke-width': self.gap_linewidth
                    }
                graphics.append(gr.draw_line((box.x, box.y), (box.x + box.dx, box.y), style))

        # Draw rectangles
        for start_x, end_x, seq_format, poswidth, heigh, fgcolor, bgcolor, label in sorted_rects:
            
            # config size for rect
            if poswidth:
                rect_width = float(poswidth) * scale
            else:
                rect_width = (end_x - start_x) * scale

            # config height for rect
            if heigh:
                h = float(heigh)
            
            box_x = (start_x - seq_start) * scale
            size_obj = Size(rect_width, h / (r * zy))
            box = make_box((box_x, 0), size_obj)
            
            if seq_format:
                if seq_format == '()':
                    box_corner_radius = 10
                elif seq_format == '[]':
                    box_corner_radius = 0
            else:
                seq_format = self.seq_format
                box_corner_radius = self.box_corner_radius
            
            graphics.append(gr.draw_rect(box, {
                'fill': fgcolor,
                'opacity': self.box_opacity,
                'stroke': bgcolor,
                'stroke-width': self.box_storke_width,
                'rx': box_corner_radius,
                'ry': box_corner_radius
            }))

            # Check if label fits
            if label:
                try:
                    ftype, fsize, fcolor, text = label.split("|")
                    fsize = int(fsize)
                except:
                    ftype = self.ftype
                    fsize = self.max_fsize
                    fcolor = self.font_color
                    text = label

                text_style = {
                    'fill': fcolor, 
                    'font-family': ftype,
                    'text-anchor': 'middle',
                    }
                text_element = gr.draw_text(box, (0.5, 0.5), text, fs_max=fsize, rotation=0, style=text_style)
                graphics.append(text_element)

        return graphics, Size(w, h / (r * zy))

def draw_pfam_domains(tree, len_alg=None):
    def parse_dom_arq_string(dom_arq_string):
        rects = []
        if not dom_arq_string:
            return rects
        for dom_arq in dom_arq_string.split("||"):
            parts = dom_arq.split("@")
            if len(parts) == 3:
                domain_name, start, end = parts
                color = colormap.get(domain_name, "lightgray")
                rects.append([
                    int(start), int(end), 
                    "()", None, None, color, color, 
                    f"arial|15|black|{domain_name}"
                ])
        return rects

    # Get the maximum sequence length end from all nodes
    if len_alg is None:
        for node in tree.traverse():
            rects = parse_dom_arq_string(node.props.get("dom_arq", None))
            for r in rects:
                if len(r) >= 2:
                    len_alg = max(len_alg, r[1])
    
    def draw_node(node, collapsed):
        rects = parse_dom_arq_string(node.props.get("dom_arq", None))
        if rects:
            if collapsed:
                yield SeqMotifFaceNew(rects, len_alg=len_alg)
            elif node.is_leaf:
                yield SeqMotifFaceNew(rects, len_alg=len_alg)

    return draw_node


def draw_kegg_ko(node, collapsed):
    
    kko = node.props.get('kegg_ko', '')
    
    if node.is_leaf or collapsed:

        return [TextFace(kko, style={'fill': 'purple'}, column=3, position='aligned'),
                RectFace(wmax= 5, style={'fill': 'white'},  column=4, position = 'aligned')]

def draw_kegg_path(node, collapsed):
    
    kpath = node.props.get('kegg_path', '')
    
    if node.is_leaf or collapsed:

        return [TextFace(kpath, style={'fill': 'purple'}, column=5, position='aligned'),
                RectFace(wmax= 5, style={'fill': 'white'},  column=6, position = 'aligned')]


def draw_pref_name(node, collapsed):
    
    pref_name = node.props.get('pref_name', '')
    
    if node.is_leaf or collapsed:

        return [TextFace(pref_name, style={'fill': 'purple'}, column=7, position='aligned'),
                RectFace(wmax= 5, style={'fill': 'white'},  column=8, position = 'aligned')]


def draw_basal_og(node, collapsed):
    
    basal_og = node.props.get('basal_og', '')
    
    if node.is_leaf or collapsed:

        return [TextFace(basal_og, style={'fill': 'purple'}, column=9, position='aligned'),
                RectFace(wmax= 5, style={'fill': 'white'},  column=10, position = 'aligned')]