#!python
import json
from sc3D import Embryo

import napari
from magicgui import magicgui
from pathlib import Path
from napari import Viewer
from napari.utils.colormaps import ALL_COLORMAPS


color_map_tissues = {
     5: [0.7411764705882353, 0.803921568627451, 1.0],
     6: [0.19607843137254902, 0.35294117647058826, 0.6078431372549019],
     7: [0.996078431372549, 0.6862745098039216, 0.08627450980392157],
     9: [0.7686274509803922, 0.27058823529411763, 0.10980392156862745],
    10: [0.10980392156862745, 1.0, 0.807843137254902],
    12: [0.7529411764705882, 0.4588235294117647, 0.6509803921568628],
    13: [0.9647058823529412, 0.13333333333333333, 0.1803921568627451],
    14: [0.7411764705882353, 0.43529411764705883, 0.6705882352941176],
    15: [0.9686274509803922, 0.8823529411764706, 0.6274509803921569],
    16: [1.0, 0.9803921568627451, 0.9803921568627451],
    18: [0.47058823529411764, 0.16470588235294117, 0.7137254901960784],
    20: [0.5019607843137255, 0.5019607843137255, 0.5019607843137255],
    21: [0.9803921568627451, 0.0, 0.5294117647058824],
    22: [0.5098039215686274, 0.1803921568627451, 0.10980392156862745],
    23: [0.5215686274509804, 0.4, 0.050980392156862744],
    24: [0.803921568627451, 0.1607843137254902, 0.5647058823529412],
    27: [0.6588235294117647, 0.6588235294117647, 0.6588235294117647],
    29: [0.0, 0.0, 0.5450980392156862],
    30: [0.5450980392156862, 0.2784313725490196, 0.36470588235294116],
    31: [1.0, 0.7568627450980392, 0.1450980392156863],
    32: [0.8705882352941177, 0.6274509803921569, 0.9921568627450981],
    33: [0.19607843137254902, 0.5137254901960784, 0.996078431372549],
    34: [0.9725490196078431, 0.6313725490196078, 0.6235294117647059],
    35: [0.7098039215686275, 0.9372549019607843, 0.7098039215686275],
    36: [0.1803921568627451, 0.8509803921568627, 1.0],
    39: [0.10980392156862745, 0.5137254901960784, 0.33725490196078434],
    40: [1.0, 0.6470588235294118, 0.30980392156862746],
    41: [0.8470588235294118, 0.7490196078431373, 0.8470588235294118]}

def display_embryo(viewer, embryo):
    # tissues_to_plot = embryo.all_tissues
    tissues_to_plot = [18, 21, 30, 31, 34]
    cells = sorted(embryo.all_cells)
    positions = [embryo.pos_3D[c] for c in cells]
    shown = [embryo.tissue[c] in tissues_to_plot for c in cells]
    properties = {'cells': cells}
    colors_rgb = [color_map_tissues[embryo.tissue[c]] for c in cells]
    properties['gene'] = [0 for _ in cells]

    points = viewer.add_points(positions, face_color=colors_rgb,
                               properties=properties,
                               metadata={'gene': None}, shown=shown)

    @magicgui(call_button='Select tissues',
              tissues={"widget_type": "Select",
                       'choices': [embryo.corres_tissue.get(t, '')
                                   for t in embryo.all_tissues],
                       'value': [embryo.corres_tissue.get(t, '')
                                 for t in tissues_to_plot]})
    def select_tissues(viewer: Viewer, tissues: str):
        tissue_to_num = {v:k for k, v in embryo.corres_tissue.items()}
        points = viewer.layers.selection.active
        if points is None:
            return
        tissues_to_plot = [tissue_to_num[t] for t in tissues]
        shown = [embryo.tissue[c] in tissues_to_plot for c in points.properties['cells']]
        points.shown = shown
        points.features['current_view'] = shown
        if points.metadata['gene'] is None:
            show_tissues(viewer)
        else:
            show_gene(viewer, points.metadata['gene'])

    @magicgui(call_button='Apply colormap',
              cmap={'label': 'Choose cmap',
                    'choices': ALL_COLORMAPS.keys()})
    def apply_cmap(viewer: Viewer, cmap: str):
        points = viewer.layers.selection.active
        if points is None:
            return
        if len(points.properties) == 0:
            return
        if points.face_color_mode.lower() != 'colormap':
            points.face_color = 'gene'
            points.face_color_mode = 'colormap'
        points.face_colormap = cmap
        points.refresh()

    @magicgui(call_button='Show gene',
              gene={'label': 'Choose gene'})
    def show_gene(viewer: Viewer, gene: str):
        points = viewer.layers.selection.active
        if points is None or not gene in embryo.anndata.raw.var_names:
            return
        if gene != points.metadata['gene']:
            colors = embryo.anndata.raw[:, gene].X.toarray()[:, 0]
            min_c, max_c = colors[points.shown].min(), colors[points.shown].max()
            colors = (colors-min_c)/(max_c-min_c)
            points.features['gene'] = colors
            points.metadata['gene'] = gene
        points.edge_color = 'black'
        points.face_color = 'gene'
        points.face_color_mode = 'colormap'
        points.face_contrast_limits = (0, 1)
        points.refresh()

    @magicgui(call_button='Show tissue')
    def show_tissues(viewer: Viewer):
        points = viewer.layers.selection.active
        if points is None:
            return
        if points.face_color_mode.lower() == 'colormap':
            points.face_color = [color_map_tissues[embryo.tissue[c]]
                                 for c in points.properties['cells']]
            points.face_color_mode = 'direct'
        points.metadata['gene'] = None
        points.refresh()

    @magicgui(call_button='Adjust contrast',
              min_c={'widget_type': 'FloatSlider', 'max': 1, 'min': 0, 'label': ''},
              max_c={'widget_type': 'FloatSlider', 'max': 1, 'min': 0, 'label': ''})
    def adj_int(viewer: Viewer, min_c: float=0, max_c: float=1):
        points = viewer.layers.selection.active
        if points is None:
            return
        if points.face_color_mode.upper() != 'COLORMAP':
            return
        if max_c < min_c:
            max_c, min_c = min_c, max_c
        points.face_contrast_limits = (min_c, max_c)
        points.refresh()

    @magicgui(call_button='Threshold cells',
              min_c={"widget_type": "FloatSlider", 'max': 1, 'min': 0, 'label': ''},
              max_c={"widget_type": "FloatSlider", 'max': 1, 'min': 0, 'label': ''})
    def threshold(viewer: Viewer, min_c: float=0, max_c: float=1):
        points = viewer.layers.selection.active
        if points is None:
            return
        if not hasattr(points.features, 'current_view'):
            points.features['current_view'] = points.shown.copy()

        points.shown = (points.features['current_view']&
                        (min_c<=points.features['gene'])&(points.features['gene']<=max_c))
        points.refresh()

    viewer.window.add_dock_widget(select_tissues)
    viewer.window.add_dock_widget(show_tissues)
    viewer.window.add_dock_widget(show_gene)
    viewer.window.add_dock_widget(threshold)
    viewer.window.add_dock_widget(adj_int)
    viewer.window.add_dock_widget(apply_cmap)

    return viewer

def loading_embryo():
    @magicgui(call_button='Load data',
              data_path={'label': 'h5ad file',
                         'widget_type': 'FileEdit',
                         'value': Path('data/registered.h5ad'),#.home(),
                         'filter': '*.h5ad'},
              tissue_names={'label': 'Tissue name',
                            'widget_type': 'FileEdit',
                            'value': Path('data/corresptissues.json'),#.home(),
                            'filter': '*.json'})
    def load_file(viewer: Viewer, data_path: str, tissue_names: str) -> Embryo:
        with open(tissue_names) as f:
            corres_tissues = json.load(f)
            corres_tissues = {eval(k): v for k, v in corres_tissues.items()}
        embryo = Embryo(data_path, store_anndata=True, corres_tissue=corres_tissues)
        viewer.window.remove_dock_widget('all')
        display_embryo(viewer, embryo)

    viewer = Viewer(ndisplay=3)
    viewer.window.add_dock_widget(load_file)
    napari.run()
    return viewer

if __name__ == '__main__':
    loading_embryo()