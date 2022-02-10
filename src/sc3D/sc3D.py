#!python
# This file is subject to the terms and conditions defined in
# file 'LICENCE', which is part of this source code package.
# Author: Leo Guignard (leo.guignard...@AT@...univ-amu.fr)

from collections import Counter
from itertools import combinations

import numpy as np
from matplotlib import pyplot as plt
import pandas as pd

from scipy.spatial import KDTree, Delaunay
from scipy.spatial.distance import cdist
from scipy.optimize import linear_sum_assignment
from scipy.interpolate import InterpolatedUnivariateSpline, interp1d
from scipy import stats
import seaborn as sns

from sklearn.decomposition import PCA
from sklearn import linear_model

import open3d as o3d
import anndata
import transformations as tr

class Embryo(object):
    """
    Embryo class to handle samples from 3D spatial
    single cell omics. It was initially designed with
    a specific dataset in mind but it should work
    for other kinds of datasets.
    """

    def set_zpos(self, z_space=30.):
        """
        Creates the dictionary containing
        the z position of the different beads

        Args:
            z_space (float): physical space between pucks
        """
        self.z_space = z_space
        self.z_pos = {}
        cs_conversion = {b: a*z_space for a, b in enumerate(self.all_cover_slips)}
        for c in self.all_cells:
            self.z_pos[c] = cs_conversion[self.cover_slip[c]]

    def read_csv(self, path, xy_resolution=1):
        """
        Reads and loads a 3D spatial single cell
        omics dataset from a csv file.

        Args:
            path (str): path to the csv file
            xy_resolution (float): resolution of the xy coordinates
        """
        with open(path) as f:
            lines = f.readlines()
        cell_id = 0
        self.all_cover_slips = set()
        for l in lines[1:]:
            x, y, z, cat = l.split(',')[1:]
            x = eval(x)
            y = eval(y)
            z = int(z.split('_')[-1].replace('"', ''))
            cat = eval(cat)
            
            if not cat in self.tissues_to_ignore:
                self.pos[cell_id] = np.array([x, y])*xy_resolution
                self.cover_slip[cell_id] = z
                self.tissue[cell_id] = cat
                self.cells_from_cover_slip.setdefault(z, set()).add(cell_id)
                self.cells_from_tissue.setdefault(cat, set()).add(cell_id)

                self.all_tissues.add(cat)
                self.all_cover_slips.add(z)
                cell_id += 1
        self.all_cover_slips = sorted(self.all_cover_slips)
    
    def read_anndata(self, path, xy_resolution=1,
                     genes_of_interest=None,
                     tissues_to_ignore=None,
                     store_anndata=False):
        """
        Reads and loads a 3D spatial single cell
        omics dataset from an anndata file.

        Args:
            path (str): path to the csv file
            xy_resolution (float): resolution of the xy coordinates
            genes_of_interest (list of str): list of genes to load
                genes_of_interest lists the genes that can then be
                interpolated between slices
            tissues_to_ignore (list of int): list of tissue ids that
                will be ignored. The beads that have been assigned these
                tissue types will not be loaded
            store_anndata (bool): whether or not to store the anndata
                matrix. The matrix is necessary when looking for
                differentially expressed genes
        """
        from anndata import read
        data = read(path)
        data_kept = set()
        if tissues_to_ignore is not None:
            data = data[~(data.obs['predicted.id'].astype(int).isin(tissues_to_ignore))]
            data_kept.update(np.where(~(data.obs['predicted.id'].isin(tissues_to_ignore)))[0])
        if self.nb_CS_begin_ignore is not None or self.nb_CS_end_ignore is not None:
            orig = sorted(set(data.obs['orig.ident']))
            cs_to_remove = orig[:self.nb_CS_begin_ignore] + orig[-self.nb_CS_end_ignore:]
            data = data[~(data.obs['orig.ident'].isin(cs_to_remove))]
            data_kept.update(np.where(~(data.obs['orig.ident'].isin(cs_to_remove)))[0])
        data_kept = np.array(list(data_kept))
        cell_id = 0
        ids = range(len(data))
        self.all_cells = list(ids)
        self.cell_names = dict(zip(ids,
                                   map(lambda x, y: str.split(x, y)[-1],
                                       data.obs_names, '_'*len(data))))
        self.pos = dict(zip(ids,
                            np.array(list(zip(data.obs['xcoord'],
                                              data.obs['ycoord'])))*xy_resolution))
        self.tissue = dict(zip(ids,
                               data.obs['predicted.id'].astype(int)))
        cs = list(map(lambda x, y: int(str.split(x, y)[1]),
                      data.obs['orig.ident'],
                      '_'*len(data.obs['orig.ident'])))
        self.cover_slip = dict(zip(ids, cs))
        if genes_of_interest is None:
            genes_of_interest = []
        elif genes_of_interest == 'all':
            genes_of_interest = data.var_names
        self.all_genes = list(genes_of_interest)
        if 0<len(genes_of_interest):
            self.gene_expression = dict(zip(ids, np.array(data.raw[:, self.all_genes].X.A)))
        else:
            self.gene_expression = {id_:[] for id_ in ids}

        for c, cs in self.cover_slip.items():
            self.cells_from_cover_slip.setdefault(cs, set()).add(c)
        for c, T in self.tissue.items():
            self.cells_from_tissue.setdefault(T, set()).add(c)
        self.all_cover_slips = sorted(set(self.cells_from_cover_slip))
        self.all_tissues = set(self.cells_from_tissue)
        self.data = data.raw[:, genes_of_interest].X.A
        if store_anndata:
            self.anndata = data

    def rigid_transform_2D(self, A, B):
        """
        Given to lists of paired positions, computes the rigid
        transformation that minimizes between the paired positions.
        Shamefully copied from there:
        https://github.com/nghiaho12/rigid_transform_3D

        Args:
            A (2 x n ndarray): list of 2D positions
            B (2 x n ndarray): list of 2D positions
        Returns:
            M (4x4 ndarray): resulting rigid matrix
        """
        assert A.shape == B.shape

        num_rows, num_cols = A.shape
        if num_rows != 2:
            raise Exception(f"matrix A is not 2xN, it is {num_rows}x{num_cols}")

        num_rows, num_cols = B.shape
        if num_rows != 2:
            raise Exception(f"matrix B is not 2xN, it is {num_rows}x{num_cols}")

        # find mean column wise
        centroid_A = np.mean(A, axis=1)
        centroid_B = np.mean(B, axis=1)

        # ensure centroids are 3x1
        centroid_A = centroid_A.reshape(-1, 1)
        centroid_B = centroid_B.reshape(-1, 1)

        # subtract mean
        Am = A - centroid_A
        Bm = B - centroid_B

        H = Am @ np.transpose(Bm)

        # find rotation
        U, S, Vt = np.linalg.svd(H)
        R = Vt.T @ U.T

        # special reflection case
        if np.linalg.det(R) < 0:
            Vt[1,:] *= -1
            R = Vt.T @ U.T

        t = -R @ centroid_A + centroid_B
        M = np.identity(3)
        M[:2, :2] = R
        M[:2, -1:] = t

        return M
    
    def register(self, pos_ref, pos_flo, apply=False, rigid=False):
        """
        Computes and if asked, apply the transformation that minizes the
        distances between two sets of paired points. The computed transformation
        is always linear but can be rigid (rotation+translation) or
        affine (rigid+shearing)

        Args:
            pos_ref (2 x n ndarray): list of the reference 2D positions
            pos_flo (2 x n ndarray): list of 2D positions to transform
            apply (bool): if true, on top of returning the transformation
                matrix, the function returns the transformed points.
                Default: False
            rigid (bool): if true a rigid transformation is computed
                otherwise an affine function is computed
        Returns:
            M (4 x 4 ndarray): resulting rigid matrix
            new_pos (2 x n ndarray): list of transformed `pos_flo`
                positions. Only returned if `apply` is `True`
        """
        if rigid:
            M = self.rigid_transform_2D(pos_flo.T, pos_ref.T)
        else:
            try:
                M = tr.affine_matrix_from_points(pos_flo.T, pos_ref.T)
            except Exception as e:
                M = self.rigid_transform_2D(pos_flo.T, pos_ref.T)
        if apply:
            pos = np.pad(pos_flo, ((0, 0), (0, 1)), 'constant', constant_values=1).T
            new_pos = np.dot(M, pos)[:2].T
            return(M, new_pos)
        else:
            return(M)

    def center_data(self):
        """
        Centers the dataset on 0.
        Stores the result in `self.centered_pos`
        Returns:
            (dict, int:[float, float]): a dictionnary that maps beads id to
                their centered positions
        """
        self.centered_pos = {}
        for cs, cells in self.cells_from_cover_slip.items():
            pos = np.array([self.pos[c] for c in cells])
            avg = np.mean(pos, axis=0)
            self.centered_pos.update(zip(cells, pos-avg))
        return(self.centered_pos)
    
    def get_tissue_centers(self):
        """
        Computes the center of mass of the different tissues
        within each puck. Stores the result in `self.tissue_centers`
        Returns:
            (dict puck_id:(dict (tissue_id, tissue_weight): float)):
                dictionary that maps a puck id to another dictionary.
                The second dictionary maps a tissue id and its weight
                to the center of mass of the tissue in that puck
        """
        self.tissue_centers = {}
        for cs, cells in self.cells_from_cover_slip.items():
            self.tissue_centers[cs] = {}
            tissues = {t: cells.intersection(T)
                       for t, T in self.cells_from_tissue.items()}
                       # if t in self.tissue_reg}
            tissues[-1] = cells
            for tissue, c_tissue in tissues.items():
                if len(c_tissue)>2:
                    pos = [self.centered_pos[ci] for ci in c_tissue]
                    for w in range(self.tissue_weight.get(tissue, 1)):
                        self.tissue_centers[cs][(tissue, w)] = np.mean(pos, axis=0)
        return(self.tissue_centers)
    
    def build_and_apply_trsf_matrix(self, cs_ref, cs_flo):
        """
        Prepare the data, compute and apply the transformation that
        matches two pucks.

        Args:
            cs_ref (int): id of the reference puck
            cs_flo (int): if of the floating puck (that will be transformed)
        """
        # getting the shared tissue between the consecutive coverslips
        tissues_ref = set(self.tissue_centers[cs_ref].keys())
        tissues_flo = set(self.tissue_centers[cs_flo].keys())
        tissues_common = list(tissues_ref.intersection(tissues_flo))
        # getting the average position of the tissue to register
        pos_flo = np.array([self.tissue_centers[cs_flo][t] for t in tissues_common])
        # getting the average position of the reference tissue
        pos_ref = np.array([self.tissue_centers_reg[cs_ref][t] for t in tissues_common])
        # computing the transformation
        M = self.rigid_transform_2D(pos_flo.T, pos_ref.T)
        # M = self.register(pos_flo, pos_ref)
        
        # preping the floating positions for the trsf
        pos = np.pad([self.centered_pos[ci] for ci in self.cells_from_cover_slip[cs_flo]],
                     ((0, 0), (0, 1)), 'constant', constant_values=1).T
        # applying the trsf
        new_pos = np.dot(M, pos)[:2].T
        # updating the position dictionary
        self.registered_pos.update(dict(zip(self.cells_from_cover_slip[cs_flo], new_pos)))
        
        # preping the floating tissue centers
        pos = np.pad([self.tissue_centers[cs_flo][t] for t in self.tissue_centers[cs_flo]],
                     ((0, 0), (0, 1)), 'constant', constant_values=1).T
        new_pos = np.dot(M, pos)[:2].T
        self.tissue_centers_reg[cs_flo] = dict(zip(self.tissue_centers[cs_flo], new_pos))

    def register_with_tissues(self):
        """
        Register together all the pucks using tissue center of masses.
        """
        if not hasattr(self, 'centered_pos'):
            self.center_data()
        if not hasattr(self, 'tissue_centers'):
            self.get_tissue_centers()
        cs_ref = self.all_cover_slips[0]
        self.tissue_centers_reg = {}
        self.tissue_centers_reg[cs_ref] = self.tissue_centers[cs_ref]
        self.registered_pos = {c: self.centered_pos[c] for c in self.cells_from_cover_slip[cs_ref]}
        for cs_flo in self.all_cover_slips[1:]:
            self.build_and_apply_trsf_matrix(cs_ref, cs_flo)
            cs_ref = cs_flo

    def build_pairing(self, cs1, cs2, rebuild=False, refine=False, th_d=None):
        """
        Build the pairing between beads from two pucks and stores it in the
        dictionary `pairing` that maps a bead id to the id of its paired bead.

        Args:
            cs1 (int): id of the first puck
            cs2 (int): id of the second puck
            rebuild (bool): if true the previously computed pairings are erased
                Default: False (you should probably keep it that way)
            refine (bool): if true, uses the previously computed registration to
                do the pairing (usually kept at False).
                Default: False
            th_d (bool | float): threshold above which a pairing is discarded.
                If th_d is a boolean, then the threshold is the median of the
                distribution of all the distances. If th_d is a float the value
                given is used as a threshold.
                Usually used as a float.
        Returns:
            pos_ref (2 x n ndarray): list of positions that have been paired from
                the first puck (`cs1`)
            pos_flo (2 x n ndarray): list of positions that have been paired from
                the second puck (`cs2`)
        """
        if rebuild:
            self.pairing = {}
        pos_ref = []
        pos_flo = []
        tot_paired = 0
        for tissue in self.all_tissues:
            cells_cs1 = np.array([c for c in self.cells_from_cover_slip[cs1]
                                  if self.tissue[c] == tissue])
            cells_cs2 = np.array([c for c in self.cells_from_cover_slip[cs2]
                                  if self.tissue[c] == tissue])
            positions_cs1 = np.array([self.final.get(c, self.registered_pos[c])
                                      for c in cells_cs1 if self.tissue[c] == tissue])
            if refine:
                positions_cs2 = np.array([self.pos_reg_aff[c]
                                          for c in cells_cs2 if self.tissue[c] == tissue])
            else:
                positions_cs2 = np.array([self.registered_pos[c]
                                          for c in cells_cs2 if self.tissue[c] == tissue])
            if len(positions_cs1) > 0 and len(positions_cs2) > 0:
                distance = cdist(positions_cs1, positions_cs2)
                copy_d = distance.copy()
                if isinstance(th_d, bool):
                    th_d_tissue = np.max(distance)/2
                    distance[th_d_tissue<distance] = np.inf
                elif isinstance(th_d, int) or isinstance(th_d, float):
                    th_d_tissue = th_d
                    distance[th_d_tissue<distance] = np.inf
                else:
                    th_d_tissue=np.inf
                try:
                    pairing = linear_sum_assignment(distance)
                    pos_ref += list(positions_cs1[pairing[0]])
                    pos_flo += list(positions_cs2[pairing[1]])
                    self.pairing.update(zip(cells_cs1[pairing[0]], cells_cs2[pairing[1]]))
                except:
                    pairing = linear_sum_assignment(copy_d)
                    pos_ref_tmp = positions_cs1[pairing[0]]
                    pos_flo_tmp = positions_cs2[pairing[1]]
                    distance_paired = np.linalg.norm(np.array(pos_ref_tmp)-np.array(pos_flo_tmp), 
                                                     axis=1).reshape(-1, 1)
                    to_keep = (distance_paired<th_d_tissue).reshape(-1)
                    pos_ref_tmp = pos_ref_tmp[to_keep]
                    pos_flo_tmp = pos_flo_tmp[to_keep]
                    pos_ref += list(pos_ref_tmp)
                    pos_flo += list(pos_flo_tmp)
                    self.pairing.update(zip(cells_cs1[pairing[0][to_keep]], cells_cs2[pairing[1][to_keep]]))
        return pos_ref, pos_flo

    def register_cs(self, cs1, cs2, refine=False, rigid=False, final=False, th_d=None):
        """
        Registers the puck `cs2` onto the puck `cs1`.

        Args:
            cs1 (int): id of the first puck
            cs2 (int): id of the second puck
            refine (bool): if true, uses the previously computed registration to
                do the pairing (usually kept at False).
                Default: False
            rebuild (bool): if true the previously computed pairings are erased
                Default: False (you should probably keep it that way)
            final (bool): if True assumes that it is the final registration between
                the two considered pucks (legacy, always True now).
                Default: True
            th_d (bool | float): threshold above which a pairing is discarded.
                If th_d is a boolean, then the threshold is the median of the
                distribution of all the distances. If th_d is a float the value
                given is used as a threshold.
                Usually used as a float.

        """
        if not hasattr(self, 'registered_pos'):
            self.register_with_tissues()
        if not hasattr(self, 'pos_reg_aff'):
            self.pos_reg_aff = {}
        if not hasattr(self, 'pairing'):
            self.pairing = {}
        if not hasattr(self, 'final') and final:
            self.final = {c: self.centered_pos[c] for c in self.cells_from_cover_slip[cs1]}
        pos_ref, pos_flo = self.build_pairing(cs1, cs2, rebuild=False, refine=refine, th_d=th_d)
        M = self.register(np.array(pos_ref), np.array(pos_flo), apply=False, rigid=rigid)
        cells_cs1 = self.cells_from_cover_slip[cs1]
        cells_cs2 = self.cells_from_cover_slip[cs2]
        positions_cs1 = np.array([self.registered_pos[c] for c in cells_cs1])
        if refine:
            positions_cs2 = np.array([self.pos_reg_aff[c] for c in cells_cs2])
        else:
            positions_cs2 = np.array([self.registered_pos[c] for c in cells_cs2])
        pos = np.pad(positions_cs2, ((0, 0), (0, 1)), 'constant', constant_values=1).T
        new_pos = np.dot(M, pos)[:2].T
        new_pos = pos[:2].T
        self.pos_reg_aff.update(zip(cells_cs2, new_pos))
        if final:
            self.final.update(zip(cells_cs2, new_pos))
        return M
        
    def build_gabriel_graph(self, node_ids, pos):
        """
        Build the gabriel graph of a set of nodes with
        associtated positions.

        Args:
            node_ids ([int, ] (size n)): list of node ids
            pos (n x m ndarray): ndarray of the positions where n is
                the number of nodes and m is the spatial dimension
        Returns:
            final_GG (dict id: set([ids, ])): the gabriel graph as
                an adjacency list, a dictionary that maps node ids
                to the list of neighboring node ids
        """
        tmp = Delaunay(pos)
        delaunay_graph = {}

        for N in tmp.simplices:
            for e1, e2 in combinations(np.sort(N), 2):
                delaunay_graph.setdefault(e1, set()).add(e2)
                delaunay_graph.setdefault(e2, set()).add(e1)
        Gabriel_graph = {}

        for e1, neighbs in delaunay_graph.items():
            for ni in neighbs:
                if not any([np.linalg.norm((pos[ni] + pos[e1])/2 - pos[i])<np.linalg.norm(pos[ni] - pos[e1])/2
                        for i in delaunay_graph[e1].intersection(delaunay_graph[ni])]):
                    Gabriel_graph.setdefault(e1, set()).add(ni)
                    Gabriel_graph.setdefault(ni, set()).add(e1)

        final_GG = {}
        for e1, neighbs in Gabriel_graph.items():
            neighbs = np.array(list(neighbs))
            distances = np.linalg.norm(pos[e1] - [pos[ni] for ni in neighbs], axis=1)
            final_GG[node_ids[e1]] = set([node_ids[ni] for ni in neighbs[distances<=5*np.median(distances)]])
        return final_GG

    def plot_coverslip(self, cs, pos='pos', ax=None,
                       tissues_to_plot=None, legend=False,
                       color=None, cells=None, **kwargs):
        """
        Plot a puck

        Args:
            cs (int): id of the puck to plot
            pos (str): attribute defining the positions to plot.
                Probably want to use 'final' since it is the registered
                positions. Despite that, default is 'pos', the original
                positions
            ax (matplotlib.AxesSubplot): can be provided to constrain the
                plot
            tissues_to_plot ([t_id, ]): list of tissue ids to plot
            legend (bool): if True a legend is ploted.
                Default: False
            color (dict t_id: [float, float, float]): a dictionary that
                maps a tissue id to a given color. If `None`, then the default
                matplotlib colors are used.
                Default: None
            cells ([id, ]): list of bead ids to plot. If `cells` is provided
                `tissues_to_plot` and `cs` are ignored
            kwargs : the kwargs are passed down to the matplotlib.scatterplot call
        Returns:
            fig (matplotlib.Figure): the created figure
            ax (matplotlib.AxesSubplot): the working axis
        """
        if ax is None:
            fig, ax = plt.subplots()
        else:
            fig = ax.get_figure()
        if isinstance(pos, str):
            positions_attr = self.__getattribute__(pos)
        else:
            positions_attr = pos
        if tissues_to_plot is None and cells is None:
            cells = self.cells_from_cover_slip[cs]
        elif cells is None:
            cells = [c for c in self.cells_from_cover_slip[cs] if self.tissue[c] in tissues_to_plot]
        positions = np.array([positions_attr[c] for c in cells])
        tissues = [self.tissue[c] for c in cells]
        if len(positions)<1:
            return fig, ax
        scatter_args = {'marker':'.', 's':25, 'cmap':'tab20',
                        'vmin':min(self.all_tissues), 'vmax':max(self.all_tissues)}
        scatter_args.update(kwargs)
        if color is None:
            color = tissues
            # scatter_args['c'] = color
        elif isinstance(color, dict):
            color = [color.get(t, [.8,]*3) for t in tissues]
        scatter = ax.scatter(*positions.T, c=color, **scatter_args)
        if legend:
            ax.legend(handles=scatter.legend_elements()[0], labels=np.unique(tissues))
        return fig, ax

    def removing_spatial_outliers(self, th=.2, n_components=3):
        """
        Removes spatial outliers given a threshold and a number of components

        Args:
            th (float): Likelyhood below which a bead is discarded.
                Default: 0.2
            n_components (int): number of components for the gaussian mixture
                model.
                Default: 3 (probably should keep it that way. Lesser than 2 will
                            crash things)
        """
        from sklearn import mixture
        self.filtered_cells = set()
        for t in self.all_tissues:
            c_to_d = {}
            cells_final = []
            for cs, cells in self.cells_from_cover_slip.items():
                cells_t = np.array(list(cells & self.cells_from_tissue[t]))
                if len(cells_t)<2:
                    continue
                cells_final.extend(list(cells_t))
                pos = [self.pos[c] for c in cells_t]
                kdtree = KDTree(pos)
                dist = list(kdtree.query(pos, k=2, workers=-1)[0][:, 1])
                c_to_d.update(zip(cells_t, dist))
            if len(cells_final)<10:
                continue
            cells_final = np.array(cells_final)

            D = np.array([d for c, d in c_to_d.items()])
            gmm = mixture.GaussianMixture(n_components=n_components, max_iter=1000,
                                          covariance_type='full').fit(D.reshape(-1,1))
            order = np.argsort(gmm.means_, axis=0)
            proba0 = gmm.predict_proba(D.reshape(-1, 1))[:, order[0,0]]
            proba1 = gmm.predict_proba(D.reshape(-1, 1))[:, order[1,0]]
            self.filtered_cells.update(cells_final[(th<proba0)|(th<proba1)])
        self.cells.intersection_update(self.filtered_cells)
        self.all_cells = set(self.all_cells).intersection(self.filtered_cells)
        self.pos = {k:self.pos[k] for k in self.filtered_cells}
        self.tissue = {k:self.tissue[k] for k in self.filtered_cells}
        self.cover_slip = {k:self.cover_slip[k] for k in self.filtered_cells}
        self.cell_names = {k:self.cell_names[k] for k in self.filtered_cells}
        self.gene_expression = {k:self.gene_expression[k] for k in self.filtered_cells}
        for t, c in self.cells_from_cover_slip.items():
            c.intersection_update(self.filtered_cells)
        for t, c in self.cells_from_tissue.items():
            c.intersection_update(self.filtered_cells)

    def reconstruct_intermediate(self, rigid=True,
                                 th_d=True, cs=None,
                                 multicore=True, genes=None):
        """
        Register all pucks against each other and build the interpolation splines

        Args:
            rigid (bool): if True, a rigid transformation is computed and applied.
                Otherwise it is an affine transformation.
                Default: True
            th_d (bool | float): threshold above which a pairing is discarded.
                If th_d is a boolean, then the threshold is the median of the
                distribution of all the distances. If th_d is a float the value
                given is used as a threshold. Usually used as a float.
            cs ([p_id, ]): list of puck ids to treat. If None, then all the pucks
                are treated.
                Default: None
            multicore (bool): useless at the time being. Maybe one day ...
            genes ([str, ]): gene names that will be interpolated
        """
        if not hasattr(self, 'z_pos'):
            self.set_zpos()
        disapear_bounds = (.1, .5, .9)
        if cs is not None:
            cs_to_treat = cs
        else:
            cs_to_treat = self.all_cover_slips
        self.GG_cs = {}
        self.KDT_cs = {}
        for i, cs1 in enumerate(cs_to_treat[:-1]):
            cs2 = cs_to_treat[i+1]
            M1 = self.register_cs(cs1, cs2, rigid=rigid, final=True, th_d=th_d)
        for cs in cs_to_treat:
            cids = list(self.cells_from_cover_slip[cs])
            pos = [self.final[c] for c in cids]
            self.GG_cs[cs] = self.build_gabriel_graph(cids, pos)
        paths = []
        inv_pairing = {v:k for k, v in self.pairing.items()}
        roots = set(self.pairing).difference(inv_pairing)
        all_paired = set(inv_pairing).union(self.pairing)
        for c in roots:
            p = [c]
            while p[-1] in self.pairing:
                p.append(self.pairing[p[-1]])
            paths.append(p)

        unmapped_down = set(self.all_cells) - set(inv_pairing)
        unmapped_down.difference_update(self.cells_from_cover_slip[min(self.all_cover_slips)])
        unmapped_up = set(self.all_cells).difference(self.pairing)
        unmapped_up.difference_update(self.cells_from_cover_slip[max(self.all_cover_slips)])


        self.KDT_cs_down = {}
        self.paired_cs_down = {}
        for cs in cs_to_treat[1:]:
            self.paired_cs_down[cs] = (set(self.cells_from_cover_slip[cs]) &
                                         set(inv_pairing))
            self.paired_cs_down[cs] = np.array(list(self.paired_cs_down[cs]))
            pos = [self.final[c] for c in self.paired_cs_down[cs]]
            self.KDT_cs_down[cs] = KDTree(pos)

        succ = {}
        arrival_down = {}
        d_to_closest_down = {}
        for c in unmapped_down:
            cs = self.cover_slip[c]
            neighbs = self.GG_cs[cs].get(c).difference(unmapped_down)
            if len(neighbs)<1:
                neighbs = [self.paired_cs_down[cs][self.KDT_cs_down[cs].query(self.final[c], 1)[1]]]
            arrival_down[c] = np.mean([self.final[inv_pairing[ni]] for ni in neighbs], axis=0)
            d_to_closest_down[c] = np.mean([np.linalg.norm(self.final[c] - self.final[ni])
                                            for ni in neighbs])

        self.KDT_cs_up = {}
        self.paired_cs_up = {}
        for cs in cs_to_treat[:-1]:
            self.paired_cs_up[cs] = (set(self.cells_from_cover_slip[cs]) &
                                         set(self.pairing))
            self.paired_cs_up[cs] = np.array(list(self.paired_cs_up[cs]))
            pos = [self.final[c] for c in self.paired_cs_up[cs]]
            self.KDT_cs_up[cs] = KDTree(pos)

        succ = {}
        arrival_up = {}
        d_to_closest_up = {}
        for c in unmapped_up:
            cs = self.cover_slip[c]
            neighbs = self.GG_cs[cs].get(c).difference(unmapped_up)
            if len(neighbs)<1:
                neighbs = [self.paired_cs_up[cs][self.KDT_cs_up[cs].query(self.final[c], 1)[1]]]
            arrival_up[c] = np.mean([self.final[self.pairing[ni]] for ni in neighbs], axis=0)
            d_to_closest_up[c] = np.mean([np.linalg.norm(self.final[c] - self.final[ni])
                                            for ni in neighbs])

        d_to_closest_vals = list(d_to_closest_down.values()) + list(d_to_closest_up.values())
        med_to_closest = np.median(d_to_closest_vals)
        min_to_closest = np.percentile(d_to_closest_vals, 1)
        max_to_closest = np.percentile(d_to_closest_vals, 99)
        end, mid, start = disapear_bounds
        dist_to_disapear = interp1d([min_to_closest, med_to_closest, max_to_closest],
                                    [start, mid, end], bounds_error=False, fill_value=(start, end))

        cells_to_treat = set(self.all_cells)
        all_trajs = {}
        if genes is not None and isinstance(genes, list):
            all_expr = {}
        elif not isinstance(genes, list):
            print('The genes to process have to be in a `list`')
            genes = None
            all_expr = []

        nb_skipped = 0
        while 0<len(cells_to_treat):
            curr_cell = cells_to_treat.pop()
            traj = [curr_cell]
            while traj[0] in inv_pairing:
                traj.insert(0, inv_pairing[traj[0]])
            while traj[-1] in self.pairing:
                traj.append(self.pairing[traj[-1]])
            if len(traj)<=1:
                nb_skipped += 1
                continue
            pos_traj = [self.final[c] for c in traj]
            z_traj = [self.z_pos[c] for c in traj]
            if traj[-1] in arrival_up:
                pos_traj.append(arrival_up[traj[-1]])
                D = dist_to_disapear(d_to_closest_up[traj[-1]])
                z_traj.append(z_traj[-1]+D*self.z_space)
            if traj[0] in arrival_down:
                pos_traj.insert(0, arrival_down[traj[0]])
                D = dist_to_disapear(d_to_closest_down[traj[0]])
                z_traj.insert(0, z_traj[0]-D*self.z_space)
            pos_traj_x, pos_traj_y = zip(*pos_traj)
            k_interp = min(3, len(pos_traj_x)-1)
            f_traj_x = InterpolatedUnivariateSpline(z_traj, pos_traj_x, k=k_interp, ext='const')
            f_traj_y = InterpolatedUnivariateSpline(z_traj, pos_traj_y, k=k_interp, ext='const')
            if genes is not None:
                for i, g in enumerate(genes):
                    if g in self.all_genes:
                        index = self.all_genes.index(g)
                        value_traj = [self.gene_expression[c][index] for c in traj]
                        z_traj_g = [self.z_pos[c] for c in traj]
                        k_interp = min(3, len(z_traj_g)-1)
                        f_traj_v = InterpolatedUnivariateSpline(z_traj_g,
                                                                value_traj,
                                                                k=1,
                                                                ext='const')
                    all_expr.setdefault(g, {}).update({traj[0]: [min(z_traj), max(z_traj), f_traj_v]})
            # f_traj_x = Rbf(z_traj, pos_traj_x, smooth=.5)#, function='gaussian', smooth=2)
            # f_traj_y = Rbf(z_traj, pos_traj_y, smooth=.5)#, function='gaussian', smooth=2)

            all_trajs[traj[0]] = [min(z_traj), max(z_traj), f_traj_x, f_traj_y]
            cells_to_treat -= set(traj)
        self.all_trajs = all_trajs
        self.all_expr = all_expr

    def plot_slice(self, angle, color_map=None, rot_orig=[0, 0, 1], origin=[0, 0, 0],
                   thickness=30, tissues=None, angle_unit='degree',
                   nb_interp=5, output_path=None, gene=None,
                   min_g1=None, min_g2=None, max_g1=None, max_g2=None,
                   main_bi_color='g', figsize=(5, 5), path_scale=None, **kwargs):
        """
        Plot an arbitrarly oriented slice according to an angle a direction and an origin.

        Args:
            angle (float): angle of the rotation of the slice
            color_map (matplotlib.cmap): color map that will be applied
            rot_origin ([int, int, int]): 3D vector of the normal of the
                rotation plan. If [0, 0, 1] is given the rotation will be
                around the z axis
            origin ([int, int, int]): coordinates of center of the rotation
            thickness (float): thickness of the slice
            tissues ([t_id, ]): list of tissue ids to plot
            angle_unit (str): if `'degree'` the angle is treated as degrees.
                Otherwise it is treated a radii
            nb_interp (int): number of pucks to interpolate in between
                existing pucks
            output_path (str): path to the output figure
            gene (str | [str, str]): gene name to interpolate. If a list
                of 2 strings is given gene colocalization is plotted
            min_g1/g2 (float): minimum threshold value for the first and
                second genes when colocalization. If `None`, the 2nd
                percentile of the gene expression is used as a threshold
            max_g1/g2 (float): maximum threshold value for the first and
                second genes when colocalization. If `None`, the 98th
                percentile of the gene expression is used as a threshold
            main_bi_color ('g' | 'r' | 'b'): when colocalization, two
                colors are used green and red+blue ('g'), etc
            figsize (float, float): width and height of the figure given
                to the function plt.figure
            path_scale (str): path to the figure that will contain the
                scale for colocalization figures
            kwargs : the keyword args are forwarded to the scatterplot function
        Returns:
            points_to_plot (n x 2 ndarray): list of the positions of the points
                that have been plotted
        """
        if tissues is None:
            tissues = self.all_tissues
        if angle_unit == 'degree':
            angle = np.deg2rad(angle)
        x_angle, y_angle, z_angle = angle
        rot_x = tr.rotation_matrix_py(x_angle, [1, 0, 0], origin)
        rot_y = tr.rotation_matrix_py(y_angle, [0, 1, 0], origin)
        rot_z = tr.rotation_matrix_py(z_angle, [0, 0, 1], origin)
        rot_composed = rot_x@rot_y@rot_z
        new_axis = (np.hstack([rot_orig, 1])@rot_composed)[:-1]
        equation = lambda pos: np.sum(new_axis*pos, axis=1)-origin@new_axis
        if gene is not None and not isinstance(gene, str):
            if len(gene)==1:
                gene = gene[0]
                points, color = self.produce_em(nb_interp, tissues, gene=gene)
                color = np.array(color)
            else:    
                colors = []
                for g in gene:
                    points, color = self.produce_em(nb_interp, tissues, gene=g)
                    colors.append(color)
                C = np.array(colors)
                if min_g1 is None:
                    min_g1 = np.percentile(C, 2, axis=1)[0]
                if min_g2 is None:
                    min_g2 = np.percentile(C, 2, axis=1)[1]
                if max_g1 is None:
                    max_g1 = np.percentile(C, 98, axis=1)[0]
                if max_g2 is None:
                    max_g2 = np.percentile(C, 98, axis=1)[1]
                norm = lambda C: (C-[[min_g1], [min_g2]]) / [[max_g1-min_g1], [max_g2-min_g2]]
                V = norm(C)
                V[V<0] = 0
                V[1<V] = 1
                final_C = np.zeros((len(colors[0]), 3))
                on_channel = (np.array(['r', 'g', 'b'])==main_bi_color.lower()).astype(int)
                final_C[:,0] = V[on_channel[0]]
                final_C[:,1] = V[on_channel[1]]
                final_C[:,2] = V[on_channel[2]]
                if path_scale:
                    scale_square = np.zeros((256, 256, 3))
                    V1 = np.linspace(0, max_g1, 256)
                    V2 = np.linspace(0, max_g2, 256)
                    VS = np.array([V1, V2])
                    VS = norm(VS)
                    VS[VS<0] = 0
                    VS[1<VS] = 1
                    scale_square[...,np.where(on_channel)[0][0]] = VS[0]
                    for ax in np.where(1-on_channel)[0]:
                        scale_square[...,ax] = VS[1].reshape(-1, 1)
                    # scale_square[...,np.where(1-on_channel)[0]]
                    fig, ax = plt.subplots(figsize=(5, 5))
                    ax.imshow(scale_square.swapaxes(1, 0), origin='lower')
                    recap_g1 = lambda x: x*255/max_g1
                    recap_g2 = lambda x: x*255/max_g2
                    vals_g1 = np.arange(np.floor(max_g1)+1, dtype=int)
                    vals_g2 = np.arange(np.floor(max_g2)+1, dtype=int)
                    ax.set_xticks(recap_g1(vals_g1))
                    ax.set_yticks(recap_g2(vals_g2))
                    ax.set_xticklabels(vals_g1)
                    ax.set_yticklabels(vals_g2)
                    ax.set_xlabel(gene[0])
                    ax.set_ylabel(gene[1])
                    fig.tight_layout()
                    fig.savefig(path_scale)
        else:
            points, color = self.produce_em(nb_interp, tissues, gene=gene)
            color = np.array(color)
        points = np.array(points)
        dist_to_plan = equation(points)
        plan = (np.abs(dist_to_plan)<thickness)
        dist_to_plan = dist_to_plan[plan]
        points_to_plot = points[plan]
        points_to_plot = (np.hstack([points_to_plot, [[1]]*points_to_plot.shape[0]])@rot_composed)[:, :-1]
        if gene is None:
            color_to_plot = np.array([color_map[c] for c in color[plan]])
        elif not isinstance(gene, str):
            color_to_plot = final_C[plan]
        else:
            color_to_plot = color[plan]
        p_order = np.argsort(dist_to_plan)
        points_to_plot = points_to_plot[p_order]
        color_to_plot = color_to_plot[p_order]
        fig = plt.figure(figsize=figsize)
        ax = fig.add_subplot(1, 1, 1)#, projection='3d')
        if gene is None:
            kwargs_scatter = { 's':5, 'color':color_to_plot}
        else:
            kwargs_scatter = { 's':5, 'c':color_to_plot}
        kwargs_scatter.update(kwargs)
        ax.scatter(*(points_to_plot.T[:-1]), **kwargs_scatter)
        ax.axis('equal')
        if output_path is not None:
            fig.savefig(output_path)
        return points_to_plot
    
    def ply_slice(self, file_name, angle, color_map, rot_orig=[0, 0, 1],
                  origin=[0, 0, 0], thickness=30, tissues=None,
                  tissues_colored=None, angle_unit='degree',
                  gene=None, nb_interp=5):
        """
        Build a ply file containing a slice

        Args:
            file_name (str): path to the output ply file
            angle (float): angle of the rotation of the slice
            color_map (matplotlib.cmap): color map that will be applied
            rot_origin ([int, int, int]): 3D vector of the normal of the
                rotation plan. If [0, 0, 1] is given the rotation will be
                around the z axis
            origin ([int, int, int]): coordinates of center of the rotation
            thickness (float): thickness of the slice
            tissues ([t_id, ]): list of tissue ids to plot
            angle_unit (str): if `'degree'` the angle is treated as degrees.
                Otherwise it is treated a radii
            nb_interp (int): number of pucks to interpolate in between
                existing pucks
            gene (str): gene name to interpolate
        Returns:
            points_to_plot (n x 2 ndarray): list of the positions of the points
                that have been plotted
        """
        if tissues is None:
            tissues = self.all_tissues
        if tissues_colored is None:
            tissues = self.all_tissues
        if angle_unit == 'degree':
            angle = np.deg2rad(angle)
        x_angle, y_angle, z_angle = angle
        rot_x = tr.rotation_matrix_py(x_angle, [1, 0, 0], origin)
        rot_y = tr.rotation_matrix_py(y_angle, [0, 1, 0], origin)
        rot_z = tr.rotation_matrix_py(z_angle, [0, 0, 1], origin)
        rot_composed = rot_x@rot_y@rot_z
        new_axis = (np.hstack([rot_orig, 1])@rot_composed)[:-1]
        equation = lambda pos: np.sum(new_axis*pos, axis=1)-origin@new_axis
        points, color = self.produce_em(nb_interp, tissues, gene=gene)
        points = np.array(points)
        color = np.array(color)
        plan = (np.abs(equation(points))<thickness)
        points_to_plot = points[plan]
        color_to_plot = color[plan]
        if gene is None:
            mapping = np.array([[.8, .8, .8] for t in range(max(color_map)+1)])
            for t in tissues_colored:
                mapping[t] = color_map.get(t, [.8, .8, .8])
        else:
            min_v = np.percentile(color_to_plot, 1)
            max_v = np.percentile(color_to_plot, 99)
            color_to_plot = plt.cm.viridis((color_to_plot - min_v)/(max_v - min_v))[:,:-1]
        pcd = o3d.geometry.PointCloud()
        pcd.points = o3d.utility.Vector3dVector(points_to_plot)
        if gene is None:
            pcd.colors = o3d.utility.Vector3dVector(mapping[color_to_plot])
        else:
            pcd.colors = o3d.utility.Vector3dVector(color_to_plot)
        o3d.io.write_point_cloud(file_name, pcd)
        
        return points_to_plot

    def anndata_slice(self, file_name, angle, gene_list, rot_orig=[0, 0, 1],
                      origin=[0, 0, 0], thickness=30, tissues=None,
                      angle_unit='degree', nb_interp=5):
        """
        Build a anndata file containing a slice

        Args:
            file_name (str): path to the output ply file
            angle (float): angle of the rotation of the slice
            color_map (matplotlib.cmap): color map that will be applied
            rot_origin ([int, int, int]): 3D vector of the normal of the
                rotation plan. If [0, 0, 1] is given the rotation will be
                around the z axis
            origin ([int, int, int]): coordinates of center of the rotation
            thickness (float): thickness of the slice
            tissues ([t_id, ]): list of tissue ids to plot
            angle_unit (str): if `'degree'` the angle is treated as degrees.
                Otherwise it is treated a radii
            nb_interp (int): number of pucks to interpolate in between
                existing pucks
            gene_list ([str, ]): list of the gene names to interpolate
                (only selected genes can be inputed)
        Returns:
            points_to_plot (n x 2 ndarray): list of the positions of the points
                that have been plotted
        """
        if tissues is None:
            tissues = self.all_tissues
        if angle_unit == 'degree':
            angle = np.deg2rad(angle)
        x_angle, y_angle, z_angle = angle
        rot_x = tr.rotation_matrix_py(x_angle, [1, 0, 0], origin)
        rot_y = tr.rotation_matrix_py(y_angle, [0, 1, 0], origin)
        rot_z = tr.rotation_matrix_py(z_angle, [0, 0, 1], origin)
        rot_composed = rot_x@rot_y@rot_z
        new_axis = (np.hstack([rot_orig, 1])@rot_composed)[:-1]
        equation = lambda pos: np.sum(new_axis*pos, axis=1)-origin@new_axis
        points, colors, genes = self.produce_em(5, tissues_to_plot=None, gene_list=gene_list)
        points = np.array(points)
        colors = np.array(colors)
        genes = np.array(genes)
        plan = (np.abs(equation(points))<thickness)
        points_to_plot = points[plan]
        points_to_plot = (np.hstack([points_to_plot, [[1]]*points_to_plot.shape[0]])@rot_composed)[:, :-1]
        color_to_plot = colors[plan]
        genes_to_plot = genes.T[plan]
        df = pd.DataFrame(genes_to_plot, columns=gene_list)
        D = anndata.AnnData(df)
        D.obsm['X_Spatial'] = points_to_plot
        D.obs['predicted.id'] = [str(k) for k in color_to_plot]
        D.write(file_name)

        return points_to_plot

    def anndata_no_extra(self, file_name, angle, rot_orig=[0, 0, 1],
                         origin=[0, 0, 0], thickness=30, angle_unit='degree'):
        """
        Build a anndata file containing a slice without doing interpolation
        but any gene can be requested

        Args:
            file_name (str): path to the output ply file
            angle (float): angle of the rotation of the slice
            color_map (matplotlib.cmap): color map that will be applied
            rot_origin ([int, int, int]): 3D vector of the normal of the
                rotation plan. If [0, 0, 1] is given the rotation will be
                around the z axis
            origin ([int, int, int]): coordinates of center of the rotation
            thickness (float): thickness of the slice
            tissues ([t_id, ]): list of tissue ids to plot
            angle_unit (str): if `'degree'` the angle is treated as degrees.
                Otherwise it is treated a radii
            nb_interp (int): number of pucks to interpolate in between
                existing pucks
            gene_list ([str, ]): list of the gene names to interpolate
        Returns:
            points_to_plot (n x 2 ndarray): list of the positions of the points
                that have been plotted
        """
        if angle_unit == 'degree':
            angle = np.deg2rad(angle)
        x_angle, y_angle, z_angle = angle
        rot_x = tr.rotation_matrix_py(x_angle, [1, 0, 0], origin)
        rot_y = tr.rotation_matrix_py(y_angle, [0, 1, 0], origin)
        rot_z = tr.rotation_matrix_py(z_angle, [0, 0, 1], origin)
        rot_composed = rot_x@rot_y@rot_z
        new_axis = (np.hstack([rot_orig, 1])@rot_composed)[:-1]
        equation = lambda pos: np.sum(new_axis*pos, axis=1)-origin@new_axis
        cells = np.array(sorted(self.all_cells))
        pos = np.array([list(self.final[c])+[self.z_pos[c]] for c in cells])
        tissue = np.array([self.tissue[c] for c in cells])
        kept = cells[(np.abs(equation(pos))<thickness)]
        data_tmp = self.anndata.copy()
        data_tmp = data_tmp[kept]
        pos_final = np.array([list(self.final[c])+[self.z_pos[c]] for c in kept])
        pos_final = (np.hstack([pos_final, [[1]]*pos_final.shape[0]])@rot_composed)[:, :-1]
        data_tmp.obsm['X_spatial2'] = pos_final
        data_tmp.write(file_name)


    def produce_em(self, nb_intra=5, tissues_to_plot=None,
                   gene=None, gene_list=None):
        """
        Interpolates beads from the previously computed splines and returns
        the list of the interpolated positions together with a list of values
        for each position corresponding either to the tissue id of the position
        or to the gene expression value if a gene name is provided.

        Args:
            nb_intra (int): number of interpolated slices to add between
                real slices
            tissues_to_plot ([t_id, ]): list of tissue ids to interpolate
                if `None` all tissues are interpolated
            gene (str): name of a gene to output its interpolated value
                for each bead
            gene_list ([str, ]): list of gene names to interpolate. If
                a gene list is given, the list gene_expr is returned.
                The list contains list of gene expressions for the interpolated
                beads. Only pre-selected genes can be inputed
        Returns:
            points (n x 3 ndarray): ndarray containing `n` bead positions
            colors (ndarray of length n): list of bead values. Tissue id
                by default gene expression value if `gene` is not `None`.
            gene_expr (`len(gene_list)` x n ndarray): array of `colors` like
                arrays containing gene expression of the genes queried in
                `gene_list`
        """
        old_spacing = sorted(set(self.z_pos.values()))
        new_spacing = np.linspace(min(old_spacing), max(old_spacing),
                                  len(old_spacing)+(len(old_spacing)-1)*nb_intra)
        points = []
        colors = []
        if gene_list is not None:
            gene_expr = [[] for _ in range(len(gene_list))]
        for c, (min_z, max_z, traj_x, traj_y) in self.all_trajs.items():
            if tissues_to_plot is None or self.tissue[c] in tissues_to_plot:
                spacing = new_spacing[(min_z<=new_spacing)&(new_spacing<=max_z)]
                points.extend(zip(traj_x(spacing), traj_y(spacing), spacing))
                if self.all_expr=={} or gene is None:
                    colors.extend([self.tissue[c]]*len(spacing))
                else:
                    min_z, max_z, traj_expr = self.all_expr[gene][c]
                    colors.extend(traj_expr(spacing))
                if gene_list is not None:
                    for g, L in zip(gene_list, gene_expr):
                        min_z, max_z, traj_expr = self.all_expr[g][c]
                        L.extend(traj_expr(spacing))
        if gene_list is not None:
            return points, colors, gene_expr
        else:
            return points, colors

    def get_differential_genes_main_axis(self, t,
                                         plot=True,
                                         th_diff=1.2,
                                         th_bins=1,
                                         nb_bins=5,
                                         th_expr=.5,
                                         exclusion_func=None):
        """
        Compute some differential gene expression along the main axis
        of a given tissue. Deprecated
        """
        pca = PCA(n_components=3)
        if exclusion_func is not None:
            cells = np.array([c for c in self.all_cells if self.tissue[c]==t if self.final[c][0]<0])
        else:
            cells = np.array([c for c in self.all_cells if self.tissue[c]==t])
        pos = np.array([list(self.final[c]) + [self.z_pos[c]] for c in cells])
        new_pos = pca.fit_transform(pos)
        order = [np.argsort(new_pos[:,0])]
        order.append(np.argsort(new_pos[:,1]))
        order.append(np.argsort(new_pos[:,2]))
        data = self.anndata.copy()
        expressing_genes = th_expr<np.mean(data[cells,:].X, axis=0)
        if plot:
            plot_data = {}
            plot_data['mean'] = np.mean(data[cells,:].X, axis=0)
            plot_data['std'] = np.std(data[cells,:].X, axis=0)
            sns.displot(plot_data, x='mean', y='std', binwidth=(.5, .25), cbar=True)
        data.X[data.X<th_expr]=np.nan
        gene_names = set()
        gene_names_axe = {}
        for axe in range(3):
            pos_compressed = new_pos[:,axe]
            bins = np.linspace(np.percentile(pos_compressed, 1),
                               np.percentile(pos_compressed, 99),
                               nb_bins+1)
            bin_cells = [cells[(bins[i]<pos_compressed) & (pos_compressed<bins[i+1])] for i in range(nb_bins)]
            val_bin = []
            for v in bin_cells:
                mean = np.nanmean(data[v,:].X, axis=0).toarray()
                mean[np.isnan(mean)] = 0
                val_bin.append(mean)
            val_bin = np.array(val_bin).T
            val_bin[np.isnan(val_bin)]=0
            expr_diff = np.array([np.abs(val_bin[:,c1]-val_bin[:,c2])
                                  for c1, c2 in list(combinations(range(nb_bins), 2))])

            candidate_genes = np.where((th_bins<np.sum(th_diff<expr_diff, axis=0))&expressing_genes)[0]
            gene_names.update(data.var_names[candidate_genes])
            gene_names_axe[axe] = data.var_names[candidate_genes]
        return gene_names, gene_names_axe

    def threshold_otsu(self, values, nbins=256):
        """Return threshold value based on Otsu's method.
            Parameters
            ----------
            image : array
            Input image.
            nbins : int
            Number of bins used to calculate histogram. This value is ignored for
            integer arrays.
            Returns
            -------
            threshold : float
            Threshold value.
            References
            ----------
            .. [1] Wikipedia, http://en.wikipedia.org/wiki/Otsu's_Method
        """
        map(np.histogram, values)
        hist, bin_edges = np.histogram(values, nbins)
        bin_centers = (bin_edges[:-1] + bin_edges[1:]) / 2.

        hist = hist.astype(float)

        # class probabilities for all possible thresholds
        weight1 = np.cumsum(hist)
        weight2 = np.cumsum(hist[::-1])[::-1]
        # class means for all possible thresholds
        mean1 = np.cumsum(hist * bin_centers) / weight1
        mean2 = (np.cumsum((hist * bin_centers)[::-1]) / weight2[::-1])[::-1]

        # Clip ends to align class 1 and class 2 variables:
        # The last value of `weight1`/`mean1` should pair with zero values in
        # `weight2`/`mean2`, which do not exist.
        variance12 = weight1[:-1] * weight2[1:] * (mean1[:-1] - mean2[1:]) ** 2

        idx = np.argmax(variance12)
        threshold = bin_centers[:-1][idx]
        return threshold

    def compute_expr_thresholds(self, all_genes=True):
        """
        Compute the expression threshold for all genes

        Returns:
            th ([float, ] ndarray): list of thresholds for each genes
                following the same order as the gene order in `self.anndata`
        """
        if all_genes:
            out = map(self.threshold_otsu, self.anndata.raw.X.todense().T)
        else:
            out = map(self.threshold_otsu, self.anndata.X.T)
        th = []
        for o in out:
            th += [o]
        th = np.array(th)
        return th

    def neighbs(self, gene, sub_data, cells):
        """
        Compute the average number of positive neighbors for the positive cells
        within a given tissue, given a gene

        Args:
            gene (int): gene id (position in the `self.anndata` array)
            sub_data (ndarray): sliced version of `self.anndata` only containing
                the beads corresponding to the tissue to analyse
            cells (ndarray): ids of the cells in `Embryo` ordered similarly to
                the `self.anndata` array (to do correspondancy)
        Returns:
            avg_nb_neighbs (float): average number of positive neighbors per
                positive cells
        """

        # Position of positive cells in `self.anndata`
        positive_cells = np.where(self.gene_expr_th[gene]<sub_data[:,gene])[0]

        # Ids of positive cells
        positive_cells = cells[positive_cells]

        nb_neighbs = []
        for p in positive_cells:
            nb_neighbs.append(len([n for n in self.full_GG[p] if n in positive_cells]))
        avg_nb_neighbs = np.mean(nb_neighbs)
        return avg_nb_neighbs

    def cell_groups(self, t, th_vol=.025,
                    all_genes=True):
        """
        Compute the local expression metric for each gene in a given tissue `t`

        Args:
            t (int): tissue id to process
            th_vol (float 0<th_vol<1): high and low volume threshold.
                Any gene expression that covers more that 1-th_vol
                fraction of the tissue volume or less that th_vol fraction
                of the tissue volume is discarded.
            all_genes (bool): True if all the genes should be considered.
                Otherwise only the previously computed variable genes are
                considered
        Returns:
            data_plot (pandas.DataFrame): pandas DataFrame containing most
                of the computed information for gene localization of the tissue
                `t`. The main value is in the column `Distance_to_reg`
        """
        if all_genes:
            data = self.anndata.raw.X
        else:
            data = self.anndata.copy().X
        cells = np.array([c for c in self.all_cells if self.tissue[c]==t])

        # Spliting the array to only have tissue *t* cells
        sub_data = data[cells]
        if all_genes:
            sub_data = np.array(sub_data.todense())

        # Occupied volume for the cells of tissue *t*
        volume_total = len(cells)

        # Volume ratio for cell expressing in a tissue for each gene
        sub_volumes = np.sum(self.gene_expr_th<sub_data, axis=0) / volume_total

        # Mask for the list of genes that are expressing enough within the tissue
        mask_expr = (th_vol<sub_volumes)&(sub_volumes<1-th_vol)

        # List of genes that are expressing enough within the tissue
        interesting_genes = np.where(mask_expr)[0]

        # Copmuting the number of cells expressing
        avg_nb_neighbs = []
        for g in interesting_genes:
            nb_N_for_g = self.neighbs(g, sub_data, cells)
            avg_nb_neighbs.append(nb_N_for_g / self.whole_tissue_nb_N[t])
        avg_nb_neighbs = np.array(avg_nb_neighbs)

        # Build a dataframe with the previously computed metrics
        data_plot = {
            'volume': sub_volumes[mask_expr],
            'avg_nb_neighbs': avg_nb_neighbs,
        }

        # Compute the linear regression
        # Value against which the linear regression is done
        # It is important that the relationship between x and y is linear!!!
        regression_x = 'avg_nb_neighbs'
        regression_y = 'volume'

        regr = linear_model.LinearRegression()
        data_x_reshaped = data_plot[regression_x].reshape(-1,1)
        data_y_reshaped = data_plot[regression_y].reshape(-1,1)
        regr.fit(data_x_reshaped, data_y_reshaped)
        b = regr.intercept_[0]
        a = regr.coef_[0][0]
        line = lambda x: (a*np.array(x)+b)

        data_plot['Distance_to_reg'] = np.abs((data_y_reshaped-
                                               regr.predict(data_x_reshaped))[:,0])
        data_plot['Interesting genes'] = interesting_genes
        if all_genes:
            data_plot['Gene names'] = np.array(self.anndata.raw[:,data_plot['Interesting genes']].var_names)
        else:
            data_plot['Gene names'] = np.array(self.anndata[:,data_plot['Interesting genes']].var_names)
        data_plot = pd.DataFrame(data_plot)

        return data_plot

    def get_3D_differential_expression(self, tissues_to_process, th_vol=.025,
                                       all_genes=True):
        """
        Compute the 3D spatial differential expression for a list of tissues and
        stores it in `self.diff_expressed_3D`.

        Args:
            tissues_to_process ([t_ids, ]): list of tissue ids to process
            th_vol (float 0<th_vol<1): high and low volume threshold.
                Any gene expression that covers more that 1-th_vol
                fraction of the tissue volume or less that th_vol fraction
                of the tissue volume is discarded.
            all_genes (bool): True if all the genes should be considered.
                Otherwise only the previously computed variable genes are
                considered
        Returns:
            self.diff_expressed_3D (dict t_id: pandas.DataFrame):
                dictionary that maps a tissue to a pandas DataFrame containing
                most of the computed information for gene localization of
                the tissue `t_id`. The main value is in the column `Distance_to_reg`
        """
        cells = list(self.all_cells)
        pos_3D = [np.array(list(self.final[c])+[self.z_pos[c]]) for c in cells]

        if not hasattr(self, 'full_GG'):
            self.full_GG = self.build_gabriel_graph(cells, pos_3D)

        if not hasattr(self, 'gene_expr_th'):
            self.gene_expr_th = self.compute_expr_thresholds(all_genes=all_genes)

        if not hasattr(self, 'whole_tissue_nb_N'):
            self.whole_tissue_nb_N = {}
            for t in self.all_tissues:
                cells = np.array([c for c in self.all_cells if self.tissue[c]==t])
                self.whole_tissue_nb_N[t] = np.median([len(self.full_GG[c]) for c in cells])

        if not hasattr(self, 'diff_expressed_3D'):
            self.diff_expressed_3D = {}
            
        for t in tissues_to_process:
            if not t in self.diff_expressed_3D:
                self.diff_expressed_3D[t] = self.cell_groups(t, th_vol=th_vol,
                                                             all_genes=all_genes)

        if hasattr(self, 'tissues_diff_expre_processed'):
            self.tissues_diff_expre_processed.extend(tissues_to_process)
        else:
            self.tissues_diff_expre_processed = tissues_to_process
        self.all_genes = all_genes
        return self.diff_expressed_3D

    def plot_top_3D_diff_expr_genes(self, tissues_to_process, nb_genes=20,
                                   repetition_allowed=False, compute_z_score=True,
                                   fig=None, ax=None, output_path=None):
        """
        Plot the top `nb_genes` genes for 3D differentially expressed genes for a
        list of tissues.

        Args:
            tissues_to_process ([t_ids, ]): list of tissue ids to process
            nb_genes (int): number of genes in the top gene list
            repetition_allowed (bool): if true, a gene can be in the top
                `nb_genes` of multiple tissues. Otherwise it is only kept
                for the tissue it has the highest localization score.
                Default: False
            compute_z_score (bool): if true, the z-score of gene expression is computed
                for each gene independently, otherwise the initial value from `self.anndata`
                is kept
                Default: True
            fig (matplotlib.Figure): figure onto which ploting the output. If fig
                is given ax should be given too. If None, a new figure is created
                or recovered from ax.
                Default: None
            ax (matplotlib.AxesSubplot): the working figure axis
                Default: None
            output_path (str): path to the desired output figure. If None, the figure
                is not saved
                Default: None
        """
        tmp_T = set(tissues_to_process).difference(self.tissues_diff_expre_processed)
        if len(tmp_T) != 0:
            print("You asked to plot tissue(s) that were not already processed")
            print("The following tissue(s) will be ignored:")
            for t in tmp_T:
                print(f"\t - id: {t}, name: {self.corres_tissue[t]}")
        tissues_to_process = list(set(tissues_to_process).intersection(self.tissues_diff_expre_processed))
        genes_of_interest = []
        gene_dict = {}
        tissue_genes = {}
        genes_in = {}
        added_genes = 1 if repetition_allowed else 4
        for t in tissues_to_process:
            data_t = self.diff_expressed_3D[t]
            G_N = data_t.sort_values('Distance_to_reg')['Interesting genes'][:-nb_genes*added_genes-1:-1]
            G_V = data_t.sort_values('Distance_to_reg')['Distance_to_reg'][:-nb_genes*added_genes-1:-1]
            genes_of_interest.extend(G_N[:nb_genes])
            for g, v in zip(G_N, G_V):
                tissue_genes.setdefault(g, []).append(t)
                gene_dict[(t, g)] = v
            genes_in[t] = list(G_N)

        if not repetition_allowed:
            dict_counter = Counter(genes_of_interest)
            acc = 0
            while any([1<k for k in dict_counter.values()]):
                t = tissues_to_process[acc%len(tissues_to_process)]
                for g in genes_in[t]:
                    if 1<dict_counter[g]:
                        tissues = np.array(tissue_genes[g])
                        values = [gene_dict[(t, g)] for t in tissues]
                        if tissues[np.argsort(values)][-1]!=t:
                            genes_in[t].remove(g)
                genes_of_interest = []
                for t in tissues_to_process:
                    genes_of_interest.extend(genes_in[t][:nb_genes])
                dict_counter = Counter(genes_of_interest)
                acc += 1
        values = np.zeros((nb_genes*len(tissues_to_process), len(tissues_to_process)))
        tissue_order = []
        for i, g in enumerate(genes_of_interest):
            for j, t in enumerate(tissues_to_process):
                data_t = self.diff_expressed_3D[t]
                if g in data_t['Interesting genes'].values:
                    values[i, j] = data_t[data_t['Interesting genes']==g]['Distance_to_reg']
                if i==0:
                    tissue_order.append(t)
        # z_score = (values - np.mean(values, axis=1).reshape(-1, 1))/np.std(values, axis=1).reshape(-1, 1)
        if compute_z_score:
            values = stats.zscore(values, axis=0)
        if ax is None:
            fig, ax = plt.subplots(figsize=(5,max(5, round(1.5*nb_genes))))
        if fig is None:
            fig = ax.get_figure()
        ax.imshow(values, interpolation='nearest', cmap='Reds')
        ax.set_xticks(range(len(tissue_order)))
        ax.set_xticklabels([self.corres_tissue[t] for t in tissue_order], rotation=90)
        ax.set_yticks(range(values.shape[0]))
        if self.all_genes:
            ax.set_yticklabels(list(self.anndata.raw[:,genes_of_interest].var_names))
        else:
            ax.set_yticklabels(list(self.anndata[:,genes_of_interest].var_names))
        fig.tight_layout()
        if output_path is not None:
            fig.savefig(output_path)
        return fig, ax

    def plot_volume_vs_neighbs(self, t, print_top=None,
                               print_genes=None, fig=None, ax=None,
                               output_path=None, **kwargs):
        """
        Plot volume of expressing cells versus the average number of expressing neighbors
        for a given tissue.

        Args:
            t (int): tissue id to treat
            print_top (int): number of gene names to plot onto the figure (slows down)
                the function significantly
            fig (matplotlib.Figure): figure onto which ploting the output. If fig
                is given ax should be given too. If None, a new figure is created
                or recovered from ax.
                Default: None
            ax (matplotlib.AxesSubplot): the working figure axis
                Default: None
            output_path (str): path to the desired output figure. If None, the figure
                is not saved
                Default: None
            kwargs are forwarded to the seaborn.scatterplot
        """
        data_plot = self.diff_expressed_3D[t]
        if ax is None:
            fig, ax = plt.subplots(figsize=(10, 8))
        if fig is None:
            fig = ax.get_figure()
        x = 'avg_nb_neighbs'
        y = 'volume'
        print_text = True
        g = sns.scatterplot(data=data_plot, x=x, y=y, ax=ax, hue='Distance_to_reg', **kwargs)
        legend = g.axes.get_legend()
        legend.set_title('Localization score')
        ax.set_ylabel('Relative volume (to total tissue volume)')
        ax.set_xlabel('Relative cell density (to the average cell density within the tissue)')
        if print_top is not None:
            top_X = data_plot.sort_values('Distance_to_reg', ascending=False)[:print_top]
            x_values = top_X[x]
            y_values = top_X[y]
            names = top_X['Gene names']
            for name, x, y in zip(names, x_values, y_values):
                plt.text(x=x,y=y,s=name,
                         fontdict=dict(color='red',size=8, fontweight='bold'), va='baseline')
        if print_genes is not None:
            for gene in print_genes:
                gene_num_all = np.where(self.anndata.var_names==gene)[0][0]
                gene_num = np.where(data_plot['Interesting genes']==gene_num_all)[0][0]
                plt.text(x=data_plot[x][gene_num],y=data_plot[y][gene_num],s=gene,
                         fontdict=dict(color='red',size=8, fontweight='bold'), va='baseline')
        ax.set_title(self.corres_tissue[t])
        fig.tight_layout()
        if output_path is not None:
            fig.savefig(output_path)

    def print_diff_expr_genes(self, tissue, nb):
        """
        Extract from the DataFrame `self.diff_expressed_3D`
        the genes that are locally expressed for a given tissue

        Args:
            tissue (int): id of the tissue to look at
            nb (int): number of genes to extract
        Returns
            order (`nb` x m pandas.DataFrame): DataFrame containing
                the top `nb` localized genes.
        """
        data_plot = self.diff_expressed_3D[tissue]
        gene_values = []
        order = data_plot.sort_values('Distance_to_reg', ascending=False)[:nb]
        return order

    def __init__(self, data_path, tissues_to_ignore=None,
                 corres_tissue=None, tissue_weight=None,
                 xy_resolution=1, genes_of_interest=None,
                 nb_CS_begin_ignore=0, nb_CS_end_ignore=0,
                 store_anndata=False):
        """
        Initialize an spatial single cell embryo

        Args:
            data_path (str): path to the file containing the sc data (h5ad format)
            tissues_to_ignore ([t_ids, ]): list of tissues to ignore. Beads belonging
                to these tissues will be discarded
            corres_tissue ({t_id: str}): dictionary that maps a tissue id to a tissue
                name
            tissue_weight ({t_id: int}): dictionary that maps a tissue id to a weight
                that will be used for the puck registration. The higher the value is
                the more aligned the tissue will be. The default value is 1
            xy_resolution (float): resolution in x and y (assumed to be isotrope)
            gene_of_interest ([str, ]): list of gene names to be selected. For some
                applications, they will be the only ones that can be processed.
            nb_CS_begin_ignore (int): number of pucks to ignore at the begining of
                the stack.
                Default: 0
            nb_CS_end_ignore (int): number of pucks to ignore at the end of the stack
                Default: 0
            store_anndata (bool): if true the anndata array is stored. Necessary when
                doing 3D differential expression analysis
                Default: False
        """
        self.cells = set()
        self.pos = {}
        self.cover_slip = {}
        self.tissue = {}
        self.all_tissues = set()
        self.cells_from_cover_slip = {}
        self.cells_from_tissue = {}
        self.all_cover_slips = []
        self.nb_CS_begin_ignore=nb_CS_begin_ignore
        self.nb_CS_end_ignore=nb_CS_end_ignore
        self.tissues_to_ignore = [] if tissues_to_ignore is None else tissues_to_ignore
        self.corres_tissue = {} if corres_tissue is None else corres_tissue
        self.tissue_weight = {} if tissue_weight is None else tissue_weight
        if data_path.split('.')[-1] == 'h5ad':
            self.read_anndata(data_path, xy_resolution=xy_resolution,
                              genes_of_interest=genes_of_interest,
                              store_anndata=store_anndata)

    