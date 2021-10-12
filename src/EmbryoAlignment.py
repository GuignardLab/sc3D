import numpy as np
import transformations as tr
from scipy.spatial import KDTree, Delaunay
from scipy.spatial.distance import cdist
from scipy.optimize import linear_sum_assignment
from scipy.interpolate import InterpolatedUnivariateSpline, interp1d
from itertools import combinations
from matplotlib import pyplot as plt
import open3d as o3d
import pandas as pd
import anndata

class Embryo(object):
    def set_zpos(self, z_space=30.):
        self.z_space = z_space
        self.z_pos = {}
        cs_conversion = {b: a*z_space for a, b in enumerate(self.all_cover_slips)}
        for c in self.all_cells:
            self.z_pos[c] = cs_conversion[self.cover_slip[c]]

    def read_csv(self, path, xy_resolution=1):
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

        # sanity check
        #if linalg.matrix_rank(H) < 3:
        #    raise ValueError("rank of H = {}, expecting 3".format(linalg.matrix_rank(H)))

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
        self.centered_pos = {}
        for cs, cells in self.cells_from_cover_slip.items():
            pos = np.array([self.pos[c] for c in cells])
            avg = np.mean(pos, axis=0)
            self.centered_pos.update(zip(cells, pos-avg))
        return(self.centered_pos)
    
    def get_tissue_centers(self):
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
                    # hull = ConvexHull(pos)
                    # self.tissue_centers[cs][tissue] = np.mean(hull.points[hull.vertices], axis=0)
                    # self.tissue_centers[cs][tissue] = np.median(pos, axis=0)
                    for w in range(self.tissue_weight.get(tissue, 1)):
                        self.tissue_centers[cs][(tissue, w)] = np.mean(pos, axis=0)
        return(self.tissue_centers)
    
    def build_and_apply_trsf_matrix(self, cs_ref, cs_flo):
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
        if not hasattr(self, 'registered_pos'):
            self.register_with_tissues()
        if not hasattr(self, 'pos_reg_aff'):
            self.pos_reg_aff = {}
        if not hasattr(self, 'pairing'):
            self.pairing = {}
        if not hasattr(self, 'final') and final:
            self.final = {c: self.centered_pos[c] for c in self.cells_from_cover_slip[cs1]}
        # pos_ref, pos_flo = self.build_pairing(cs1, cs2, rebuild=True, refine=refine, th_d=th_d)
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

    def transfer_cells(self, cs1, cs2):
        if not hasattr(self, 'final'):
            self.register_cs(cs1, cs2, rigid=True, final=True)
        not_paired = self.cells_from_cover_slip[cs1].difference(self.pairing)
        
    def build_gabriel_graph(self, node_ids, pos):
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

    def get_KDTree_GG(self, cs, source):
        if not hasattr(self, 'cells_paired_cs'):
            self.cells_paired_cs = {}
        if not hasattr(self, 'cells_not_paired_cs'):
            self.cells_not_paired_cs = {}
        if not hasattr(self, 'pos_paired_cs'):
            self.pos_paired_cs = {}
        if not hasattr(self, 'pos_not_paired_cs'):
            self.pos_not_paired_cs = {}
        if not hasattr(self, 'kdt_cs'):
            self.kdt_cs = {}
        if not hasattr(self, 'GG_cs'):
            self.GG_cs = {}
        if source==cs:
            cells_not_paired_cs = self.cells_from_cover_slip[cs].difference(self.pairing.keys())
        else:
            cells_not_paired_cs = self.cells_from_cover_slip[cs].difference(self.pairing.values())        
        cells_paired_cs = list(self.cells_from_cover_slip[cs].difference(cells_not_paired_cs))
        cells_not_paired_cs = list(cells_not_paired_cs)
        pos_not_paired_cs = [self.final[c] for c in cells_not_paired_cs]
        pos_paired_cs = np.array([self.final[c] for c in cells_paired_cs])
        kdt_cs = KDTree(pos_paired_cs)
        all_cells = list(self.cells_from_cover_slip[cs])
        GG_cs = self.build_gabriel_graph(all_cells, [self.final[c] for c in all_cells])
        self.cells_paired_cs[cs] = cells_paired_cs
        self.cells_not_paired_cs[cs] = cells_not_paired_cs
        self.pos_paired_cs[cs] = pos_paired_cs
        self.pos_not_paired_cs[cs] = pos_not_paired_cs
        self.kdt_cs[cs] = kdt_cs
        self.GG_cs[cs] = GG_cs


    def build_trajectories(self, cs, source, disapear_bounds=None):
        if not hasattr(self, 'trajectories'):
            self.trajectories = {}
        if disapear_bounds is None:
            disapear_bounds = (.1, .5, .9)
        
        nb_genes = len(self.all_genes)
        trajectories = {}
        gene_traj = {}
        if source==cs:
            cells = list(self.cells_paired_cs[cs])
            pairs = [(c, self.pairing[c]) for c in cells]
            to_map = lambda pair: list(map(interp1d, [[0, 1],]*nb_genes, zip(*self.data[pair, genes].X)))
            gene_traj = dict(zip(cells, map(to_map, pairs)))
            for c in self.cells_paired_cs[cs]:
                start_x, start_y = self.final[c]
                end_x, end_y = self.final[self.pairing[c]]
                interp_x = interp1d([0, 1], [start_x, end_x])
                interp_y = interp1d([0, 1], [start_y, end_y])
                trajectories[c] = interp_x, interp_y
        start, mid, end = disapear_bounds
        d_to_closest = {}
        arrival = {}
        gene_arrival = {}
        if source==cs:
            pairing_to_use = self.pairing
        else:
            pairing_to_use = {v: k for k, v in self.pairing.items()}
        for c in self.cells_not_paired_cs[cs]:
            neighbs = self.GG_cs[cs].get(c).difference(self.cells_not_paired_cs[cs])
            if len(neighbs)<1:
                neighbs = [self.cells_paired_cs[cs][self.kdt_cs[cs].query(self.final[c], 1)[1]]]
            arrival[c] = np.mean([self.final[pairing_to_use[ni]] for ni in neighbs], axis=0)
            d_to_closest[c] = np.mean([np.linalg.norm(self.final[c] - self.final[ni]) for ni in neighbs])
        d_to_closest_vals = list(d_to_closest.values())
        med_to_closest = np.median(d_to_closest_vals)
        min_to_closest = np.min(d_to_closest_vals)
        max_to_closest = np.max(d_to_closest_vals)
        if source==cs:
            dist_to_disapear = interp1d([max_to_closest, med_to_closest, min_to_closest],
                                        [start, mid, end])
        else:
            dist_to_disapear = interp1d([min_to_closest, med_to_closest, max_to_closest],
                                        [start, mid, end])
        for c in self.cells_not_paired_cs[cs]:
            start_x, start_y = self.final[c]
            end_x, end_y = arrival[c]
            D = dist_to_disapear(d_to_closest[c])
            if source==cs:
                interp_x = interp1d([0, D], [start_x, end_x], bounds_error=False)
                interp_y = interp1d([0, D], [start_y, end_y], bounds_error=False)
            else:
                interp_x = interp1d([1, D], [start_x, end_x], bounds_error=False)
                interp_y = interp1d([1, D], [start_y, end_y], bounds_error=False)
            trajectories[c] = interp_x, interp_y
        if source in self.trajectories:
            self.trajectories[source].update(trajectories)
        else:
            self.trajectories[source] = trajectories
        
    def get_new_positions(self, P, cs):
        pos_tmp = {c: np.array([self.trajectories[cs][c][0](P), self.trajectories[cs][c][1](P)])
                   for c in self.trajectories[cs]}
        cells = [c for c in pos_tmp if not np.isnan(pos_tmp[c][0])]
        pos = np.array([pos_tmp[c] for c in cells])
        return cells, pos
            
    # def reconstruct_intermediate(self, rigid=True, th_d=True, cs=None):
    #     if cs is not None:
    #         cs_to_treat = cs
    #     else:
    #         cs_to_treat = self.all_cover_slips
    #     for i, cs1 in enumerate(cs_to_treat[:-1]):
    #         cs2 = cs_to_treat[i+1]
    #         M1 = self.register_cs(cs1, cs2, rigid=rigid, final=True, th_d=th_d)
    #         self.get_KDTree_GG(cs1, source=cs1)
    #         self.get_KDTree_GG(cs2, source=cs1)
    #         # self.build_trajectories(cs1, source=cs1)
    #         # self.build_trajectories(cs2, source=cs1)
        
    def plot_coverslip(self, cs, pos='pos', ax=None,
                       tissues_to_plot=None, legend=False,
                       color=None, **kwargs):
        if ax is None:
            fig, ax = plt.subplots()
        else:
            fig = ax.get_figure()
        if isinstance(pos, str):
            positions_attr = self.__getattribute__(pos)
        else:
            positions_attr = pos
        if tissues_to_plot is None:
            positions = np.array([positions_attr[c] for c in self.cells_from_cover_slip[cs]])
            tissues = [self.tissue[c] for c in self.cells_from_cover_slip[cs]]
        else:
            cells = [c for c in self.cells_from_cover_slip[cs] if self.tissue[c] in tissues_to_plot]
            positions = np.array([positions_attr[c] for c in cells])
            tissues = [self.tissue[c] for c in cells]
        if len(positions)<1:
            return fig, ax
        scatter_args = {'marker':'.', 's':25, 'cmap':'tab20',
                        'vmin':min(self.all_tissues), 'vmax':max(self.all_tissues)}
        scatter_args.update(kwargs)
        if color == None:
            color = tissues
        scatter = ax.scatter(*positions.T, c=color, **scatter_args)
        if legend:
            ax.legend(handles=scatter.legend_elements()[0], labels=np.unique(tissues))
        return fig, ax

    def removing_spatial_outliers(self, th=.2):
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
            dpgmm = mixture.GaussianMixture(n_components=3, max_iter=1000,
                                            covariance_type='full').fit(D.reshape(-1,1))
            proba0 = dpgmm.predict_proba(D.reshape(-1, 1))[:, 0]
            proba1 = dpgmm.predict_proba(D.reshape(-1, 1))[:, 1]
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

        start_points = {}
        cells_to_treat = set(self.all_cells)
        all_trajs = []
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

            all_trajs.append([traj[0], min(z_traj), max(z_traj), f_traj_x, f_traj_y])
            cells_to_treat -= set(traj)
        self.all_trajs = all_trajs
        self.all_expr = all_expr

    def plot_slice(self, angle, color_map, rot_orig=[0, 0, 1], origin=[0, 0, 0],
                   thickness=30, tissues=None, angle_unit='degree',
                   nb_interp=5, output_path=None, **kwargs):
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
        points, color = self.produce_em(nb_interp, tissues)
        points = np.array(points)
        color = np.array(color)
        dist_to_plan = equation(points)
        plan = (np.abs(dist_to_plan)<thickness)
        dist_to_plan = dist_to_plan[plan]
        points_to_plot = points[plan]
        points_to_plot = (np.hstack([points_to_plot, [[1]]*points_to_plot.shape[0]])@rot_composed)[:, :-1]
        color_to_plot = np.array([color_map[c] for c in color[plan]])
        p_order = np.argsort(dist_to_plan)
        points_to_plot = points_to_plot[p_order]
        color_to_plot = color_to_plot[p_order]
        fig = plt.figure(figsize=(8, 8))
        ax = fig.add_subplot(1, 1, 1)#, projection='3d')
        kwargs_scatter = { 's':5, 'color':color_to_plot}
        kwargs_scatter.update(kwargs)
        ax.scatter(*(points_to_plot.T[:-1]), **kwargs_scatter)
        ax.axis('equal')
        if output_path is not None:
            fig.savefig(output_path)
        return points_to_plot

    def ply_slice(self, file_name, angle, color_map, rot_orig=[0, 0, 1],
                  origin=[0, 0, 0], thickness=30, tissues=None,
                  tissues_colored=None, angle_unit='degree',
                  nb_interp=5):
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
        points, color = self.produce_em(nb_interp, tissues)
        points = np.array(points)
        color = np.array(color)
        plan = (np.abs(equation(points))<thickness)
        points_to_plot = points[plan]
        color_to_plot = color[plan]
        mapping = np.array([[.8, .8, .8] for t in range(max(color_map)+1)])
        for t in tissues_colored:
            mapping[t] = color_map.get(t, [.8, .8, .8])
        pcd = o3d.geometry.PointCloud()
        pcd.points = o3d.utility.Vector3dVector(points_to_plot)
        pcd.colors = o3d.utility.Vector3dVector(mapping[color_to_plot])
        o3d.io.write_point_cloud(file_name, pcd)
        
        return points_to_plot

    def anndata_slice(self, file_name, angle, gene_list, rot_orig=[0, 0, 1],
                      origin=[0, 0, 0], thickness=30, tissues=None,
                      angle_unit='degree', nb_interp=5):
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
        old_spacing = sorted(set(self.z_pos.values()))
        new_spacing = np.linspace(min(old_spacing), max(old_spacing),
                                  len(old_spacing)+(len(old_spacing)-1)*nb_intra)
        points = []
        colors = []
        if gene_list is not None:
            gene_expr = [[] for _ in range(len(gene_list))]
        for c, min_z, max_z, traj_x, traj_y in self.all_trajs:
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

    def __init__(self, data_path, tissues_to_ignore=None,
                 corres_tissue=None, tissue_weight=None,
                 xy_resolution=1, genes_of_interest=None,
                 nb_CS_begin_ignore=0, nb_CS_end_ignore=0,
                 store_anndata=False):
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

    