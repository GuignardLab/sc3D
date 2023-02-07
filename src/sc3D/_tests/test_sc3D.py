from sc3D import Embryo
import numpy as np


def test_sc3D():
    em = Embryo("data/data_test.h5ad", store_anndata=True)
    assert len(em.all_cells) == 120
    em.smooth_data()
    em.plot_coverslip(7)
    em.get_3D_differential_expression([21])
    em.plot_top_3D_diff_expr_genes([21, 23], repetition_allowed=True)
    em.plot_top_3D_diff_expr_genes([21, 23], repetition_allowed=False)
    em.plot_volume_vs_neighbs(21)
    em.removing_spatial_outliers()
    em.compute_volumes()
    em.set_zpos()

    em = Embryo(
        f"data/DLPFC.h5ad",
        tissue_id="layer_guess",
        pos_id="spatial",
        array_id="z",
        z_space=30,
        store_anndata=True,
    )
    em.produce_em()
    em.registration_3d()
    origin = np.mean([em.final[c] for c in em.all_cells], axis=0)
    origin = np.hstack([origin, 80])
    angles = np.array([-5.0, 5.0, 0.0])
    points_to_plot = em.plot_slice(
        angles, color_map="viridis", origin=origin, thickness=30, nb_interp=5
    )
