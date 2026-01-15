from sc3D import SpatialOmicArray
import numpy as np

em = SpatialOmicArray("data/data_test.h5ad", store_anndata=True)

em2 = SpatialOmicArray(
    "data/DLPFC.h5ad",
    tissue_id="layer_guess",
    pos_id="spatial",
    array_id="z",
    z_space=30,
    store_anndata=True,
)


def test_sc3D():
    assert len(em.all_cells) == 120


def test_smooth():
    em.smooth_data()


def test_plot_coverslip():
    em.plot_coverslip(7)


def test_3D_diff():
    em.get_3D_differential_expression([21])


def test_plot_diff():
    em.plot_top_3D_diff_expr_genes([21, 23], repetition_allowed=True)
    em.plot_top_3D_diff_expr_genes([21, 23], repetition_allowed=False)


def test_vol_neighbs():
    em.plot_volume_vs_neighbs(21)


def test_spatial_outliers():
    em.removing_spatial_outliers()


def test_volumes():
    em.compute_volumes()


def test_z_pos():
    em.set_zpos()


def test_produce():
    em2.produce_em()


def test_registration():
    em2.registration_3d()


def test_plot_slices():
    origin = np.mean([em2.final[c] for c in em2.all_cells], axis=0)
    origin = np.hstack([origin, 80])
    angles = np.array([-5.0, 5.0, 0.0])
    _ = em2.plot_slice(
        angles, color_map="viridis", origin=origin, thickness=30, nb_interp=5
    )
