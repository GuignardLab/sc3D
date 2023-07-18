"""
This file is a copy of the useful function
from Christophe Gohlke's amazing package transformations
that can be found there: https://pypi.org/project/transformations/

For a reason that I don't really understand, the package
does not install on newer version of Python therefore
I decided to, for the time being at least, copy some of its work here.

Léo
"""

import numpy as np
import math


class transformations:
    _EPS = np.finfo(float).eps * 4.0

    @classmethod
    def unit_vector(clf, data, axis=None, out=None):
        """Return ndarray normalized by length, i.e. Euclidean norm, along axis.

        >>> v0 = numpy.random.random(3)
        >>> v1 = unit_vector(v0)
        >>> numpy.allclose(v1, v0 / numpy.linalg.norm(v0))
        True
        >>> v0 = numpy.random.rand(5, 4, 3)
        >>> v1 = unit_vector(v0, axis=-1)
        >>> v2 = v0 / numpy.expand_dims(numpy.sqrt(numpy.sum(v0*v0, axis=2)), 2)
        >>> numpy.allclose(v1, v2)
        True
        >>> v1 = unit_vector(v0, axis=1)
        >>> v2 = v0 / numpy.expand_dims(numpy.sqrt(numpy.sum(v0*v0, axis=1)), 1)
        >>> numpy.allclose(v1, v2)
        True
        >>> v1 = numpy.empty((5, 4, 3))
        >>> unit_vector(v0, axis=1, out=v1)
        >>> numpy.allclose(v1, v2)
        True
        >>> list(unit_vector([]))
        []
        >>> list(unit_vector([1]))
        [1.0]

        """
        if out is None:
            data = np.array(data, dtype=np.float64, copy=True)
            if data.ndim == 1:
                data /= math.sqrt(np.dot(data, data))
                return data
        else:
            if out is not data:
                out[:] = np.array(data, copy=False)
            data = out
        length = np.atleast_1d(np.sum(data * data, axis))
        np.sqrt(length, length)
        if axis is not None:
            length = np.expand_dims(length, axis)
        data /= length
        if out is None:
            return data
        return None

    @classmethod
    def rotation_matrix(clf, angle, direction, point=None):
        """Return matrix to rotate about axis defined by point and direction.

        >>> R = rotation_matrix(math.pi/2, [0, 0, 1], [1, 0, 0])
        >>> numpy.allclose(numpy.dot(R, [0, 0, 0, 1]), [1, -1, 0, 1])
        True
        >>> angle = (random.random() - 0.5) * (2*math.pi)
        >>> direc = numpy.random.random(3) - 0.5
        >>> point = numpy.random.random(3) - 0.5
        >>> R0 = rotation_matrix(angle, direc, point)
        >>> R1 = rotation_matrix(angle-2*math.pi, direc, point)
        >>> is_same_transform(R0, R1)
        True
        >>> R0 = rotation_matrix(angle, direc, point)
        >>> R1 = rotation_matrix(-angle, -direc, point)
        >>> is_same_transform(R0, R1)
        True
        >>> I = numpy.identity(4, numpy.float64)
        >>> numpy.allclose(I, rotation_matrix(math.pi*2, direc))
        True
        >>> numpy.allclose(2, numpy.trace(rotation_matrix(math.pi/2,
        ...                                               direc, point)))
        True

        """
        import math

        sina = math.sin(angle)
        cosa = math.cos(angle)
        direction = clf.unit_vector(direction[:3])
        # rotation matrix around unit vector
        R = np.diag([cosa, cosa, cosa])
        R += np.outer(direction, direction) * (1.0 - cosa)
        direction *= sina
        R += np.array(
            [
                [0.0, -direction[2], direction[1]],
                [direction[2], 0.0, -direction[0]],
                [-direction[1], direction[0], 0.0],
            ]
        )
        M = np.identity(4)
        M[:3, :3] = R
        if point is not None:
            # rotation not around origin
            point = np.array(point[:3], dtype=np.float64, copy=False)
            M[:3, 3] = point - np.dot(R, point)
        return M

    @classmethod
    def vector_norm(clf, data, axis=None, out=None):
        """Return length, i.e. Euclidean norm, of ndarray along axis.

        >>> v = numpy.random.random(3)
        >>> n = vector_norm(v)
        >>> numpy.allclose(n, numpy.linalg.norm(v))
        True
        >>> v = numpy.random.rand(6, 5, 3)
        >>> n = vector_norm(v, axis=-1)
        >>> numpy.allclose(n, numpy.sqrt(numpy.sum(v*v, axis=2)))
        True
        >>> n = vector_norm(v, axis=1)
        >>> numpy.allclose(n, numpy.sqrt(numpy.sum(v*v, axis=1)))
        True
        >>> v = numpy.random.rand(5, 4, 3)
        >>> n = numpy.empty((5, 3))
        >>> vector_norm(v, axis=1, out=n)
        >>> numpy.allclose(n, numpy.sqrt(numpy.sum(v*v, axis=1)))
        True
        >>> vector_norm([])
        0.0
        >>> vector_norm([1])
        1.0

        """
        data = np.array(data, dtype=np.float64, copy=True)
        if out is None:
            if data.ndim == 1:
                return math.sqrt(np.dot(data, data))
            data *= data
            out = np.atleast_1d(np.sum(data, axis=axis))
            np.sqrt(out, out)
            return out
        data *= data
        np.sum(data, axis=axis, out=out)
        np.sqrt(out, out)
        return None

    @classmethod
    def quaternion_matrix(clf, quaternion):
        """Return homogeneous rotation matrix from quaternion.

        >>> M = quaternion_matrix([0.99810947, 0.06146124, 0, 0])
        >>> numpy.allclose(M, rotation_matrix(0.123, [1, 0, 0]))
        True
        >>> M = quaternion_matrix([1, 0, 0, 0])
        >>> numpy.allclose(M, numpy.identity(4))
        True
        >>> M = quaternion_matrix([0, 1, 0, 0])
        >>> numpy.allclose(M, numpy.diag([1, -1, -1, 1]))
        True

        """
        q = np.array(quaternion, dtype=np.float64, copy=True)
        n = np.dot(q, q)
        if n < clf._EPS:
            return np.identity(4)
        q *= math.sqrt(2.0 / n)
        q = np.outer(q, q)
        return np.array(
            [
                [
                    1.0 - q[2, 2] - q[3, 3],
                    q[1, 2] - q[3, 0],
                    q[1, 3] + q[2, 0],
                    0.0,
                ],
                [
                    q[1, 2] + q[3, 0],
                    1.0 - q[1, 1] - q[3, 3],
                    q[2, 3] - q[1, 0],
                    0.0,
                ],
                [
                    q[1, 3] - q[2, 0],
                    q[2, 3] + q[1, 0],
                    1.0 - q[1, 1] - q[2, 2],
                    0.0,
                ],
                [0.0, 0.0, 0.0, 1.0],
            ]
        )

    @classmethod
    def affine_matrix_from_points(
        clf, v0, v1, shear=True, scale=True, usesvd=True
    ):
        """Return affine transform matrix to register two point sets.

        v0 and v1 are shape (ndims, -1) arrays of at least ndims non-homogeneous
        coordinates, where ndims is the dimensionality of the coordinate space.

        If shear is False, a similarity transformation matrix is returned.
        If also scale is False, a rigid/Euclidean transformation matrix
        is returned.

        By default the algorithm by Hartley and Zissermann [15] is used.
        If usesvd is True, similarity and Euclidean transformation matrices
        are calculated by minimizing the weighted sum of squared deviations
        (RMSD) according to the algorithm by Kabsch [8].
        Otherwise, and if ndims is 3, the quaternion based algorithm by Horn [9]
        is used, which is slower when using this Python implementation.

        The returned matrix performs rotation, translation and uniform scaling
        (if specified).

        >>> v0 = [[0, 1031, 1031, 0], [0, 0, 1600, 1600]]
        >>> v1 = [[675, 826, 826, 677], [55, 52, 281, 277]]
        >>> affine_matrix_from_points(v0, v1)
        array([[  0.14549,   0.00062, 675.50008],
            [  0.00048,   0.14094,  53.24971],
            [  0.     ,   0.     ,   1.     ]])
        >>> T = translation_matrix(numpy.random.random(3)-0.5)
        >>> R = random_rotation_matrix(numpy.random.random(3))
        >>> S = scale_matrix(random.random())
        >>> M = concatenate_matrices(T, R, S)
        >>> v0 = (numpy.random.rand(4, 100) - 0.5) * 20
        >>> v0[3] = 1
        >>> v1 = numpy.dot(M, v0)
        >>> v0[:3] += numpy.random.normal(0, 1e-8, 300).reshape(3, -1)
        >>> M = affine_matrix_from_points(v0[:3], v1[:3])
        >>> numpy.allclose(v1, numpy.dot(M, v0))
        True

        More examples in superimposition_matrix()

        """
        v0 = np.array(v0, dtype=np.float64, copy=True)
        v1 = np.array(v1, dtype=np.float64, copy=True)

        ndims = v0.shape[0]
        if ndims < 2 or v0.shape[1] < ndims or v0.shape != v1.shape:
            raise ValueError("input arrays are of wrong shape or type")

        # move centroids to origin
        t0 = -np.mean(v0, axis=1)
        M0 = np.identity(ndims + 1)
        M0[:ndims, ndims] = t0
        v0 += t0.reshape(ndims, 1)
        t1 = -np.mean(v1, axis=1)
        M1 = np.identity(ndims + 1)
        M1[:ndims, ndims] = t1
        v1 += t1.reshape(ndims, 1)

        if shear:
            # Affine transformation
            A = np.concatenate((v0, v1), axis=0)
            u, s, vh = np.linalg.svd(A.T)
            vh = vh[:ndims].T
            B = vh[:ndims]
            C = vh[ndims : 2 * ndims]
            t = np.dot(C, np.linalg.pinv(B))
            t = np.concatenate((t, np.zeros((ndims, 1))), axis=1)
            M = np.vstack((t, ((0.0,) * ndims) + (1.0,)))
        elif usesvd or ndims != 3:
            # Rigid transformation via SVD of covariance matrix
            u, s, vh = np.linalg.svd(np.dot(v1, v0.T))
            # rotation matrix from SVD orthonormal bases
            R = np.dot(u, vh)
            if np.linalg.det(R) < 0.0:
                # R does not constitute right handed system
                R -= np.outer(u[:, ndims - 1], vh[ndims - 1, :] * 2.0)
                s[-1] *= -1.0
            # homogeneous transformation matrix
            M = np.identity(ndims + 1)
            M[:ndims, :ndims] = R
        else:
            # Rigid transformation matrix via quaternion
            # compute symmetric matrix N
            xx, yy, zz = np.sum(v0 * v1, axis=1)
            xy, yz, zx = np.sum(v0 * np.roll(v1, -1, axis=0), axis=1)
            xz, yx, zy = np.sum(v0 * np.roll(v1, -2, axis=0), axis=1)
            N = [
                [xx + yy + zz, 0.0, 0.0, 0.0],
                [yz - zy, xx - yy - zz, 0.0, 0.0],
                [zx - xz, xy + yx, yy - xx - zz, 0.0],
                [xy - yx, zx + xz, yz + zy, zz - xx - yy],
            ]
            # quaternion: eigenvector corresponding to most positive eigenvalue
            w, V = np.linalg.eigh(N)
            q = V[:, np.argmax(w)]
            q /= clf.vector_norm(q)  # unit quaternion
            # homogeneous transformation matrix
            M = clf.quaternion_matrix(q)

        if scale and not shear:
            # Affine transformation; scale is ratio of RMS deviations from centroid
            v0 *= v0
            v1 *= v1
            M[:ndims, :ndims] *= math.sqrt(np.sum(v1) / np.sum(v0))

        # move centroids back
        M = np.dot(np.linalg.inv(M1), np.dot(M, M0))
        M /= M[ndims, ndims]
        return M
