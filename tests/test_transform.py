import pytest
from pytest import approx
import numpy as np
from autovps.transform import Transform

SQRT2 = np.sqrt(2)/2

def test_rotx():

    I = Transform(np.eye(4))
    I.rotx(45)
    assert (approx(I.tform) == [[1, 0, 0, 0],
        [0, SQRT2, -SQRT2, 0],
        [0, SQRT2, SQRT2, 0],
        [0, 0, 0, 1]])


def test_roty():
    I = Transform(np.eye(4))
    I.roty(45)
    assert approx(I.tform) == [[SQRT2, 0, SQRT2, 0],
        [0, 1, 0, 0],
        [-SQRT2, 0, SQRT2, 0],
        [0, 0, 0, 1]]


def test_rotz():
    I = Transform(np.eye(4))
    I.rotz(45)
    assert approx(I.tform) == [[SQRT2, -SQRT2, 0, 0],
        [SQRT2, SQRT2, 0, 0],
        [0, 0, 1, 0],
        [0, 0, 0, 1]]


def test_scale():
    I = Transform(np.eye(4))
    I.scaleXYZ(2, 3, 4)

    assert (I.tform == np.diag([2, 3, 4, 1])).all()
    assert (I.get_scale() == [2, 3, 4]).all()

    I2 = Transform(np.eye(4))
    I2.scale([2, 3, 4])
    assert (I.get_scale() == I2.get_scale()).all()


def test_concat():
    targ = [[SQRT2, 0, SQRT2, 0],
        [.5, SQRT2, -.5, 0],
        [-.5, SQRT2, .5, 0],
        [0, 0, 0, 1]]

    I = Transform(np.eye(4))

    rx = Transform(np.eye(4))
    rx.rotx(45)

    ry = Transform(np.eye(4))
    ry.roty(45)

    I.concatenate_matrix(rx.tform)
    I.concatenate_matrix(ry.tform)

    assert approx(I.tform) == targ

    # overload
    I2 = Transform(np.eye(4))
    X = I2*rx*ry
    assert approx(X.tform) == targ
