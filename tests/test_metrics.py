import torch

from perturbmodel.evaluation import e_distance


def test_e_distance_zero_for_identical_distributions():
    g = torch.Generator().manual_seed(0)
    x = torch.randn(500, 10, generator=g)
    y = torch.randn(500, 10, generator=g)
    assert abs(e_distance(x, y)) < 0.05


def test_e_distance_positive_for_shifted_distribution():
    g = torch.Generator().manual_seed(0)
    x = torch.randn(500, 10, generator=g)
    y = torch.randn(500, 10, generator=g) + 2.0
    assert e_distance(x, y) > 1.0
