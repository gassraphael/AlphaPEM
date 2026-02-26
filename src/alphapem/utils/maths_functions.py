# -*- coding: utf-8 -*-

"""This module contains mathematical functions which are used for modeling the PEM fuel cell."""

# ________________________________________________Mathematical functions________________________________________________

def hmean(terms, weights=None):
    """
    Calculate the weighted harmonic mean of a list of terms with corresponding weights.
    It is more efficient to express this function in the code than calling hmean from scipy.stats.

    Parameters
    ----------
    terms (list of float):
        The terms to calculate the harmonic mean for.
    weights (list of float):
        The weights corresponding to each term. If None, uniform weights are assumed.

    Returns
    -------
    float:
        The weighted harmonic mean.
    """

    n = len(terms)
    if weights is None:
        weights = [1] * n  # Assign equal weights if not provided

    if len(weights) != n:
        raise ValueError("The length of terms and weights must be the same.")

    # Calculate the weighted harmonic mean
    weighted_sum = 0
    total_weight = 0
    for w, t in zip(weights, terms):
        if t != 0:
            weighted_sum += w / t
        total_weight += w

    if weighted_sum == 0:
        return float('inf')  # Avoid division by zero

    return total_weight / weighted_sum


def average(terms, weights=None):
    """
    Calculate the weighted arithmetic mean of a list of terms with corresponding weights.
    It is more efficient to express this function in the code than calling average from numpy.

    Parameters
    ----------
    terms (list of float):
        The terms to calculate the average for.
    weights (list of float, optional):
        The weights corresponding to each term. If None, uniform weights are assumed.

    Returns
    -------
    float:
        The weighted arithmetic mean.
    """
    n = len(terms)
    if weights is None:
        total_weight = n
        weighted_sum = 0.0
        for t in terms:
            weighted_sum += t
    else:
        if n != len(weights):
            raise ValueError("The length of terms and weights must be the same.")
        total_weight = 0.0
        weighted_sum = 0.0
        for i in range(n):
            w = weights[i]
            total_weight += w
            weighted_sum += w * terms[i]

    if total_weight == 0:
        return float('nan')
    return weighted_sum / total_weight


def interpolate(terms, distances):
    """
    Fast inverse distance interpolation for exactly 2 points.

    Parameters
    ----------
    terms : list of float
        The values at each node ([y1, y2]).
    distances : list of float
        The distances from each node to the interpolation point ([d1, d2]).

    Returns
    -------
    float
        The interpolated value.
    """
    if len(terms) != 2 or len(distances) != 2:
        raise ValueError("This function only supports interpolation with 2 points.")
    y1, y2 = terms
    d1, d2 = distances
    if d1 == 0: return y1
    if d2 == 0: return y2
    return (d2 * y1 + d1 * y2) / (d1 + d2)


def d_dx(y_minus, y_plus, dx = None, dx_minus = None, dx_plus = None):
    """
    Computes the centered first derivative (second order) with different steps to the left and right.

    Parameters
    ----------
    y_minus : float
        Value at the left point (i-1).
    y_plus : float
        Value at the right point (i+1).
    dx : float
        Half step between (i-1) and (i+1) when dx_minus = dx_plus.
    dx_minus : float
        Step between (i-1) and i.
    dx_plus : float
        Step between i and (i+1).

    Returns
    -------
    float
        Approximation of the first derivative at i.
    """

    # Case of uniform grid spacing
    if dx is None:
        if dx_minus is None or dx_plus is None:
            raise ValueError("Either dx or both dx_minus and dx_plus must be provided.")
    else:
        if dx == 0:
            raise ValueError("dx must be non-zero.")
        return (y_plus - y_minus) / (2.0 * dx)

    # Case of non-uniform grid spacing (dx is None and dx_minus and dx_plus are provided)
    if dx_minus <= 0 or dx_plus <= 0:
        raise ValueError("dx_minus and dx_plus must be positive non-zero values.")
    y_0 = interpolate([y_minus, y_plus], [dx_minus, dx_plus])
    return (y_plus * dx_minus**2 + y_0 * (dx_plus**2 - dx_minus**2) - y_minus * dx_plus**2)  / \
           (dx_minus * dx_plus * (dx_minus + dx_plus))
