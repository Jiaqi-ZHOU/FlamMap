from __future__ import annotations

import numpy as np
from scipy.constants import R


def fit_nasa7_poly(temps, cp_array, hf_array, s_array):
    res = np.polyfit(temps, cp_array / R, deg=4)
    a5, a4, a3, a2, a1 = res
    a6_over_t = hf_array / (R * temps) - (
        a1 + a2 * temps / 2 + a3 * temps**2 / 3 + a4 * temps**3 / 4 + a5 * temps**4 / 5
    )
    x = 1 / temps
    a6 = np.sum(x * a6_over_t) / np.sum(x * x)
    a7 = np.mean(
        s_array / R
        - (
            a1 * np.log(temps)
            + a2 * temps
            + a3 * temps**2 / 2
            + a4 * temps**3 / 3
            + a5 * temps**4 / 4
        )
    )
    return [a1, a2, a3, a4, a5, a6, a7]
