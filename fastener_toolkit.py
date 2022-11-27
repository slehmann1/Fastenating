# Author: Sam Lehmann
# Network with him at: https://www.linkedin.com/in/samuellehmann/
# Date: 2022-11-27
# Description:
# A toolkit of common fastener calculations, based on the methodology outlined in
# Norton, Robert L. Machine Design: An Integrated Approach. Pearson, 2020.

import math
import unittest
import pandas as pd
import numpy as np

# Parameters for the Cornwell equation for joints with two plates of the same material as described by Norton for
# equation 15.19 and presented in table 15-8
CORNWELL_PARAMS = pd.DataFrame(
    np.array([[0.1, 0.4389, -0.9197, 0.8901, -0.3187],
              [0.2, 0.6118, -1.1715, 1.0875, -0.3806],
              [0.3, 0.6932, -1.2426, 1.1177, -0.3845],
              [0.4, 0.7351, -1.2612, 1.1111, -0.3779],
              [0.5, 0.758, -1.2632, 1.0979, -0.3708],
              [0.6, 0.7709, -1.2600, 1.0851, -0.3647],
              [0.7, 0.7773, -1.2543, 1.0735, -0.3595],
              [0.8, 0.78, -1.2503, 1.0672, -0.3571],
              [0.9, 0.7797, -1.2458, 1.062, -0.3552],
              [1.0, 0.7774, -1.2413, 1.0577, -0.3537],
              [1.25, 0.7667, -1.2333, 1.0548, -0.3535],
              [1.5, 0.7518, -1.2264, 1.0554, -0.3550],
              [1.75, 0.735, -1.2202, 1.0581, -0.3574],
              [2.00, 0.7175, -1.2133, 1.0604, -0.3596]]),
    columns=["j", "p0", "p1", "p2", "p3"])


def get_tensile_stress_area(d_major, d_minor, pitch=None, num_threads=None):
    """
    Returns the tensile stress area of the bolt rounded to 3 decimal places
    :param d_major: Major diameter of the bolt [mm or in]
    :param d_minor: Major diameter of the bolt [mm or in]
    :param pitch: Pitch of the bolt [mm or in]
    :param num_threads: Number of threads per inch [npi]
    :return: the tensile stress area of the bolt [mm^2}
    """

    # Reference eq 15.1 from Norton
    if pitch is not None:
        d_pitch = d_major - 0.649519 * pitch

    else:
        d_pitch = d_major - 0.649519 / num_threads

    return math.pi / 4 * ((d_pitch + d_minor) / 2) ** 2

def get_bolt_stiffness(a_ts, a_cs, l_unthreaded, l_threaded, E_b):
    """

    :param a_ts: Tensile stress area of the bolt [mm^2 or in^2]
    :param a_cs: Total cross sectional area of the bolt [mm or in]
    :param l_unthreaded: Length of the unthreaded shank within the grip zone[mm or in]
    :param l_threaded: Length of the threaded shank within the grip zone [mm or in]
    :param E_b: Young's modulus of the bolt [MPa or psi]
    :return:The spring constant of the bolt [N/mm or lbf/in]
    """
    return a_ts * a_cs * E_b / (a_cs * l_threaded + a_ts * l_unthreaded)


def get_joint_constant(d_b, l, E_m, E_b):
    """
    Determines the stiffness of clamped members based on the Cornwell method. Makes the assumption that both members
    are of an equivalent material. Reference Eq 15.19 from Norton
    :param d_b: The diameter of the bolt [mm or in]
    :param l: The clamped length of the joint [mm or in]
    :param E_m: The member material Young's modulus [MPa or psi]
    :param E_b: The bolt material Young's modulus [MPa or psi]
    :return: The member stiffness [N/mm or lbf/in]
    """
    j = d_b / l

    # Linearly interpolate through the table to obtain values for p0, p1, p2, and p3

    j_2 = CORNWELL_PARAMS[CORNWELL_PARAMS['j'] >= j - 0.01].iloc[0, 0]

    # Check if j_2 is the max value in the table, if yes, then there is no need to linearly interpolate
    if (j_2 > 1.75):
        j_1 = j_2
    else:
        j_1 = CORNWELL_PARAMS[CORNWELL_PARAMS['j'] <= j + 0.01].iloc[-1, 0]

    # Bound the values for the range of J values published
    j_1 = bound_val(j_1, [0.1, 2])
    j_2 = bound_val(j_2, [0.1, 2])

    p_0 = linterp(j_1, j_2, float(CORNWELL_PARAMS.loc[CORNWELL_PARAMS["j"] == j_1]['p0']),
                  float(CORNWELL_PARAMS.loc[CORNWELL_PARAMS["j"] == j_2]['p0']), j)

    p_1 = linterp(j_1, j_2, float(CORNWELL_PARAMS.loc[CORNWELL_PARAMS["j"] == j_1]['p1']),
                  float(CORNWELL_PARAMS.loc[CORNWELL_PARAMS["j"] == j_2]['p1']), j)

    p_2 = linterp(j_1, j_2, float(CORNWELL_PARAMS.loc[CORNWELL_PARAMS["j"] == j_1]['p2']),
                  float(CORNWELL_PARAMS.loc[CORNWELL_PARAMS["j"] == j_2]['p2']), j)

    p_3 = linterp(j_1, j_2, float(CORNWELL_PARAMS.loc[CORNWELL_PARAMS["j"] == j_1]['p3']),
                  float(CORNWELL_PARAMS.loc[CORNWELL_PARAMS["j"] == j_2]['p3']), j)

    # Plate to modulus ratio, r
    r = E_m / E_b

    c = p_3 * r ** 3 + p_2 * r ** 2 + p_1 * r + p_0

    return c


def segregate_loads(c, load):
    """
    Identifies the quantity of a load carried by the bolt and by the members
    :param c: Joint stiffness [N/mm or lbf/in]
    :param load: Load applied to the bolt [N or lbf]
    :return: load borne by the bolt, load borne by the members [N or lbf]
    """

    p_b = c * load
    p_m = (1 - c) * load

    return p_b, p_m


def bolt_yield_safety_factor(c, load, preload, a_ts, b_ys):
    """
    Determines the factor of safety against yielding the bolt under statically applied tension load
    :param c: Joint constant [N/mm or lbf/in]
    :param load: Load applied to the joint [N or lbf]
    :param preload: Preload applied to the bolt [N or lbf]
    :param a_ts: Tensile stress area of the bolt [mm^2 or in^2]
    :param b_ys: Yield strength of the bolt [MPa or psi]
    :return:Factor of safety against yielding under statically applied tension load
    """
    # Portion of load carried by the bolt
    f_b = segregate_loads(c, load)[0] + preload

    # Stress in bolt
    sigma_b = f_b / a_ts

    # Factor of safety against yield
    n_y = b_ys / sigma_b

    return n_y


def joint_separation_safety_factor(c, load, preload):
    """
    Determines the factor of safety against joint separation
    :param c:Joint constant [N/mm or lbf/in]
    :param load:Load applied to the joint [N or lbf]
    :param preload:Preload applied to the join [N or lbf]
    :return:The factor of safety
    """
    p_0 = preload / (1 - c)

    return p_0 / load


def fatigue_safety_factor(d, c, max_load, min_load, preload, a_ts, s_ut, s_end, s_y, ISO=True):
    """
    Determines factor of safety in fatigue for a bolted joint loaded with alternating tension loads through the use of a
    modified goodman diagram
    :param d: Bolt diameter [mm or in]
    :param c: Joint constant [N/mm or lbf/in]
    :param max_load: Maximum tensile load applied to the joint [N or lbf]
    :param min_load: Minimum tensile load applied to the joint [N or lbf]
    :param preload: Preload applied to the bolt [N or lbf]
    :param a_ts: Tensile stress area [mm^2 or in^2]
    :param s_ut: Ultimate tensile strength of the bolt [MPa or psi]
    :param s_y: Yield strength of the bolt [MPa or psi]
    :param s_end: Endurance limit of the bolt [MPa or psi]
    :param ISO: Boolean - Is the fastener metric?
    :return:The factor of safety
    """
    f_b_max = segregate_loads(c, max_load)[0] + preload
    f_b_min = segregate_loads(c, min_load)[0] + preload

    f_mean = (f_b_max + f_b_min) / 2
    f_alt = (f_b_max - f_b_min) / 2

    sigma_mean = f_mean / a_ts
    sigma_alt = f_alt / a_ts
    sigma_preload = preload / a_ts

    # Determine stress concentration factors
    if ISO:
        # Stress concentration factor as per eq 15.15C
        k_f = 5.7 + 0.02682 * d
    else:
        k_f = 5.7 + 0.6812 * d

    if isinstance(sigma_mean, float):
        k_fm = det_mean_stress_concentration_factor(sigma_mean, sigma_alt, s_y, k_f)
    else:
        k_fm = np.zeros(len(sigma_mean))
        for i in range(0, len(sigma_mean)):
            k_fm[i] = det_mean_stress_concentration_factor(sigma_mean[i], sigma_alt[i], s_y, k_f)

    sigma_mean *= k_fm
    sigma_preload *= k_fm
    sigma_alt *= k_f

    # Factor of safety in fatigue as per the goodman methodology
    n_f = s_end * (s_ut - sigma_preload) / (s_end * (sigma_mean - sigma_preload) + s_ut * sigma_alt)
    return n_f


def det_mean_stress_concentration_factor(sigma_mean, sigma_alt, s_y, k_f):
    """
    Determines the mean stress concentration factor when a body is subjected to alternating loads, Ref eq 6.17 [Norton]
    :param sigma_mean: Mean stress [MPa or psi]
    :param sigma_alt: Alternating stress [MPa or psi]
    :param s_y: Material yield strength [MPa or psi]
    :param k_f: Fatigue stress concentration factor
    :return: Mean stress concentration factor
    """
    if (k_f * (sigma_mean + sigma_alt)) < s_y:
        k_fm = k_f
    elif k_f * (sigma_mean + sigma_alt) > s_y:
        k_fm = (s_y - k_f * sigma_alt) / sigma_mean
    else:
        k_fm = 0

    return k_fm


def bound_val(val, limits):
    """
    Bounds a value within limits, setting its value to the upper or lower limit if it is out of bounds
    :param val: The value to be limited
    :param limits: An array of two values, a lower and upper limit
    :return: The bounded value
    """
    if val < limits[0]:
        val = limits[0]
    elif val > limits[1]:
        val = limits[1]

    return val


def linterp(x1, x2, y1, y2, x3):
    """
    Linearly interpolates to determine y3 for a range of values defined by [x1,y1], [x2, y2], [x3, y3 = ?]
    return: y3
    """

    if x2 == x1:
        return y1

    m = (y2 - y1) / (x2 - x1)
    return (x3 - x1) * m + y1


class TestFastenerToolkit(unittest.TestCase):
    # Description:
    # Unit test cases allowing verification of the fastener toolkit using the examples provided by Norton

    def test_get_tensile_stress_area(self):
        """
        Test the tensile stress area function with standard values for an M3 Bolt
        :return:
        """
        actual = round(get_tensile_stress_area(0.5, 3, 2.39),3)
        expected = 5.038
        self.assertEqual(actual, expected)

    def test_get_bolt_spring_constant(self):
        """
        Test the bolt spring constant using the test case provided by Norton in example 15-2
        :return:
        """
        actual = round(get_bolt_stiffness(0.052431, (math.pi * (0.3125 / 2) ** 2), 1.625, 0.375, 30E6),3)
        norton = 1058613.179
        self.assertEqual(norton, actual)

    def test_safety_factors(self):
        """
        Test the safety factors: yield, fatigue, and joint separation using the values of example 15-3 from Norton
        :return:
        """
        min_load_app = 0
        max_load_app = 1000
        s_end = 25726
        s_ut = 120000
        b_ys = 92000
        d = 5 / 16
        preload = 4011
        a_ts = get_tensile_stress_area(d, 0.24033, num_threads=18)

        # This example utilizes a different methodology to determine the joint constant. Use the value provided by
        # Norton
        c = 0.09056

        norton_y =1.18
        norton_f = 1.37
        norton_j = 4.4

        n_f = fatigue_safety_factor(d,c,max_load_app, min_load_app,preload,a_ts,s_ut,s_end,b_ys,False)
        n_y = bolt_yield_safety_factor(c,max_load_app,preload,a_ts,b_ys)
        n_j = joint_separation_safety_factor(c,max_load_app,preload)

        self.assertEqual(round(n_f,2), norton_f)
        self.assertEqual(round(n_y,2), norton_y)
        self.assertEqual(round(n_j,1), norton_j)

