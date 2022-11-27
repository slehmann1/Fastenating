# Author: Sam Lehmann
# Network with him at: https://www.linkedin.com/in/samuellehmann/
# Date: 2022-11-27

import math
import fastener_toolkit as ft
import visualizer as vis
import numpy as np


def sample_case():
    """
    Generates graphs and performs calculations for a sample case for a 5/16-18 UNC-2A steel bolt with class 5.2
    clamping steel members with a grip length of 3 inch whilst subjected to axial cyclic loads ranging between 0 lbf
    and 1000lb.
    :return: None
    """

    # Test Case parameters
    max_preload = 4456
    num_samples = 1000
    min_load_app = 0 # [lbf]
    max_load_app = 1000 # [lbf]
    s_end = 25726 # Endurance limit of bolt [psi]
    s_ut = 120000 # Ultimate tensile strength of bolt [psi]
    b_ys = 92000 # Yield strength of bolt [psi]
    s_p = 80000 # Proof load of bolt [psi]
    d = 5 / 16 # [in]

    E_m = 30E6 # Member Young's modulus [psi]
    E_b = 30E6 # Bolt Young's modulus [psi]
    l = 3 # Grip length [in]
    l_threaded = 2 # Threaded length within the grip zone [in]

    # Determine joint constant and bolt stiffness
    a_ts = ft.get_tensile_stress_area(d, 0.24033, num_threads=18)  # Minor diameter is from a lookup chart [in]
    c = ft.get_joint_constant(d,l,E_m, E_b)
    k_b = ft.get_bolt_stiffness(a_ts, (math.pi * (d / 2) ** 2), l-l_threaded, l_threaded, E_b)

    # Create results, an np array with columns [preload, yield FOS, joint separation FOS, fatigue FOS, minimum FOS]
    results = np.vstack((np.arange(0, max_preload, (max_preload / num_samples)), np.zeros([1, num_samples]),
                         np.zeros([1, num_samples]), np.zeros([1, num_samples]),
                         np.zeros([1, num_samples]))).transpose()
    results[:, 1] = ft.bolt_yield_safety_factor(c, max_load_app, results[:, 0], a_ts, b_ys)
    results[:, 2] = ft.joint_separation_safety_factor(c, max_load_app, results[:, 0])
    results[:, 3] = ft.fatigue_safety_factor(d, c, max_load_app, min_load_app, results[:, 0], a_ts, s_ut, s_end, b_ys,
                                             False)
    results[:, 4] = np.amin(results[:, 1:4], 1)

    best_preload = results[np.where(results[:, 4] == np.amax(results[:, 4]))[0][0], 0]

    # Plot a joint diagram for the optimal preload, where the maximum load is applied
    vis.gen_joint_diagram(c, k_b, best_preload, max_load_app, ISO=False)
    # Generate charts
    vis.gen_preload_plot(results, ISO=False)
    vis.gen_proof_percentage_plot(results, s_p, a_ts)



# Press the green button in the gutter to run the script.
if __name__ == '__main__':
    sample_case()
