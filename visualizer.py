# Author: Sam Lehmann
# Network with him at: https://www.linkedin.com/in/samuellehmann/
# Date: 2022-11-27
# Description:
# Allows the creation of graphs for analysis of bolted joints

import numpy as np
import matplotlib.pyplot as plt
import fastener_toolkit as ft

MAX_SAFETY_FACTOR = 5


def gen_preload_plot(results, ISO=True):
    """
    Creates a plot of the safety factor of the bolted joint as a function of preload
    :param results: An np array with columns [preload, yield FOS, joint separation FOS, fatigue FOS, minimum FOS]
    :param ISO: Boolean - Is the fastener metric?
    :return:
    """
    # Determine the optimal preload
    max_safety_factor = np.amax(results[:, 4])
    best_preload = results[np.where(results[:, 4] == np.amax(results[:, 4]))[0][0], 0]

    print(f"Max safety factor is {round(max_safety_factor, 3)} with a preload value of {round(best_preload, 3)} N")

    plt.plot(results[:, 0], results[:, 1], label="Factor Of Safety Against Yield")
    plt.plot(results[:, 0], results[:, 2], label="Factor Of Safety Against Joint Separation")
    plt.plot(results[:, 0], results[:, 3], label="Factor Of Safety Against Fatigue Failure")
    plt.scatter(best_preload, max_safety_factor, marker="*", color="black", label="Highest Factor Of Safety", zorder=2)
    plt.gca().axhspan(0, 1, color='red', zorder=1, alpha=0.1)
    plt.gca().axhline(1, color='red', zorder=2)
    if ISO:
        plt.xlabel("Preload [N]")
    else:
        plt.xlabel("Preload [lbf]")
    plt.ylabel("Factor of Safety")
    plt.title("Bolt Factor of Safety as a Function of Preload")
    plt.gca().set_ylim(ymin=0, ymax=MAX_SAFETY_FACTOR)
    plt.gca().set_xlim(xmin=0, xmax=np.max(results[:,0]))
    plt.legend()
    plt.show()


def gen_proof_percentage_plot(results, proof_strength, a_ts):
    """
    Creates a plot of the safety factor of the bolted joint as a function of the preload percentage of proof
    :param results: An np array with columns [preload, yield FOS, joint separation FOS, fatigue FOS, minimum FOS]
    :param proof_strength: The proof strength of the fastener [MPa or psi]
    :param a_ts: The tensile stress area of the bolt [MPa or psi]
    :return:
    """

    # Convert results into percentage of proof strength
    results[:, 0] /= (a_ts * proof_strength)
    results[:, 0] *= 100

    # Determine the optimal safety factor
    max_safety_factor = np.amax(results[:, 4])
    best_preload = results[np.where(results[:, 4] == np.amax(results[:, 4]))[0][0], 0]

    print(f"Max safety factor is {round(max_safety_factor, 3)} with a preload percentage of proof of "
          f"{round(best_preload, 3)}%")

    plt.plot(results[:, 0], results[:, 1], label="Factor Of Safety Against Yield")
    plt.plot(results[:, 0], results[:, 2], label="Factor Of Safety Against Joint Separation")
    plt.plot(results[:, 0], results[:, 3], label="Factor Of Safety Against Fatigue Failure")
    plt.scatter(best_preload, max_safety_factor, marker="*", color="black", label="Highest Factor Of Safety",
                zorder=2)
    plt.gca().axhspan(0, 1, color='red', zorder=1, alpha=0.1)
    plt.gca().axhline(1, color='red', zorder=2)
    plt.xlabel("Percentage of Proof Load %")
    plt.ylabel("Factor of Safety")
    plt.title("Bolt Factor of Safety as a Function of Preload Percentage Of Proof Load")
    plt.gca().set_ylim(ymin=0, ymax=MAX_SAFETY_FACTOR)
    plt.gca().set_xlim(xmin=0, xmax=100)
    plt.legend()
    plt.show()


def gen_joint_diagram(c, k_b, preload, load, ISO=True):
    """
    Creates a standard fastener joint diagram showing changes in deflection and force as load and preload are applied
    :param c: The joint constant
    :param k_b: The bolt stiffness [N/mm or lbf/in]
    :param preload: The preload applied to the fastener [N or lbf]
    :param load: The load applied to the joint [N or lbf]
    :param ISO: Boolean - Is the fastener metric?
    :return:
    """

    bolt_load, member_load = ft.segregate_loads(c, load)

    # Determine member stiffness with reference to EQ. 15.13c [Norton]
    k_m = k_b / c - k_b

    # Determine deflection values
    defl_member_1 = -preload / k_m
    defl_member_2 = -(preload - member_load) / k_m
    defl_bolt_1 = preload / k_b
    defl_bolt_2 = (bolt_load + preload) / k_b

    plt.plot([0, defl_bolt_1, defl_bolt_2], [0, preload, preload + bolt_load], label="Bolt Line")
    plt.plot([0, defl_member_1, defl_member_2], [0, preload, preload - member_load], label="Member Line")

    plt.annotate("Preload Point", [defl_bolt_1, preload], color="black", weight="bold", xytext=(-100, -20),
                 textcoords="offset points", fontsize=10, va="center", ha='center',
                 arrowprops=dict(arrowstyle="->"))
    plt.annotate("Applied Load Point", [defl_bolt_2, preload + bolt_load], weight="bold", xytext=(-70, 0),
                 textcoords="offset points",
                 color="black", fontsize=10, va="center", ha='center',
                 arrowprops=dict(arrowstyle="->"))

    plt.annotate("Preload Point", [defl_member_1, preload], weight="bold", color="black", xytext=(70, 0),
                 textcoords="offset points", fontsize=10, va="center", ha='center',
                 arrowprops=dict(arrowstyle="->"))
    plt.annotate("Applied Load Point", [defl_member_2, preload - member_load], weight="bold", xytext=(100, -20),
                 textcoords="offset points",
                 color="black", fontsize=10, va="center", ha='center',
                 arrowprops=dict(arrowstyle="->"))

    plt.legend()
    plt.title("Joint Diagram")
    if ISO:
        plt.xlabel("Displacement [mm]")
        plt.ylabel("Load [N]")
    else:
        plt.xlabel("Displacement [in]")
        plt.ylabel("Load [lbf]")

    plt.gca().set_ylim(ymin=0)

    # Add extension and compression arrows
    label_y = preload / 50
    prop = dict(arrowstyle="-|>,head_width=0.4,head_length=0.8", shrinkA=0, shrinkB=0, color = "red")
    plt.annotate("", xy=(defl_member_1, 0), xytext=(0, 0), arrowprops=prop)
    prop = dict(arrowstyle="-|>,head_width=0.4,head_length=0.8", shrinkA=0, shrinkB=0, color="green")
    plt.annotate("", xy=(defl_bolt_2, 0), xytext=(0, 0), arrowprops=prop)
    plt.annotate("Extension", ha = "right", xy = (defl_bolt_2,label_y), weight="bold")
    plt.annotate("Compression", xy=(defl_member_1, label_y), weight="bold")

    plt.show()
