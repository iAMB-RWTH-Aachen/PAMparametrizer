import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from matplotlib.lines import Line2D
import matplotlib.container as mcontainer
import pandas as pd
import os
import seaborn as sns

from Scripts.i3_analysis.metabolic_flux_distribution_vs_exp import main_iJN1463 as plot_intracell_flux_distribution_pp
from Scripts.i3_analysis.metabolic_flux_distribution_vs_exp import main_iCGB21FR as plot_intracell_flux_distribution_cgb
from Scripts.Data_requirements.analyse_data_reduction_results import plot_progression_of_errors, plot_deviation_of_error


def extract_color_and_marker(handle):
    """Extract line color and marker style from different matplotlib handle types."""
    color, marker, mface, medge = None, None, None, None

    if isinstance(handle, mcontainer.ErrorbarContainer):
        # Extract the main line and markerline
        if handle.lines[0]:  # main line
            color = handle.lines[0].get_color()
        if handle.lines[1]:  # markerline
            marker = handle.lines[1][0].get_marker()
            mface = handle.lines[1][0].get_markerfacecolor()
            medge = handle.lines[1][0].get_markeredgecolor()

    elif isinstance(handle, Line2D):
        # Regular Line2D (from plot or scatter)
        color = handle.get_color()
        marker = handle.get_marker()
        mface = handle.get_markerfacecolor()
        medge = handle.get_markeredgecolor()

    else:
        # fallback (e.g., Patch, Collection)
        try:
            color = handle.get_facecolor()[0]
        except Exception:
            color = 'black'

    return color, marker, mface, medge



FONTSIZE = 11
def main():
    fig = plt.figure(figsize=(15, 8))  # Adjust overall size as needed
    gs = gridspec.GridSpec(4, 6, height_ratios=[1.2, 0.6,1,0.3],
                           width_ratios=[0.05, 1, 1.5, 0.2, 1.4, 0.05], hspace=0.1, wspace=0.3)

    ax_b = plot_intracell_flux_distribution_pp(gs = gs[0, 1:], fig = fig, vrange=(-6,13))

    ax_c = plot_intracell_flux_distribution_cgb(gs = gs[2, 1:3], fig = fig, vrange=(-6,13), cbar = False)
    # ax_c.set_ylabel('')

    ax_a = fig.add_subplot(gs[2, 4:])
    final_errors = pd.read_excel(os.path.join('Results', 'data_reduction_results', 'r_squared_for_analysis.xlsx'))
    plot_progression_of_errors(final_errors, ax = ax_a, legend=False)

    ax_legend = fig.add_subplot(gs[3, 4:])
    ax_legend.axis("off")
    h, l = ax_a.get_legend_handles_labels()

    ax_legend.legend(h, l, loc="upper left",
                     fontsize=FONTSIZE-3,
                     ncol=2, frameon=False, bbox_to_anchor=(-0.1, 0.5))

    annotations = ["A",'','','', '',"B",'','',"C"]
    fontsize = FONTSIZE  # Adjust as needed

    for ax, label in zip(fig.axes, annotations):
        ax.annotate(label, xy=(-0.05, 1), xycoords="axes fraction",
                    fontsize=fontsize, fontweight='bold',
                    xytext=(-5, 5), textcoords="offset points",
                    ha="right", va="bottom")

    fig.tight_layout()
    fig.savefig(os.path.join('Figures', 'Figure4_data_reduction.png'))

def main_datareduc():
    fig = plt.figure(figsize=(10.5/2.56, 15/2.56))  # Adjust overall size as needed
    gs = gridspec.GridSpec(3, 1, hspace = 0.2, height_ratios=[1,1,0.05])

    ax_a = fig.add_subplot(gs[0])
    ax_b = fig.add_subplot(gs[1])
    final_errors = pd.read_excel(os.path.join('Results', 'data_reduction_results', 'r_squared_for_analysis.xlsx'))
    plot_progression_of_errors(final_errors, ax = ax_a, legend=False)
    plot_deviation_of_error(final_errors, ax = ax_b, legend=False)

    # ax_a = fig.add_subplot(gs[0])

    ax_legend = fig.add_subplot(gs[-1])
    ax_legend.axis("off")
    h_a, l_a = ax_a.get_legend_handles_labels()  # lineplot
    h_b, l_b = ax_b.get_legend_handles_labels()  # scatter

    # Combine by matching labels
    combined_handles = []
    combined_labels = []

    for label in set(l_a) & set(l_b):
        handle_a = h_a[l_a.index(label)]
        handle_b = h_b[l_b.index(label)]

        color_a, _, _, _ = extract_color_and_marker(handle_a)
        _, marker_b, mface_b, medge_b = extract_color_and_marker(handle_b)

        combined_handles.append(Line2D(
            [0], [0],
            color=color_a,
            marker=marker_b or 'o',
            markerfacecolor=mface_b or color_a,
            markeredgecolor=medge_b or 'black',
            linestyle='-',
            linewidth=2,
            markersize=6
        ))
        combined_labels.append(label)

    # Draw the combined legend
    ax_legend.legend(
        combined_handles, combined_labels,
        loc="upper left",
        fontsize=FONTSIZE,
        ncol=2,
        frameon=False,
        bbox_to_anchor=(-0.1, 0.5)
    )

    annotations = ["A","B"]
    fontsize = FONTSIZE  # Adjust as needed

    for ax, label in zip(fig.axes, annotations):
        ax.annotate(label, xy=(-0.15, 1), xycoords="axes fraction",
                    fontsize=fontsize, fontweight='bold',
                    xytext=(-5, 5), textcoords="offset points",
                    ha="right", va="bottom")

    # fig.tight_layout()
    fig.savefig(os.path.join('Figures', 'Figure4_data_reduction.png'))

if __name__ == "__main__":
    main_datareduc()
