import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import pandas as pd
import os
import seaborn as sns

from Scripts.i3_analysis.metabolic_flux_distribution_vs_exp import main_iJN1463 as plot_intracell_flux_distribution_pp
from Scripts.i3_analysis.metabolic_flux_distribution_vs_exp import main_iCGB21FR as plot_intracell_flux_distribution_cgb
from Scripts.Data_requirements.analyse_data_reduction_results import plot_progression_of_errors

FONTSIZE = 16
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

if __name__ == "__main__":
    main()
