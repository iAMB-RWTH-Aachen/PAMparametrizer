from cobra.io.sbml import read_sbml_model
import pandas as pd
import os
import numpy as np
from matplotlib import pyplot as plt
from matplotlib.lines import Line2D
from matplotlib import gridspec
import seaborn as sns
from typing import Union, List, Dict, Literal
from PAModelpy import PAModel
from PAModelpy.utils import set_up_pam
from scipy.stats import gaussian_kde, pearsonr
from matplotlib.gridspec import GridSpec, GridSpecFromSubplotSpec

from cobra.flux_analysis import flux_variability_analysis
from Modules.PAMparametrizer.utils.pam_generation import create_pamodel_from_diagnostics_file
from Figures.Scripts.Figure1_iml1515_kcat_analysis import add_pathway_annotation_to_proteins


PARAM_FILE_OLD = os.path.join('Results', '1_preprocessing','proteinAllocationModel_iML1515_EnzymaticData_250912.xlsx')
PARAM_FILE_PREPROC = os.path.join('Results','2_parametrization','proteinAllocationModel_iML1515_EnzymaticData_multi.xlsx')
MODEL_FILE = os.path.join('Models', 'iML1515.xml')
SUBSTRATE_ID = 'EX_glc__D_e'
GEROSA_GLC_UPTAKE = -9.654

RESULT_FLUX_PATH = os.path.join('Results', '3_analysis', 'iML1515_alternative_models_predictions.xlsx')

COG_DESCRIPTION2LETTER = {
    'Translation, ribosomal structure and biogenesis':'J',
    'RNA processing and modification': 'A',
    'Transcription': 'K',
    'Replication, recombination and repair': 'L',
    'Chromatin structure and dynamics':'B',
    'Cell cycle control, cell division, chromosome partitioning': 'D',
    'Nuclear structure':'Y',
    'Defense mechanisms': 'V',
    'Signal transduction mechanisms': 'T',
    'Cell wall/membrane/envelope biogenesis':'M',
    'Cell motility':'N',
    'Cytoskeleton':'Z',
    'Extracellular structures':'W',
    'Intracellular trafficking, secretion, and vesicular transport':'U',
   'Posttranslational modification, protein turnover, chaperones':'O',
    'Energy production and conversion':'C',
    'Carbohydrate transport and metabolism':'G',
    'Amino acid transport and metabolism': 'E',
    'Nucleotide transport and metabolism':'F',
    'Coenzyme transport and metabolism':'H',
    'Lipid transport and metabolism':'I',
    'Inorganic ion transport and metabolism':'P',
    'Secondary metabolites biosynthesis, transport and catabolism':'Q',
    'General function prediction only':'R',
    'Function unknown':'S',
    'Overall': 'Overall'
}

COG_MAPPER = {'Amino acid transport and metabolism': 'Amino acid metabolism',
       'Carbohydrate transport and metabolism': 'Carbon metabolism',
       'Cell cycle control, cell division, chromosome partitioning': 'Cell cycle',
       'Cell wall/membrane/envelope biogenesis':"Cell membrane biogenesis",
       'Coenzyme transport and metabolism': 'Coenzyme metabolism',
              'Defense mechanisms':'Defense mechanisms',
       'Energy production and conversion': 'Energy generation', 'Function unknown': 'Function unknown',
       'General function prediction only': 'General function',
       'Inorganic ion transport and metabolism': 'Inorganic ion metabolism',
       'Intracellular trafficking, secretion, and vesicular transport': 'Intracellular transport',
       'Lipid transport and metabolism': 'Lipid metabolism',
       'Nucleotide transport and metabolism': 'Nucleotide metabolism',
       'Posttranslational modification, protein turnover, chaperones':'Protein modification',
       'Replication, recombination and repair': 'Replication',
       'Secondary metabolites biosynthesis, transport and catabolism': 'Secondary metabolites metabolism',
       'Signal transduction mechanisms': 'Signal transduction', 'Transcription': 'Transcription',
       'Translation, ribosomal structure and biogenesis': 'Translation',
    'Overall': 'Overall'}

def _merge_on_rxn(df_left, df_right, suffix):
    """Helper – left‑outer merge on ``rxn_id`` and keep the COG column."""
    merged = pd.merge(
        df_left,
        df_right,
        on=["rxn_id", "COG description"],
        how="inner",
        suffixes=("", suffix),
    )
    return merged


def get_all_kcat_values(data_file_paths: list[pd.DataFrame],
                                     label_names:list[str]):
    all_kcat_values = pd.DataFrame(columns = ['rxn_id', 'enzyme_id', 'direction'])
    for label, data_file_path in zip(label_names,data_file_paths):
        aes_parameter_df = pd.read_excel(data_file_path, sheet_name='ActiveEnzymes')[['rxn_id', 'enzyme_id', 'direction', 'kcat_values']]
        all_kcat_values = pd.merge(all_kcat_values, aes_parameter_df, on = ['rxn_id', 'enzyme_id', 'direction'],
                                    how ='outer', suffixes=['', label])
    return all_kcat_values

def get_all_protein_concentrations():
    protein_concs = pd.read_excel(RESULT_FLUX_PATH, sheet_name=None)
    all_protein_concs = pd.DataFrame(columns = ['protein'])
    for model_id, df in protein_concs.items():
        if not model_id.endswith('_protein'): continue
        all_protein_concs = pd.merge(all_protein_concs, df, on = 'protein',
                                    how ='outer', suffixes=['', '+_'+model_id.replace('_protein', '')])
    return all_protein_concs

def get_flux_predictions_and_fva(model_files: List[str],
                                 diagnostic_files: List[str],
                                     label_names:List[str],
                                     simulated_fluxes_result_file: str = RESULT_FLUX_PATH,
                                     ):
    pam = set_up_pam(pam_info_file = PARAM_FILE_PREPROC, model = MODEL_FILE, sensitivity = False)
    models = [set_up_pam(pam_info_file = file, model = MODEL_FILE, sensitivity = False)
              if 'xml' not in file else read_sbml_model(file)
              for file in model_files
    ]
    alt_pams = [create_pamodel_from_diagnostics_file(diag_file, pam.copy_with_pickle()) for diag_file in diagnostic_files]

    for label, model in zip(label_names, models+alt_pams):
        print('------------------------------------------------------------------------')
        print('Analyzing flux distribution of model', label, '\n')

        model.reactions.get_by_id(SUBSTRATE_ID).lower_bound = GEROSA_GLC_UPTAKE
        model.optimize()

        try:
            protein_concentrations = {enz.id: [enz.concentration] for enz in model.enzymes}
        except:
            continue

        all_reactions = [rxn for rxn in model.reactions if 'CE_' not in rxn.id]#CE to ignore the catalytic events
        fluxes = pd.DataFrame.from_dict({rxn.id: rxn.flux for rxn in all_reactions},
                              columns =['flux'], orient='index',)
        fva_results = flux_variability_analysis(model,all_reactions,)
        all_results = pd.merge(fluxes, fva_results, left_index=True, right_index=True, how='outer')

        kwargs = {'mode':'a', 'if_sheet_exists':'replace'} \
            if os.path.exists(simulated_fluxes_result_file) else {'mode':'w'}
        with pd.ExcelWriter(simulated_fluxes_result_file,
                            engine='openpyxl', **kwargs) as writer:
            pd.DataFrame.from_dict(protein_concentrations, orient='index',
                                   columns = ['concentration'
                                              ]).reset_index().rename(
                {'index':'protein'}, axis = 1).to_excel(writer, sheet_name=label + '_protein', index = False)
            all_results.to_excel(writer, sheet_name=label)

def merge_flux_dfs_and_calc_cv(model_flux_dict:Dict[str,pd.DataFrame],
                           flux_col="flux",
                           id_col="rxn_id",
                               cog_col:str = "COG description",):
    """
    Merge the predicted fluxes of several models (given as a dict) and
    compute the coefficient‑of‑variation (CV) of those fluxes for each reaction.

    Parameters
    ----------
    model_flux_dict : dict[str, pd.DataFrame]
        ``{model_id: df}`` – each ``df`` must contain at least ``id_col`` and
        ``flux_col`` (the *predicted* flux, **not** the maximum/minimum values).
    flux_col : str, optional
        Name of the column that holds the predicted flux in each frame.
    id_col   : str, optional
        Column that identifies the reaction (common to all frames).
    cv_name  : str, optional
        Column name for the resulting CV.

    Returns
    -------
    pd.DataFrame
        A data‑frame with ``id_col`` and a column ``cv_name`` that holds the
        CV of the fluxes across the supplied models.  All model‑specific flux
        columns are also kept (prefixed with the model id) for downstream use.
    """

    wide = None
    for model_id, df in model_flux_dict.items():
        tmp = df.copy()
        tmp = tmp.rename(columns={flux_col: f"flux_{model_id}"})
        if wide is None:
            wide = tmp[[id_col, f"flux_{model_id}", cog_col]]
        else:
            wide = pd.merge(wide, tmp[[id_col, f"flux_{model_id}"]], on=id_col, how="outer").drop_duplicates()

    return determine_coefficient_of_variation_for_entity(df = wide, entity = flux_col)

def determine_flux_coefficient_of_variation():
    def calculate_cv(row):
        mean = (row['maximum'] + row['minimum'])/2
        variation = row['maximum'] - row['minimum']
        if mean == 0 or variation == 0: return 0
        return variation/mean

    fva_results = pd.read_excel(RESULT_FLUX_PATH, sheet_name=None)
    parsed_fva_results ={}
    for model_id, df in fva_results.items():
        if model_id.endswith('_protein'): continue
        # if value is smaller than feasibility tolerance, the value can be interpreted as 0
        df.loc[:, ["minimum", "maximum"]] = df[["minimum", "maximum"]].mask(
            df[["minimum", "maximum"]].abs() < 1e-6, 0
        )
        df['flux_cv'] = df.apply(calculate_cv, axis=1)
        parsed_fva_results[model_id] = df
    return parsed_fva_results

def determine_coefficient_of_variation_for_entity(df:pd.DataFrame, entity:str = 'kcat'):
    def calculate_cv(row):
        numeric_row = pd.to_numeric(row, errors='coerce').dropna()
        mean, std = numeric_row.mean(),numeric_row.std()
        if any([mean == 0, std==0]): return 0
        return std/mean
    df[f'{entity}_cv'] = df.apply(calculate_cv, axis=1)
    return df


def get_variability_statistics(df_flux,
                         df_kcat, kcat_cv,
                         df_prot, protein_cv,
                         cog_list=None):
        """
        Convenience wrapper – creates a figure with two rows per COG:
            row 1 – flux‑vs‑kcat density
            row 2 – protein‑vs‑kcat density
        """

        def get_stats(df, col_start='flux', cog=None):
            df = df[df["COG description"] == cog] if cog is not None else df
            columns = [c for c in df.columns
                       if c.startswith(col_start) and
                       not any([n in c for n in ['GotEnzymes', 'After preprocessing', 'cv']])
                       ]
            # Convert to a NumPy array – this flattens the 2‑D block into 1‑D
            values = df[columns].to_numpy().ravel()
            # Drop NaNs (if any) before the statistics
            values = abs(values[~np.isnan(values)])
            mean = values.mean()
            col_start = col_start if col_start != 'concentration' else 'protein'
            mean_cv = df[col_start + '_cv'].to_numpy().ravel().mean()
            std_cv = df[col_start + '_cv'].to_numpy().ravel().std()
            std = df[columns].std().mean()
            return mean, std, mean_cv, std_cv

        def get_pearson(df_x, df_y, cog=None):
            data = _merge_on_rxn(df_y, df_x, "_x")
            data = data[data["COG description"] == cog] if cog is not None else data
            cv_col_y = [col for col in df_y.columns if col.endswith("_cv")][0]
            cv_col_x = [col for col in df_x.columns if col.endswith("_cv")][0]
            r, p = pearsonr(x=data[cv_col_x],
                            y=data[cv_col_y])
            return r, p

        def get_stat_row(cog=None):
            mean_flux, std_flux, mean_flux_cv, std_flux_cv = get_stats(df_flux, cog=cog)
            mean_kcat, std_kcat, mean_kcat_cv, std_kcat_cv = get_stats(df_kcat, col_start='kcat', cog=cog)
            mean_protein, std_protein, mean_protein_cv, std_protein_cv = get_stats(df_prot, col_start='concentration', cog=cog)
            r_f_k, p_f_k = get_pearson(df_flux, kcat_cv, cog=cog)
            r_p_k, p_p_k = get_pearson(protein_cv, kcat_cv, cog=cog)
            cog = 'Overall' if cog is None else cog
            return {'flux_mean': mean_flux, 'flux_std': std_flux, 'flux_std_norm': std_flux / mean_flux,
                    'mean_flux_cv': mean_flux_cv, 'std_flux_cv': std_flux_cv,
                    'kcat_mean': mean_kcat, 'kcat_std': std_kcat, 'kcat_std_norm': std_kcat / mean_kcat,
                    'mean_kcat_cv': mean_kcat_cv, 'std_kcat_cv': std_kcat_cv,
                    'protein_mean': mean_protein, 'protein_std': std_protein,
                    'protein_std_norm': std_protein / mean_protein, 'mean_protein_cv': mean_protein_cv,
                    'std_protein_cv': std_protein_cv,
                    'pearson_flux_kcat': r_f_k, 'pval_pearson_flux_kcat': p_f_k,
                    'pearson_prot_kcat': r_p_k, 'pval_pearson_prot_kcat': p_p_k,
                    'cog': cog,
                    }

        if cog_list is None:
            cog_list = (
                df_flux["COG description"]
                # .dropna()
                .unique()
            )
        stat_rows = [get_stat_row()]
        for i, cog in enumerate(cog_list):
            try:
                stat_rows.append(get_stat_row(cog))
            except:
                continue

        return pd.DataFrame(stat_rows).dropna()


def annotate_df_with_cog(df:pd.DataFrame,
                         rxn2protein_mapping: pd.DataFrame,
                         ) -> pd.DataFrame:
    if 'protein' in df.columns:
        df = pd.merge(df, rxn2protein_mapping[['rxn_id', 'enzyme_id']],
                      left_on='protein', right_on='enzyme_id',
                      how='inner')

    if 'Unnamed: 0' in df.columns:
        df = df.rename({'Unnamed: 0': 'rxn_id'}, axis=1)
    return add_pathway_annotation_to_proteins(
        df_with_proteins=df,
        merge_on='rxn_id')

def classify_enzyme_cv_for_all_models(flux_cv,
                                      protein_cv,
                                      kcat_cv_no_cog):
    def enzyme_type_flux_kcat(flux_cv, kcat_cv):
        low_flux = flux_cv <= 0.5
        low_kcat = kcat_cv <= 0.5
        if low_flux and low_kcat: return '(i) constrained coupled'
        if low_kcat: return '(iii) kcat constrained uncoupled'
        if low_flux: return '(ii) flux constrained uncoupled'
        return '(iv) unconstrained uncoupled'

    def enzyme_type_protein_kcat(protein_cv, kcat_cv):
        low_protein = protein_cv <= 0.5
        low_kcat = kcat_cv <= 0.5
        if low_protein and low_kcat: return '(i) stable coupled'
        if low_kcat: return '(iii) stable kcat uncoupled'
        if low_protein: return '(ii) stable protein uncoupled'
        return '(iv) variable uncoupled'

    all_cvs = pd.merge(
        flux_cv[~flux_cv.filter(regex=r"^flux").eq(0).all(axis=1)][['rxn_id', 'flux_cv', 'COG description']],
        protein_cv[['rxn_id', 'enzyme_id', 'protein_cv']],
        on='rxn_id', how='right')
    all_cvs = pd.merge(all_cvs, kcat_cv_no_cog.reset_index()[['rxn_id', 'enzyme_id', 'direction', 'kcat_cv']],
                       on=['rxn_id', 'enzyme_id'])

    all_cvs = all_cvs[all_cvs['flux_cv'] > 0]
    all_cvs['enzyme_type_flux_kcat'] = all_cvs.apply(lambda row: enzyme_type_flux_kcat(row.flux_cv, row.kcat_cv), axis=1)
    all_cvs['enzyme_type_protein_kcat'] = all_cvs.apply(lambda row: enzyme_type_protein_kcat(row.protein_cv, row.kcat_cv), axis=1)
    return all_cvs

def determine_perc_per_type(flux_cv,
                            protein_cv,
                            kcat_cv_no_cog,
                            enzyme_type:Literal['flux_kcat', 'protein_kcat'] = 'flux_kcat'
                            ):
    enzyme_type = 'enzyme_type_' + enzyme_type
    all_cvs = classify_enzyme_cv_for_all_models(flux_cv=flux_cv,
                                                protein_cv=protein_cv,
                                                kcat_cv_no_cog=kcat_cv_no_cog)
    counts = (all_cvs.groupby(['COG description', enzyme_type])
              .size()
              .unstack(fill_value=0))
    perc = counts.div(counts.sum(axis=1), axis=0) * 100
    perc = perc.round(1).reset_index()

    overall_counts = (all_cvs
                      .groupby(enzyme_type)
                      .size())  # keep same order

    overall_perc = (overall_counts
                    / overall_counts.sum() * 100)
    overall_row = overall_perc.to_frame().T
    overall_row.insert(0, 'COG description', 'Overall')
    overall_row = overall_row[perc.columns]
    perc = pd.concat([perc, overall_row], ignore_index=True)
    print(perc.to_latex())
    return perc

def get_cvs_for_boxplot(df_flux, df_kcat, df_prot):
    all_zero_reactions = df_flux[df_flux.filter(regex=r"^flux").eq(0).all(axis=1)]['rxn_id']

    df_flux_box = df_flux[['rxn_id','COG description', 'flux_cv']].rename(
        columns={'flux_cv': 'cv'}).assign(metric='flux')


    df_kcat_box = df_kcat[['rxn_id','COG description', 'kcat_cv']] \
        .rename(columns={'kcat_cv': 'cv'}).assign(metric='kcat')

    df_prot_box = df_prot[['rxn_id','COG description', 'protein_cv']].rename(
        columns={'protein_cv': 'cv'}).assign(metric='protein')

    df_overall = pd.concat([df_flux_box, df_kcat_box, df_prot_box],
                       ignore_index=True)
    df_overall['COG description'] = 'Overall'

    all_cvs_per_cog = pd.concat([df_flux_box, df_kcat_box, df_prot_box, df_overall],
                       ignore_index=True)
    return all_cvs_per_cog[~all_cvs_per_cog['rxn_id'].isin(all_zero_reactions)]

def create_barplot_classified_enzymes_and_flux_mean(df_flux,
                     df_kcat, kcat_cv,
                     df_prot, protein_cv,kcat_cv_no_cog):
    stat_df = get_variability_statistics(df_flux = df_flux,
                                         df_kcat = df_kcat, kcat_cv = kcat_cv,
                                         df_prot = df_prot, protein_cv = protein_cv)
    perc_flux_kcat = determine_perc_per_type(flux_cv, protein_cv, kcat_cv_no_cog)
    perc_protein_kcat = determine_perc_per_type(flux_cv, protein_cv, kcat_cv_no_cog, enzyme_type='protein_kcat')

    df_boxplot = get_cvs_for_boxplot(df_flux, df_kcat, df_prot)

    cog_to_plot = ['Overall'] + [cog for cog in stat_df.sort_values('flux_mean', ascending=False)['cog'].iloc[:11] if
                                 cog != 'Overall']
    stat_df = stat_df[stat_df.cog.isin(cog_to_plot)].set_index('cog')
    stat_df = stat_df.loc[cog_to_plot].reset_index()

    df_boxplot = df_boxplot[df_boxplot['COG description'].isin(cog_to_plot)].set_index('COG description')
    df_boxplot = df_boxplot.loc[cog_to_plot].reset_index()

    plt.rcParams.update({'font.size': 11})
    fig = plt.figure(figsize=(21 / 2.54, 30 / 2.54))
    gs = gridspec.GridSpec(2, 1, height_ratios=[1,1.5], hspace=0.85)

    colors_flux_kcat = sns.color_palette("colorblind", n_colors=len(perc_flux_kcat.columns))
    colors_protein_kcat = sns.color_palette("Set2", n_colors=len(perc_flux_kcat.columns))
    print(colors_protein_kcat[1])

    color_map = {**{'flux': 'darkgrey', 'protein': 'slategrey'},
                 **{l: c for l, c in
                                                          zip([col
                                                               for col in list(perc_flux_kcat.columns)
                                                               if 'COG' not in col],
                                                              colors_flux_kcat)},

                 **{l: c for l, c in
                    zip([col
                         for col in list(perc_protein_kcat.columns)
                         if 'COG' not in col],
                        [colors_protein_kcat[2],colors_protein_kcat[3],colors_protein_kcat[0],colors_protein_kcat[1]])}
    }

    n_cogs = len(cog_to_plot)
    x_pos = np.arange(n_cogs)  # centre of each COG group
    group_width = 0.8  # total width occupied by a COG
    stack_width = 0.2  # width of the *stacked* bar
    ax = fig.add_subplot(gs[0])

    for perc, xposition in zip([perc_flux_kcat, perc_protein_kcat],
                               [x_pos - group_width / 2 + stack_width / 2, x_pos]):
        perc = perc[perc['COG description'].isin(cog_to_plot)].set_index('COG description')
        perc = perc.loc[cog_to_plot].reset_index()
        bottom = np.zeros(len(perc))
        for col in perc.columns:
            if col == 'COG description': continue
            heights = perc[col].values
            p = ax.bar(xposition,
                       heights,
                       width=stack_width,
                       label=col,
                       color=color_map[col],
                       bottom=bottom)
            bottom += heights

    # for stat_col, xposition in zip(['flux', 'protein'], [x_pos, x_pos + group_width / 2 - stack_width / 2]):
    for stat_col, xposition in zip(['flux'], [x_pos + group_width / 2 - stack_width / 2]):

        ax2 = ax.twinx()
        conv_factor = 1 if stat_col == 'flux' else 1e5
        ax2.errorbar(
            xposition,  # centre of the group
            stat_df[f"{stat_col}_mean"] *conv_factor,
            yerr=stat_df[f"{stat_col}_std"]*conv_factor,
            fmt='o',  # marker only (no line)
            color='black' if stat_col == 'flux' else 'midnightblue',
            capsize=4,
            zorder=5
        )
        ax2.bar(
            xposition,
            stat_df[f"{stat_col}_mean"]*conv_factor,
            width=stack_width,
            color=color_map[stat_col],
            label=f'Mean {stat_col}',
        )
        unit = r"rate [mmol/g$_{CDW}$/h]" if stat_col == 'flux' else r"concentration \n [$\cdot 10^{-5}$ g_${protein}$/g$_{CDW}$]"
        ax2.set_ylabel(f'Mean {stat_col} {unit}')
        stat_df['ymax'] = stat_df[f"{stat_col}_mean"] + stat_df[f"{stat_col}_std"]
        ax2.set_ylim([0, stat_df["ymax"].max() * 1.1 *conv_factor])
        if stat_col == 'protein': ax2.spines['right'].set_position(('outward', 60))

    ax.grid(visible=True, alpha=0.2, linewidth=0.7)
    ax.set_ylabel('Percentage of flux-carrying enzymes')
    ax.set_xticks(x_pos)
    ax.set_xticklabels([COG_MAPPER[c] for c in cog_to_plot], rotation=45, ha='right',)

    handles, labels = [], []
    for ax in fig.axes:
        h, l = ax.get_legend_handles_labels()
        for handle, label in zip(h, l):
            if label in labels: continue
            handles.append(handle)
            labels.append(label)

    #add dummy legend handles and labels to fill columns
    handles += [Line2D([],[], color = 'white')]*2
    labels += ['', '']

    ax.legend(handles, labels, bbox_to_anchor=(0.5, -0.97), loc='lower center', borderaxespad=0., ncols=3)
    gs_violin = gridspec.GridSpecFromSubplotSpec(3,1,subplot_spec=gs[1], hspace=0)
    color_map = {'kcat': 'purple', 'flux': 'blue', 'protein': 'red'}
    data_per_cog = {}
    for  cog, sub in df_boxplot.groupby('COG description'):
        data_per_cog[cog] = {
            'kcat': sub.loc[sub.metric == 'kcat', 'cv'].values,
            'flux': sub.loc[sub.metric == 'flux', 'cv'].values,
            'protein': sub.loc[sub.metric == 'protein', 'cv'].values
        }

    axs = [fig.add_subplot(gs_violin[i]) for i in range(3)]
    for type, ax_violin in zip(['kcat', 'flux', 'protein'], axs):
        data = [abs(sub_df[type]) for cog, sub_df in data_per_cog.items()]
        violin = ax_violin.violinplot(data,
                                      showmeans=True)
        ax_violin.set_ylabel(type)
        ax_violin.grid(visible=True, alpha=0.2, linewidth=0.7)
        ax_violin.set_ylim([-0.2, 3.9])

        for partname in ('cbars', 'cmins', 'cmaxes', 'cmeans'):
            vp = violin[partname]
            vp.set_edgecolor(color_map[type])
            vp.set_linewidth(1)

        for pc in violin['bodies']:
            pc.set_facecolor(color_map[type])
            pc.set_edgecolor(color_map[type])

    axs[1].set_ylabel(r'Coefficient of variation [$\frac{stdev}{mean}$]'+'\n' + r'k$_{\text{cat}}$')
    for ax in axs[:-1]:
        ax.set_xticks([])
        ax.set_xticklabels([])
    axs[-1].set_xticks([pos+1 for pos in x_pos])
    axs[-1].set_xticklabels([COG_MAPPER[c] for c in cog_to_plot], rotation=45, ha='right', )

    for ax, annotation in zip([fig.axes[0], fig.axes[3]], ['A', 'B']):
        ax.annotate(annotation, xy=(-0.1, 0.95), xycoords="axes fraction",
                    fontsize=16, fontweight='bold',
                    xytext=(-6, 4.5), textcoords="offset points",
                    ha="right", va="bottom")

    plt.tight_layout()
    plt.subplots_adjust(top=0.95, bottom=0.15)
    plt.savefig('Figures/SuppFig_flux_and_kcat_variability.png')


if __name__ == '__main__':
    NUM_ALT_MODELS = 10
    diagnostic_files = list()
    param_files = list()
    all_model_labels = ['iML1515','GotEnzymes', 'After preprocessing'] \
                                                 + [f'alternative {i}' for i in range(1, NUM_ALT_MODELS + 1)]
    for file_nmbr in range(1, NUM_ALT_MODELS + 1):
        diagnostic_files += [os.path.join('Results', '2_parametrization', 'diagnostics',
                                     f'pam_parametrizer_diagnostics_{file_nmbr}.xlsx')]
        param_files += [os.path.join('Results', '3_analysis', 'parameter_files',
                                     f'proteinAllocationModel_EnzymaticData_iML1515_{file_nmbr}.xlsx')]


    # get_flux_predictions_and_fva(model_files=[MODEL_FILE, PARAM_FILE_OLD,
    #                                   PARAM_FILE_PREPROC], diagnostic_files=diagnostic_files,
    #                                  label_names=all_model_labels)


    protein_concs = get_all_protein_concentrations()
    kcat_values = get_all_kcat_values(param_files,
                                      label_names= [f'alternative {i}' for i in range(1, NUM_ALT_MODELS + 1)]
                                      ).set_index(['rxn_id','enzyme_id', 'direction'])

    flux_cv_per_model_dict = determine_flux_coefficient_of_variation()

    kcat_cv_no_cog = determine_coefficient_of_variation_for_entity(kcat_values)
    protein_cv = determine_coefficient_of_variation_for_entity(protein_concs, entity='protein')
    rxn2protein_mapping = kcat_cv_no_cog.reset_index()[['rxn_id', 'enzyme_id', 'direction']].copy()

    kcat_cv = annotate_df_with_cog(kcat_cv_no_cog, rxn2protein_mapping)
    protein_cv = annotate_df_with_cog(protein_cv, rxn2protein_mapping)
    flux_cv_per_model_dict = {
        model_id: annotate_df_with_cog(df, rxn2protein_mapping)
        for model_id, df in flux_cv_per_model_dict.items()
    }
    flux_cv = merge_flux_dfs_and_calc_cv(flux_cv_per_model_dict)

    kcat_values =annotate_df_with_cog(kcat_values, rxn2protein_mapping)
    protein_concs =annotate_df_with_cog(protein_concs, rxn2protein_mapping)



    create_barplot_classified_enzymes_and_flux_mean(df_flux=flux_cv,
                     df_kcat=kcat_values,
                     df_prot=protein_concs,
                     kcat_cv=kcat_cv,
                     protein_cv=protein_cv,
                     kcat_cv_no_cog = kcat_cv_no_cog)







