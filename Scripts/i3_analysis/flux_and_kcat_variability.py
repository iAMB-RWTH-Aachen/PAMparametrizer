from cobra.io.sbml import read_sbml_model
import pandas as pd
import os
import numpy as np
from matplotlib import pyplot as plt
from matplotlib.lines import Line2D
from matplotlib import gridspec
import seaborn as sns
from typing import Union, List, Dict
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

def plot_flux_cv_by_model(flux_dfs, model_names, ax=None):
    """
    Parameters
    ----------
    flux_dfs : list[pd.DataFrame]
        Each element contains a column named ``flux_cv``.
    model_names : list[str]
        Human‑readable names for the models (same order as `flux_dfs`).
    ax : matplotlib.axes.Axes, optional
        Axes onto which the plot is drawn.

    Returns
    -------
    matplotlib.axes.Axes
    """
    # 1️⃣  Assemble long‑form table
    long = pd.concat(
        [
            df.assign(model=name)[["flux_cv", "model"]]
            for df, name in zip(flux_dfs, model_names)
        ],
        ignore_index=True,
    )
    long['log_cv'] = long.apply(lambda row: np.log10(abs(row['flux_cv'])) if row['flux_cv']!=0 else 0, axis=1)
    # 2️⃣  Plot
    if ax is None:
        fig, ax = plt.subplots(figsize=(max(8, len(model_names) * 1.2), 6))
    sns.violinplot(
        data=long.dropna(),
        x="model",
        y="log_cv",
        inner="quartile",
        cut=0,
        # log_scale=True,
        ax=ax,
    )
    ax.set_yscale("log")
    ax.set_ylabel("Absolute flux CV (log scale)")
    ax.set_xlabel("")
    plt.setp(ax.get_xticklabels(), rotation=45, ha="right")
    return ax


def plot_overall_cv(df_flux, df_kcat, df_prot, ax=None):
    """
    Combine the three CV columns into one long table and plot as violins.
    """
    # 1️⃣  Build long table
    dfs = [
        df_flux.rename(columns={"flux_cv": "cv"}).assign(type="flux"),
        df_kcat.rename(columns={"kcat_cv": "cv"}).assign(type="kcat"),
        df_prot.rename(columns={"protein_cv": "cv"}).assign(type="protein"),
    ]
    long = pd.concat(dfs, ignore_index=True)[["cv", "type"]]

    # 2️⃣  Plot
    if ax is None:
        fig, ax = plt.subplots(figsize=(6, 6))

    ax.boxplot(long.dropna())
    ax.set_yscale("log")
    ax.set_ylabel("Absolute CV (log scale)")
    ax.set_xlabel("")
    return ax

# ----------------------------------------------------------------------
# 3.  COG‑specific density / scatter plots
# ----------------------------------------------------------------------
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


def plot_cog_density_distributions(df_y: pd.DataFrame,
                                  df_x: pd.DataFrame,
                                  cog:str = None,
                                  gs=None,
                                   fig = None
                                  ):
    """
    Density plot (hexbin) of flux CV (y) vs kcat CV (x) for a single COG.
    """
    data = _merge_on_rxn(df_y, df_x, "_x")
    data = data[data["COG description"] == cog] if cog is not None else data

    cv_col_y = [col for col in df_y.columns if col.endswith("_cv")][0]
    cv_col_x = [col for col in df_x.columns if col.endswith("_cv")][0]
    if gs is None:
        fig = plt.figure(figsize=(6, 5))
        gs = GridSpec(ncols=1, nrows=1, wspace=0.05, hspace=0.05)

    fig, (ax_joint, ax_margx, ax_margy) = kde_joint_plot(x = data[[cv_col_x]],
                       y = data[[cv_col_y]], fig=fig, gs =gs,
                   xlog=True, ylog=True,)


    xlabel = f"{cv_col_x.replace('_cv', '')}"
    ylabel = f"{cv_col_y.replace('_cv', '')}"

    cog = cog if cog is not None else ''
    ax_joint.set_xlabel(f"{xlabel} CV (log)"+cog)
    ax_joint.set_ylabel(f"{ylabel} CV (log)")
    return ax_joint

def kde_joint_plot(
    x,
    y,
    *,
    bw_method="scott",
    gridsize=100,
    fill_color="#4c72b0",
        density_without_color = 1e-2,
        cmap_name: str ='viridis',
    fill_alpha=0.4,
    line_color="#4c72b0",
    line_width=1.5,
    marginal_color="#4c72b0",
    marginal_alpha=0.8,
    marginal_linewidth=1.5,
    xlog=False,
    ylog=False,
    ax=None,
    figsize=(6, 6),
        fig=None,
        gs =None,
    title=None,
):
    """
    Joint KDE (filled contour) with marginal 1‑D KDEs on the top and right side.

    Parameters
    ----------
    x, y : array‑like, shape (n_samples,)
        Data to be plotted.
    bw_method : str or float, optional
        Band‑width selector passed to `gaussian_kde`.  Default: ``'scott'``.
    gridsize : int, optional
        Number of grid points per dimension for the 2‑D density.
    fill_color, line_color, marginal_color : colour specifications
    fill_alpha, marginal_alpha : float
        Transparency for the filled contour and marginal curves.
    xlog, ylog : bool, optional
        Plot the axes on a symlog scale (useful when data contain zeros or
        negative values).  Set to ``True`` to apply the transformation.
    ax : matplotlib.axes.Axes, optional
        Central axis to embed the plot into.  If ``None`` a new figure is
        created.
    figsize : tuple, optional
        Figure size when a new figure is created.
    title : str, optional
        Plot title.

    Returns
    -------
    fig, (ax_joint, ax_marg_x, ax_marg_y) : tuple
        Figure and the three axes (joint, top‑margin, right‑margin).
    """
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)

    mask = np.isfinite(x) & np.isfinite(y)
    x, y = x[mask], y[mask]

    if xlog:
        x = np.sign(x) * np.log10(np.abs(x) + 1e-12)
    if ylog:
        y = np.sign(y) * np.log10(np.abs(y) + 1e-12)

    if fig is None and ax is None:
        fig = plt.figure(figsize=figsize)
        gs = GridSpec(2, 2, figure=fig,
                      width_ratios=[4, 1],
                      height_ratios=[1, 4],
                      wspace=0.05,
                      hspace=0.05)
    if fig is not None and gs is not None:
        gs = GridSpecFromSubplotSpec(2, 2, subplot_spec=gs,
                                     width_ratios=[4, 1],
                                     height_ratios=[1, 4],
                                     wspace=0.05,
                                     hspace=0.05
                                     )
        ax_joint   = fig.add_subplot(gs[1, 0])
        ax_marg_x = fig.add_subplot(gs[0, 0],sharex=ax_joint)
        ax_marg_y = fig.add_subplot(gs[1, 1], sharey=ax_joint)
    elif ax is not None:
        # embed into already‑existing axes (central one)
        fig = ax.figure
        # create two small axes aligned with the supplied one
        pos = ax.get_position()
        width, height = pos.width, pos.height
        pad = 0.02  # relative padding between panels

        ax_marg_x = fig.add_axes(
            [pos.x0, pos.y0 + height + pad, width, height * 0.25], sharex=ax
        )
        ax_marg_y = fig.add_axes(
            [pos.x0 + width + pad, pos.y0, width * 0.25, height], sharey=ax
        )
        ax_joint = ax

    # ------------------------------------------------------------------
    # 2️⃣  2‑D KDE (joint)
    # ------------------------------------------------------------------
    kde2d = gaussian_kde(np.vstack([x, y]), bw_method=bw_method)
    xmin, xmax = x.min(), x.max()
    ymin, ymax = y.min(), y.max()
    # add a small margin so the contour does not touch the frame
    xr = xmax - xmin
    yr = ymax - ymin
    r = 0.05
    xmin -= r * xr
    xmax += r * xr
    ymin -= r * yr
    ymax += r * yr

    xx, yy = np.meshgrid(
        np.linspace(xmin, xmax, gridsize),
        np.linspace(ymin, ymax, gridsize)
    )
    zz = kde2d(np.vstack([xx.ravel(), yy.ravel()])).reshape(xx.shape)
    cmap = plt.get_cmap(cmap_name).copy()  # copy so we don't modify the global cmap
    cmap.set_under("white")  # colour used for values < vmin

    # density threshold that will be treated as “very low”
    vmin = density_without_color * zz.max()
    norm = plt.Normalize(vmin=vmin, vmax=zz.max())

    # filled contour – values below vmin become white
    ax_joint.contourf(
        xx, yy, zz,
        levels=20,
        cmap=cmap,
        norm=norm,
        alpha=fill_alpha,
    )
    # density outline
    ax_joint.contour(
        xx, yy, zz, levels=10, colors=line_color, linewidths=line_width
    )

    # ------------------------------------------------------------------
    # 3️⃣  1‑D KDEs (marginals)
    # ------------------------------------------------------------------
    kde_x = gaussian_kde(x, bw_method=bw_method)
    kde_y = gaussian_kde(y, bw_method=bw_method)

    xv = np.linspace(xmin, xmax, gridsize * 2)
    yv = np.linspace(ymin, ymax, gridsize * 2)

    ax_marg_x.plot(xv, kde_x(xv), color=marginal_color,
                   lw=marginal_linewidth, alpha=marginal_alpha)
    ax_marg_x.fill_between(xv, 0, kde_x(xv), color=marginal_color,
                           alpha=marginal_alpha)

    ax_marg_y.plot(kde_y(yv), yv, color=marginal_color,
                   lw=marginal_linewidth, alpha=marginal_alpha)
    ax_marg_y.fill_betweenx(yv, 0, kde_y(yv), color=marginal_color,
                            alpha=marginal_alpha)

    ax_marg_x.axis("off")
    ax_marg_y.axis("off")

    if xlog:
        ax_joint.set_xlabel("log10(|x|)"); ax_marg_x.set_xlabel("")
    else:
        ax_joint.set_xlabel("x")
    if ylog:
        ax_joint.set_ylabel("log10(|y|)"); ax_marg_y.set_ylabel("")
    else:
        ax_joint.set_ylabel("y")

    if title:
        ax_joint.set_title(title, pad=20)

    return fig, (ax_joint, ax_marg_x, ax_marg_y)

def plot_all_cogs_density(df_flux, df_kcat, df_prot, cog_list=None):
    """
    Convenience wrapper – creates a figure with two rows per COG:
        row 1 – flux‑vs‑kcat density
        row 2 – protein‑vs‑kcat density
    """
    def flux_stats(df, flux_columns=None):
        if flux_columns is None:
            flux_columns = [c for c in df.columns if c.startswith("flux")]
        # Convert to a NumPy array – this flattens the 2‑D block into 1‑D
        values = df[flux_columns].to_numpy().ravel()
        # Drop NaNs (if any) before the statistics
        values = abs(values[~np.isnan(values)])
        mean = values.mean()
        std = values.std(ddof=1)  # sample std, same convention as pandas .std()
        return mean, std


    if cog_list is None:
        cog_list = (
            df_flux["COG description"]
            .dropna()
            .unique()
        )
    n = len(cog_list)

    fig = plt.figure(figsize=(12, 4 * n),
        constrained_layout=True,)

    gs = GridSpec(ncols=2, nrows=n+1, wspace=0.1, hspace=0.1)

    plot_cog_density_distributions(df_flux, df_kcat, gs=gs[0, 0], fig=fig)
    plot_cog_density_distributions(df_prot, df_kcat,gs=gs[0, 1], fig=fig)

    mean_flux, std_flux = flux_stats(flux_cv)
    print(f'Flux distribution: mean: {mean_flux}, std: {std_flux}\n')
    for i, cog in enumerate(cog_list):
        print(cog)
        mean_flux, std_flux = flux_stats(flux_cv[flux_cv["COG description"]==cog])
        print(f'Flux distribution: mean: {mean_flux}, std: {std_flux}')
        try:
            plot_cog_density_distributions(df_flux, df_kcat, cog, gs = gs[i+1,0], fig=fig)
            plot_cog_density_distributions(df_prot, df_kcat, cog, gs = gs[i+1,1], fig=fig)
            y_pos = 1 - (i + 0.5) / n
            fig.text(
                0.5, y_pos,
                f"COG: {cog}",
                ha="center",
                va="center",
                fontsize=12,
                weight="bold",
            )

        except:
            continue
    return fig

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


def plot_variability(df_flux,
                     df_kcat, kcat_cv,
                     df_prot, protein_cv,
                     cog_list=None):
    """
    Convenience wrapper – creates a figure with two rows per COG:
        row 1 – flux‑vs‑kcat density
        row 2 – protein‑vs‑kcat density
    """

    if cog_list is None:
        cog_list = list(
            df_flux["COG description"]
            # .dropna()
            .unique()
        )
        cog_list.append("Overall")
    stat_df = get_variability_statistics(df_flux = df_flux,
                                         df_kcat = df_kcat, kcat_cv = kcat_cv,
                                         df_prot = df_prot, protein_cv = protein_cv,
                                         cog_list = cog_list)
    print(stat_df.to_markdown())
    size_range = (30, 350)
    plt.rcParams['font.size'] = 12
    fig = plt.Figure(figsize=(21/2.54, 17/2.54), layout='tight')
    gs = plt.GridSpec(ncols=2, nrows=2, height_ratios=[2,1])
    axs = [fig.add_subplot(gs_inner) for gs_inner in [gs[0,0], gs[0,1]]]
    # axs = [fig.add_subplot(gs_inner) for gs_inner in [gs[0,0], gs[0,1], gs[1,0],  gs[1,1],]]


    colors = sns.color_palette("Set2", n_colors=len(cog_list))
    cmap = {l: c for l, c in
           zip(cog_list, colors)}
    legend_elements = {}
    labelmapper = {'_std_norm': ' CV', '_mean':' mean', 'pearson_': 'pearson correlation coefficient ', '_cv':' CV'}

    for (x,y, size), ax in zip([("flux_std_norm", "kcat_std_norm", "flux_mean"),
                                ("kcat_std_norm", "protein_std_norm", "protein_mean")],
                               axs):
        for cog, sub_df in stat_df.groupby('cog'):
            label = f"{COG_DESCRIPTION2LETTER[cog]}: {COG_MAPPER[cog]}" if cog is not None else 'all'
            ax.scatter(
                x=sub_df[x],
                y=sub_df[y],
                s=np.interp(np.abs(sub_df[size]),
                            [stat_df[size].min(), stat_df[size].max()],
                            size_range),
                color=cmap[cog],
                marker="o",
                clip_on=False,
                label = label
                )
            ax.annotate(label.split(':')[0], (sub_df[x].iloc[0], sub_df[y].iloc[0]))
            if cog not in legend_elements:
                legend_elements[cog] = Line2D([0],[0], color='w', marker='o', clip_on=False, label = label, markerfacecolor=cmap[cog])
        xlabel, ylabel, bubble_size = x, y, size
        for old, new in labelmapper.items(): xlabel = xlabel.replace(old, new)
        for old, new in labelmapper.items(): ylabel = ylabel.replace(old, new)
        for old, new in labelmapper.items(): bubble_size = bubble_size.replace(old, new)


        ax.set_xlabel(xlabel.replace('_', ' vs. '))
        # if not 'pearson' in x:
        #     ax.set_xscale("log")
        # if not 'pearson' in y:
        #     ax.set_yscale("log")
        ax.set_ylabel(ylabel.replace('_', ' vs. '))
        ax.set_title(f"bubble size:{bubble_size}", fontsize =11)
        ax.grid(visible=True, alpha=0.2, linewidth=0.7)
    legend_ax = fig.add_subplot(gs[1,:])
    legend_ax.legend(handles=legend_elements.values(), loc='upper left', bbox_to_anchor=(0,1.7), fontsize = 11, ncols = 2)
    legend_ax.set_axis_off()
    fig.savefig('Results/3_analysis/bubble_chart.png')
    plt.show()
    return fig


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
    def enzyme_type(flux_cv, kcat_cv):
        low_flux = flux_cv <= 0.5
        low_kcat = kcat_cv <= 0.5
        if low_flux and low_kcat: return '(i) constrained coupled'
        if low_kcat: return '(iii) kcat constrained uncoupled'
        if low_flux: return '(ii) flux constrained uncoupled'
        return '(iv) unconstrained uncoupled'

    all_cvs = pd.merge(
        flux_cv[~flux_cv.filter(regex=r"^flux").eq(0).all(axis=1)][['rxn_id', 'flux_cv', 'COG description']],
        protein_cv[['rxn_id', 'enzyme_id', 'protein_cv']],
        on='rxn_id', how='right')
    all_cvs = pd.merge(all_cvs, kcat_cv_no_cog.reset_index()[['rxn_id', 'enzyme_id', 'direction', 'kcat_cv']],
                       on=['rxn_id', 'enzyme_id'])

    all_cvs = all_cvs[all_cvs['flux_cv'] > 0]
    all_cvs['enzyme_type'] = all_cvs.apply(lambda row: enzyme_type(row.flux_cv, row.kcat_cv), axis=1)


    counts = (all_cvs.groupby(['COG description', 'enzyme_type'])
              .size()
              .unstack(fill_value=0))
    perc = counts.div(counts.sum(axis=1), axis=0) * 100
    perc = perc.round(1).reset_index()

    overall_counts = (all_cvs
                      .groupby('enzyme_type')
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
    perc = classify_enzyme_cv_for_all_models(flux_cv, protein_cv, kcat_cv_no_cog)
    df_boxplot = get_cvs_for_boxplot(df_flux, df_kcat, df_prot)

    cog_to_plot = ['Overall'] + [cog for cog in stat_df.sort_values('flux_mean', ascending=False)['cog'].iloc[:10] if
                                 cog != 'Overall']
    stat_df = stat_df[stat_df.cog.isin(cog_to_plot)].set_index('cog')
    stat_df = stat_df.loc[cog_to_plot].reset_index()

    perc = perc[perc['COG description'].isin(cog_to_plot)].set_index('COG description')
    perc = perc.loc[cog_to_plot].reset_index()

    df_boxplot = df_boxplot[df_boxplot['COG description'].isin(cog_to_plot)].set_index('COG description')
    df_boxplot = df_boxplot.loc[cog_to_plot].reset_index()

    plt.rcParams.update({'font.size': 11})
    fig = plt.figure(figsize=(21 / 2.54, 30 / 2.54))
    gs = gridspec.GridSpec(2, 1, height_ratios=[1,1.5], hspace=0.7)
    bottom = np.zeros(len(perc))

    colors = sns.color_palette("colorblind", n_colors=len(perc.columns))
    color_map = {**{'flux': 'darkgrey', 'protein': 'slategrey'}, **{l: c for l, c in
                                                          zip([col for col in perc.columns if 'COG' not in col],
                                                              colors)}}
    cog_labels = perc['COG description'].tolist()
    n_cogs = len(cog_labels)
    x_pos = np.arange(n_cogs)  # centre of each COG group
    group_width = 0.8  # total width occupied by a COG
    stack_width = 0.2  # width of the *stacked* bar

    ax = fig.add_subplot(gs[0])
    for col in perc.columns:
        if col == 'COG description': continue
        heights = perc[col].values
        p = ax.bar(x_pos - group_width / 2 + stack_width / 2,
                   heights,
                   width=stack_width,
                   label=col,
                   color=color_map[col],
                   bottom=bottom)
        bottom += heights

    for stat_col, xposition in zip(['flux', 'protein'], [x_pos, x_pos + group_width / 2 - stack_width / 2]):
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
    ax.set_xticklabels([COG_MAPPER[c] for c in cog_labels], rotation=45, ha='right',)

    handles, labels = [], []
    for ax in fig.axes:
        h, l = ax.get_legend_handles_labels()
        for handle, label in zip(h, l):
            if label in labels: continue
            handles.append(handle)
            labels.append(label)

    ax.legend(handles, labels, bbox_to_anchor=(0.5, -0.77), loc='lower center', borderaxespad=0., ncols=3)
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
    for type, ax_violin in zip(df_boxplot.metric.unique(), axs):
        data = [abs(sub_df[type]) for cog, sub_df in data_per_cog.items()]
        violin = ax_violin.violinplot(data,
                                      showmeans=True)
        ax_violin.set_ylabel(type)
        ax_violin.grid(visible=True, alpha=0.2, linewidth=0.7)


        for partname in ('cbars', 'cmins', 'cmaxes', 'cmeans'):
            vp = violin[partname]
            vp.set_edgecolor(color_map[type])
            vp.set_linewidth(1)

        for pc in violin['bodies']:
            pc.set_facecolor(color_map[type])
            pc.set_edgecolor(color_map[type])

    axs[1].set_ylabel(r'Coefficient of variation [$\frac{mean}{stdev}$]'+'\n' + r'k$_{\text{cat}}$')
    for ax in axs[:-1]:
        ax.set_xticks([])
        ax.set_xticklabels([])
    axs[-1].set_xticks([pos+1 for pos in x_pos])
    axs[-1].set_xticklabels([COG_MAPPER[c] for c in cog_labels], rotation=45, ha='right', )

    for ax, annotation in zip([fig.axes[0], fig.axes[3]], ['A', 'B']):
        ax.annotate(annotation, xy=(-0.1, 0.95), xycoords="axes fraction",
                    fontsize=16, fontweight='bold',
                    xytext=(-6, 4.5), textcoords="offset points",
                    ha="right", va="bottom")

    plt.tight_layout()
    plt.subplots_adjust(top=0.95, bottom=0.15)
    plt.savefig('Results/3_analysis/SuppFig_flux_and_kcat_variability.png')


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


    # ------------------------------------------------------------------
    # 1) Flux‑CV per model
    # ------------------------------------------------------------------
    # plot_flux_cv_by_model(flux_dfs=flux_cv_per_model_dict.values(),
    #                       model_names = flux_cv_per_model_dict.keys(), )
    # plt.tight_layout()
    # plt.show()

    # ------------------------------------------------------------------
    # 2) Overall CV distribution
    # ------------------------------------------------------------------
    # plot_overall_cv(flux_cv, kcat_cv, protein_cv)
    # plt.tight_layout()
    # plt.show()
    #
    # # ------------------------------------------------------------------
    # # 3) COG‑specific density / scatter plots
    # # ------------------------------------------------------------------
    plot_variability(df_flux=flux_cv,
                     df_kcat=kcat_values,
                     df_prot=protein_concs,
                     kcat_cv=kcat_cv,
                     protein_cv=protein_cv,)

    create_barplot_classified_enzymes_and_flux_mean(df_flux=flux_cv,
                     df_kcat=kcat_values,
                     df_prot=protein_concs,
                     kcat_cv=kcat_cv,
                     protein_cv=protein_cv,
                     kcat_cv_no_cog = kcat_cv_no_cog)







