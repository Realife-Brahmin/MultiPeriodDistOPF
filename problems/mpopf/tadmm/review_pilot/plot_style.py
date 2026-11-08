"""Match the manuscript's IEEE serif typography and PGFPlots colors.

Sources in the IAS-Trans manuscript: main.tex (IEEEtran),
figures/traj_common.tex and figures/scaling_inflation_common.tex.
The paper assigns violet to tADMM, green to feasible centralized references,
and maroon to infeasible centralized iterates. System colors label systems.
All diagnostic curves here are tADMM; distinguish runs with line styles.
"""
TADMM = '#6A4C9C'
CENTRAL_REFERENCE = '#2E8B57'
INFEASIBLE_REFERENCE = '#8C1D1D'
THRESHOLD = '#505050'
SYSTEM_COLORS = {'ieee123C_1ph': '#1F5186', 'ieee2522C_1ph': '#B8860B',
                 'large10kC_1ph': '#227822'}
OBJECTIVE = PRIMAL = DUAL = RHO = TADMM


def configure_paper_style():
    import matplotlib as mpl
    mpl.rcParams.update({
        'font.family': 'serif',
        'font.serif': ['Times New Roman', 'Nimbus Roman', 'STIXGeneral'],
        'mathtext.fontset': 'stix',
        'font.size': 11,
        'axes.labelsize': 11,
        'axes.titlesize': 12,
        'axes.titleweight': 'bold',
        'xtick.labelsize': 10,
        'ytick.labelsize': 10,
        'legend.fontsize': 9,
        'axes.linewidth': 0.6,
        'xtick.direction': 'out',
        'ytick.direction': 'out',
        'lines.linewidth': 1.5,
        'grid.color': '#E0E0E0',
        'grid.linewidth': 0.5,
        'grid.alpha': 1.0,
        'legend.frameon': True,
        'legend.facecolor': 'white',
        'legend.edgecolor': '#BFBFBF',
        'legend.framealpha': 0.92,
        'pdf.fonttype': 42,
        'ps.fonttype': 42,
        'savefig.dpi': 220,
    })
