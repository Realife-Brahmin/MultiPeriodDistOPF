"""Mirror the existing tADMM convergence plot's semantic color palette.

Source: ../plots/plot_convergence_from_csv.jl, Wong (2011) palette block.
Keep metric colors stable; distinguish systems/backends/parity by line or marker.
"""
OBJECTIVE = '#0072B2'       # RGB(0, 114, 178)
PRIMAL = '#009E73'          # RGB(0, 158, 115)
DUAL = '#E69F00'            # RGB(230, 159, 0)
RHO = '#CC79A7'             # RGB(204, 121, 167)
CENTRAL_REFERENCE = '#D55E00'  # RGB(213, 94, 0)
THRESHOLD = '#505050'       # RGB(80, 80, 80)
