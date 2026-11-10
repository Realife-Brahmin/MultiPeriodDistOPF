"""Audit-only central logging; numerical DDP/Ipopt options stay inherited."""
from pathlib import Path
import sys
root=Path(sys.argv[1])
p=root/'spar_ddp/examples/power_system/centralized_ipopt_matched.jl'
p.write_text('''using LinearAlgebra
println("MATCHED_CENTRAL_RUNTIME blas=$(BLAS.get_num_threads()) backend=$(BLAS.get_config()) omp=$(get(ENV, "OMP_NUM_THREADS", "unset"))")
println("MATCHED_CENTRAL_AFFINITY ", read(`/usr/bin/taskset -apc $(getpid())`, String))
@assert BLAS.get_num_threads() == 1
'''+p.read_text())
s=Path('affinity_legacy.py').read_text()
assert s.count('64')==5, s.count('64')
Path('affinity.py').write_text(s.replace('64','budget').replace('workers=int(sys.argv[1])','budget=int(os.environ["STUDY_CORE_BUDGET"])\nworkers=int(sys.argv[1])'))

