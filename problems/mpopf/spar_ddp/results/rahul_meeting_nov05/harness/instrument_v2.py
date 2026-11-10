from pathlib import Path
import re
import sys

source, replay = (Path(p) for p in sys.argv[1:])
text = source.read_text()
operations = {
    "K .=": "assemble_matrix", "g .=": "assemble_vector",
    "C = cholesky!": "cholesky", "transpose!": "transpose", "copyto!": "copy",
    "ldiv!(U', g)": "solve_vector", "ldiv!(U', w.": "solve_matrix",
    "BLAS.syrk!": "syrk", "LinearAlgebra.copytri!": "mirror",
    "mul!(x.Mab": "gemm", "mul!(x.ma": "gemv", "mul!(x.mb": "gemv",
}
names = list(dict.fromkeys(operations.values())) + ["fallback"]
functions = re.findall(r"^function _pp_combine!\(x::PintElem, y::Pint(?:Elem|Val), w::PairWork\)\n.*?^end$", text, re.M | re.S)
assert len(functions) == 2
generated = []
for function in functions:
    lines = function.splitlines()
    generated.append(lines[0].replace("_pp_combine!", "profile_combine!").replace("w::PairWork)", "w::PairWork, ticks)"))
    for line in lines[1:]:
        if "r = pint_combine" in line:
            generated.append("        stamp = time_ns()")
        if "return " in line and line.startswith("        "):
            generated.append(f'        ticks[{names.index("fallback")+1}] += time_ns() - stamp')
        if line.strip().startswith("C = cholesky!"):
            generated.extend(("    stamp = time_ns()", line,
                              f"    ticks[{names.index('cholesky')+1}] += time_ns() - stamp"))
            continue
        if line.startswith("    ") and not line.startswith("        "):
            statements = line.strip().split("; ")
            if all(any(stmt.startswith(prefix) for prefix in operations) for stmt in statements):
                for stmt in statements:
                    name = next(name for prefix, name in operations.items() if stmt.startswith(prefix))
                    generated.extend(("    stamp = time_ns()", "    " + stmt,
                                      f"    ticks[{names.index(name)+1}] += time_ns() - stamp"))
                continue
        generated.append(line)
Path("merge_timers.jl").write_text("const MERGE_OPERATIONS = " + repr(",".join(names)).replace("'", '"') + "\n" + "\n".join(generated) + "\n")
s = replay.read_text()
s = s.replace('    const REPLAY_INPUT', '    include(joinpath(ENV["REPLAY_ROOT"], "merge_timers.jl"))\n    const REPLAY_INPUT', 1)
assert s.count('        z=_pp_combine!(x,y,REPLAY_WORK[])') == 1
s = s.replace('        z=_pp_combine!(x,y,REPLAY_WORK[])',
              '        ticks = zeros(UInt64, length(split(MERGE_OPERATIONS, char(44))))\n        z=profile_combine!(x,y,REPLAY_WORK[],ticks)'.replace('char(44)', "','"))
s = s.replace('pid=myid())', 'pid=myid(),ticks=ticks,two_sided=!val)', 1)
s = s.replace('counts=[32,48,63]', 'counts=[63]')
needle = '            println("ROUND_DIAG'
assert needle in s
i = s.index(needle)
s = s[:i] + '''            for (q,r) in enumerate(rr)
                println("MERGE_OP,$W,$scheduling,$rep,$round,$q,$(r.pid),$(r.two_sided),$(r.compute),$(join(r.ticks,','))")
            end
''' + s[i:]
s = s.replace('counts=[63]', 'println("MERGE_OP_FIELDS,workers,mode,rep,round,task,pid,two_sided,kernel_s,", MERGE_OPERATIONS)\ncounts=[63]')
Path("replay.jl").write_text(s)
