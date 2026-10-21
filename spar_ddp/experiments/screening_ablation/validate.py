import re,sys,math
text=open(sys.argv[1]).read(); tol=1e-6 if sys.argv[2]=='ieee123C_1ph' else 1e-5
for bad in ('WARMUP_STREAM_MISMATCH','FILTERDDP_PROCS_FAIL','UNHANDLED TASK ERROR','ERROR: LoadError'):
    if bad in text: raise SystemExit('Rejected: '+bad)
assert 'SCREEN_ABLATION_ORIGINAL_BOUNDS enabled' in text
assert re.search(r'solve complete: [0-9.]+ s, iterations=\d+, status=0\b',text)
m=re.search(r'final residuals: primal=([\deE.+-]+) dual=([\deE.+-]+) complementarity=([\deE.+-]+)',text)
assert m and all(math.isfinite(float(v)) and 0<=float(v)<=tol for v in m.groups()), 'final residual gate'
b=re.search(r'BOUND_CHECK .*voltages below=(\d+) above=(\d+) worst_violation=([\deE.+-]+) min_ell=([\deE.+-]+) min_P_Subs=([\deE.+-]+)',text)
assert b and int(b[1])==int(b[2])==0 and float(b[3])<=tol and float(b[4])>=-tol and float(b[5])>=-tol, 'original bound gate'
assert 'FILTERDDP_ORIGINAL_PRIMAL' in text, 'missing per-iteration original feasibility'
print('SCREEN_ABLATION_NUMERICAL_PASS')
