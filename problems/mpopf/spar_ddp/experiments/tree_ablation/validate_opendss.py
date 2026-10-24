import re,sys
text=open(sys.argv[1]).read()
checks=[s for s in text.splitlines() if s.startswith('OPENDSS_CHECK ')]
voltage=[s for s in text.splitlines() if s.startswith('OPENDSS_SUMMARY ') and '| voltage:' in s]
if len(checks)!=1 or not re.search(r'\bconverged=true\b', checks[0]):
    raise SystemExit('OpenDSS did not report converged=true')
if len(voltage)!=1: raise SystemExit('Missing or duplicate OpenDSS voltage summary')
low=re.search(r'below limit: (\d+) of',voltage[0]); high=re.search(r'above limit: (\d+),',voltage[0])
if low is None or high is None: raise SystemExit('Unknown OpenDSS voltage summary schema')
if int(low[1]) or int(high[1]): raise SystemExit('OpenDSS voltage violations')
print('TREE_ABLATION_OPENDSS_PASS '+voltage[0])
