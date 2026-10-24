import ast,csv,hashlib,json,pathlib,re,sys
def digest(path):
    h=hashlib.sha256()
    with path.open('rb') as f:
        for chunk in iter(lambda:f.read(8*1024**2),b''):h.update(chunk)
    return h.hexdigest()
C=pathlib.Path(sys.argv[1]);system=sys.argv[2];T=int(sys.argv[3]);P=C/'results'/f'{system}_T{T}_dense_capture';capture=P/'captures'
with (capture/'coverage.csv').open() as f:r=list(csv.DictReader(f))
by={x['tag']:x for x in r};assert len(r)==len(by)==5
assert set(by)=={'early','middle','late','barrier_before','barrier_after'}
assert int(by['early']['iteration'])==0 and int(by['early']['barrier_iteration'])==0
assert int(by['early']['pass'])<int(by['middle']['pass'])<int(by['late']['pass'])
a,b=by['barrier_before'],by['barrier_after']
assert a['iteration']==b['iteration'] and int(b['barrier_iteration'])==int(a['barrier_iteration'])+1
assert int(b['pass'])==int(a['pass'])+1 and float(b['mu'])<float(a['mu'])
assert int(by['late']['barrier_iteration'])>0
text=(P/'solve.log').read_text();assert re.search(r'solve complete: [\d.]+ s, iterations=\d+, status=0',text)
match=re.search(r'solve complete: [\d.]+ s, iterations=(\d+)',text);assert int(by['late']['iteration'])==int(match[1])
for name in {x['snapshot'] for x in r}:
    D=capture/name;textmeta=(D/'metadata.toml').read_text()
    assert int(re.search(r'(?m)^stages = (\d+)',textmeta)[1])==T
    assert digest(D/'leaves.jls')==re.search(r'(?m)^leaves_sha256 = "([0-9a-f]+)"',textmeta)[1]
    # This writer emits simple scalar TOML and [[stage_files]] tables; avoid
    # a Python3.11-only parser on Kamiak's Python3.9.
    records=textmeta.split('[[stage_files]]')[1:];assert len(records)==len(set((1,(T+1)//2,T)))
    for record in records:
        file=re.search(r'(?m)^file = "([^"]+)"',record)[1];sha=re.search(r'(?m)^sha256 = "([0-9a-f]+)"',record)[1]
        assert digest(D/file)==sha
    tags=ast.literal_eval(re.search(r'(?m)^tags = (\[[^\n]+\])',textmeta)[1])
    for tag in tags:assert by[tag]['snapshot']==name
assert all('FROZEN_VALUE_CAPTURE' in text and f'snapshot={x["snapshot"]}' in text for x in r)
print(json.dumps(dict(system=system,T=T,actual_tags=sorted(by),coverage_verified=True,points=r),indent=2))
