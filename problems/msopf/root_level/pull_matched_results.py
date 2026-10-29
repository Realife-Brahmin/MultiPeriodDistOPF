"""Fetch checkpointed MSOPF evidence; never run simulations or modify remote results."""
from pathlib import Path
import subprocess,tarfile,copy,json
ROOT=Path(__file__).resolve().parents[1]/'results/matched_secondary_oct25'
REMOTE='aryan.jha@kamiak.wsu.edu'
command='cd ~/msopf/jobs/matched_secondary_oct25 && paths=(results logs submitted.tsv submitted_refs.tsv queue_before.txt tasks_*.txt); for p in oracle_results oracle_tasks.txt oracle_pilot_job.txt oracle_queue_before.txt submitted_oracle.tsv; do if [[ -e "$p" ]]; then paths+=("$p"); fi; done; tar -czf partial_evidence.tgz --exclude="*.jls" "${paths[@]}"'
r=subprocess.run(['ssh','-o','BatchMode=yes','-o','ConnectTimeout=12',REMOTE,command],capture_output=True,text=True)
if r.returncode and not (r.returncode==1 and 'file changed as we read it' in r.stderr):raise RuntimeError(r.stderr)
subprocess.run(['scp',REMOTE+':msopf/jobs/matched_secondary_oct25/partial_evidence.tgz',str(ROOT/'cluster_snapshot.tgz')],check=True)
dest=ROOT/'cluster';dest.mkdir(exist_ok=True);mapping={}
with tarfile.open(ROOT/'cluster_snapshot.tgz') as archive:
 for original in archive:
  member=copy.copy(original)
  if '\r' in member.name:member.name=member.name.replace('\r','__CR');mapping[original.name]=member.name
  archive.extract(member,dest,filter='data')
(ROOT/'cluster_path_mapping.json').write_text(json.dumps(mapping,indent=2)+'\n',encoding='utf-8')
print('Updated raw evidence snapshot; encoded original CR path components, leaving numerical file contents unchanged.')
