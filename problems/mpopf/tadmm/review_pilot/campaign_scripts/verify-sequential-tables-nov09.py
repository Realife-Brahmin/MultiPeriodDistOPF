"""Independent CSV re-scoring of clocks, continuity, powers and energy budgets."""
import argparse
import csv
import gzip
import hashlib
import json
import math
from pathlib import Path
import re
import tomllib

p=argparse.ArgumentParser();p.add_argument('root',type=Path);p.add_argument('output',type=Path);args=p.parse_args()
# Current frozen DSS file is uniform95%; explicitly verify instead of assuming
# the minimum efficiency represents mixed devices in older table schemas.
storage=(args.root/'metadata/campaign/validation_source/Storage.dss').read_text()
ecs={float(v)/100 for v in re.findall(r'%EffCharge=([\d.]+)',storage,re.I)}
eds={float(v)/100 for v in re.findall(r'%EffDischarge=([\d.]+)',storage,re.I)}
assert ecs==eds=={.95}
def read_table(case,name):
    path=case/(name+'.csv.gz')
    data=gzip.decompress(path.read_bytes())
    return data,list(csv.DictReader(data.decode().splitlines()))
def same(a,b,tol=1e-6):
    assert math.isclose(a,b,rel_tol=1e-10,abs_tol=tol),(a,b)
result=[]
for case in sorted((args.root/'tables').iterdir()):
    if not case.is_dir():continue
    s=tomllib.loads((case/'summary.toml').read_text());tables={}
    for name in ('battery','bus','branch','period'):
        data,rows=read_table(case,name)
        assert hashlib.sha256(data).hexdigest()==s['table_sha256'][name+'.csv']
        tables[name]=rows
    T=s['T'];dt=s['dt_h'];ec=s['used_charge_efficiency_min'];ed=s['used_discharge_efficiency_min']
    assert ec==ed==(1.0 if s['efficiency_mode']=='model' else .95)
    assert len(tables['period'])==T and len(tables['battery'])==T*s['batteries']
    before={};requested={};modelreq={};counts={};maxtrack=maxbalance=maxbound=maxreq=maxopf=0.0;clamps=0;worst=[]
    for row in tables['battery']:
        t=int(row['t']);bus=int(row['bus']);assert t==counts.get(bus,0)+1;counts[bus]=counts.get(bus,0)+1
        cmd=float(row['command_kw']);actual=float(row['actual_kw']);e0=float(row['energy_before_kwh']);e1=float(row['energy_after_kwh']);eo=float(row['opf_energy_kwh'])
        assert all(math.isfinite(x) for x in (cmd,actual,e0,e1,eo))
        same(e0,before.get(bus,e0));before[bus]=e1
        row_ec=float(row.get('charge_efficiency',ec));row_ed=float(row.get('discharge_efficiency',ed))
        change=dt*(row_ec*max(-actual,0)-max(actual,0)/row_ed)
        balance=e1-e0-change;same(balance,float(row['engine_energy_balance_error_kwh']))
        requested[bus]=requested.get(bus,e0)+dt*(row_ec*max(-cmd,0)-max(cmd,0)/row_ed)
        modelreq[bus]=modelreq.get(bus,e0)-dt*cmd
        same(requested[bus],float(row['requested_efficiency_energy_kwh']))
        same(modelreq[bus],float(row['requested_model_energy_kwh']))
        lo=float(row['SOC_low_kwh']);hi=float(row['SOC_high_kwh'])
        bound=max(0,lo-e1,e1-hi);reqbound=max(0,lo-requested[bus],requested[bus]-hi)
        same(bound,float(row['engine_bound_error_kwh']));same(reqbound,float(row['requested_bound_error_kwh']))
        maxtrack=max(maxtrack,abs(actual-cmd));maxbalance=max(maxbalance,abs(balance));maxbound=max(maxbound,bound);maxreq=max(maxreq,reqbound);maxopf=max(maxopf,abs(e1-eo))
        clamps+=abs(balance)>.01
        worst.append(dict(t=t,bus=bus,command_kw=cmd,actual_kw=actual,engine_energy_kwh=e1,opf_energy_kwh=eo,energy_accounting_error_kwh=balance,engine_bound_error_kwh=bound,requested_bound_error_kwh=reqbound))
    assert all(v==T for v in counts.values()) and len(counts)==s['batteries']
    for key,value in [('max_power_tracking_kw',maxtrack),('max_engine_energy_balance_kwh',maxbalance),('max_engine_SOC_bound_kwh',maxbound),('max_requested_SOC_bound_kwh',maxreq),('max_engine_vs_opf_energy_kwh',maxopf)]:same(value,s[key])
    assert clamps==s['reserve_or_capacity_clamp_records']
    for t,row in enumerate(tables['period'],1):
        assert int(row['t'])==t;same(float(row['hour_start']),(t-1)*dt);same(float(row['hour_end']),t*dt)
    maxdv=max(abs(float(r['dss_voltage_pu'])-float(r['opf_voltage_pu'])) for r in tables['bus']);same(maxdv,s['max_voltage_difference_pu'],1e-10)
    for key,cols in [('max_line_P_difference_kw',('dss_P_kw','opf_P_kw')),('max_line_Q_difference_kvar',('dss_Q_kvar','opf_Q_kvar')),('max_ell_difference_pu',('dss_ell_pu','opf_ell_pu'))]:
        same(max(abs(float(r[cols[0]])-float(r[cols[1]])) for r in tables['branch']),s[key])
    # All current cases use1000kVA base, independently recorded by the original
    # exported instances/campaign; infer/check here from the unit conversion.
    base=s.get('base_kva',1000.0)
    execution=s['converged'] and s['voltage_violations']==0 and max(maxtrack,maxbalance,maxbound,maxreq)/base<=1e-5
    samepass=execution and maxopf/base<=1e-5 and max(s[k] for k in ('max_substation_difference_kw','max_line_P_difference_kw','max_line_Q_difference_kvar'))/base<=1e-5 and s['max_voltage_difference_pu']<=1e-4
    assert execution==s['dispatch_energy_execution_pass'] and samepass==s['sequential_schedule_pass']
    worst.sort(key=lambda r:max(abs(r['energy_accounting_error_kwh']),r['requested_bound_error_kwh'],r['engine_bound_error_kwh']),reverse=True)
    top=case/'worst-battery-discrepancies.csv'
    with top.open('w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=list(worst[0]));w.writeheader();w.writerows(worst[:100])
    result.append(dict(case=case.name,passed=samepass,table_rows={k:len(v) for k,v in tables.items()},
        max_tracking_kw=maxtrack,max_energy_accounting_error_kwh=maxbalance,max_engine_SOC_bound_kwh=maxbound,max_unconstrained_requested_SOC_bound_kwh=maxreq,
        max_engine_vs_OPF_energy_kwh=maxopf,clamp_records=clamps,worst_battery=worst[0],all_table_hashes_match=True,
        independent_csv_checks_pass=True))
args.output.write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result,indent=2))
