"""Publish only independently replayed OPF evidence; never infer unresolved boundaries."""
import csv,json,math,itertools
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
ROOT=Path(__file__).resolve().parents[1]
OUT=ROOT/'results/secondary_8500_oct25'
def rows(p):
    with p.open(newline='',encoding='utf-8-sig') as f:return list(csv.DictReader(f))
def validate(folder,complete=True):
    meta=(folder/'metadata.txt').read_text(encoding='utf-8')
    if complete:assert 'exit_code=0' in meta,folder
    r=rows(folder/'solves.csv');ok=[];failed=[]
    for x in r:
        if x['status'] not in ('LOCALLY_SOLVED','OPTIMAL','ALMOST_LOCALLY_SOLVED'):
            failed.append(x['status']);continue
        assert x['replay_ok']=='true',(folder,x)
        for key,lim in [('dP_kW',1e-3),('dQ_kvar',1e-3),('dV_pu',1e-7),('dB_kW',1e-3),('terminal_kWh',.01),('soc_violation_kWh',.01),('energy_balance_kWh',.01),('battery_power_violation_kW',.01),('pv_q_violation_kvar',.01)]:
            assert math.isfinite(float(x[key])) and float(x[key])<lim,(folder,key,x[key])
        ok.append(x)
    e=rows(folder/'edges.csv') if (folder/'edges.csv').exists() else []
    for x in e:
        assert float(x['vmin'])>=.95-1e-5 and float(x['vmax'])<=1.05+1e-5 and float(x['min_import_kW'])>=-.01
    return {'folder':str(folder.relative_to(ROOT)),'verified_solves':len(ok),'unsuccessful_statuses':failed,'unresolved_rays':sum(x['method']=='unresolved' for x in e),'max_dP_kW':max(float(x['dP_kW']) for x in ok),'max_dV_pu':max(float(x['dV_pu']) for x in ok)}
def sensitivity(folder,n=4):
    d=rows(next(folder.glob('sensitivity_*.csv')))
    P=np.array([[float(x['import_kW']) for x in d if x['source']==f'subs{i}'] for i in range(1,n+1)])
    K=np.array([[[float(x[f'K_subs{j}']) for j in range(2,n+1)] for x in d if x['source']==f'subs{i}'] for i in range(1,n+1)])
    return P,K

def append_viewer(folder,label,title):
    P,K=sensitivity(folder);base=rows(folder/'baseline.csv');edges=rows(folder/'edges.csv')
    sys={'id':label,'title':title,'sources':[f'subs{i}' for i in range(1,5)],'normals':(-K.sum(axis=1)).tolist(),'offsets':P.sum(axis=1).tolist(),'linearity_max':None,'study_note':f'Exploratory T = {P.shape[1]} over 24 hours; six axis directions per feasible battery case, not a reconstructed surface. Dashed-model predictions omit voltage constraints. PV reactive dispatch; hard cyclic SOC. No random linearity test.','cases':{}}
    sys['study_note'] += ' Equal-angle local feasibility: ' + '; '.join(b['case'] + ' ' + ('admissible' if b['origin_feasible']=='true' else 'negative local margin; no rays attempted') for b in base) + '. Endpoint categories use import/voltage proximity heuristics, not active-set certificates. Backflow source identities were not archived. Predictions are exact affine half-space meshes; no radial clipping is used. Sparse samples do not establish a measured surface.'
    for case in ('nobatt','cyclic'):
        points=[]
        for e in edges:
            if e['case']!=case or e['method'] not in ('bracketed','coarse_bracket'):continue
            u=np.array([float(e[f'u{i}']) for i in (1,2,3)])
            A=K@u
            candidates=np.where(A.sum(axis=1)<0,P.sum(axis=1)/np.maximum(-A.sum(axis=1),1e-300),np.inf) if case=='cyclic' else np.min(np.where(A<0,P/np.maximum(-A,1e-300),np.inf),axis=1)
            face=int(np.argmin(candidates));r=float(e['measured_deg'])
            if e['binding']=='backflow':
                # Determine the binding source from the recorded endpoint solve, not the predictor.
                # Per-source endpoint imports are not archived in this screen: mark identity unknown.
                bound=-2
            else:bound=-1 if e['binding'] in ('undervoltage','overvoltage') else -2
            points.append([*(u*r).tolist(),face,r,float(e['predicted_deg']),False,bound,e['binding']+' (source identity not archived)' if bound==-2 else e['binding'],float(e['upper_deg'])])
        sys['cases'][case]={'points':points}
    return sys

def main():
    reports=[]
    import shutil,hashlib
    shutil.copy2(ROOT/'Manifest.toml',OUT/'Manifest.snapshot.toml')
    decks={}
    for f in sorted((ROOT.parents[1]/'test_networks/smartds_p1u_xf_1ph').glob('*.dss')):
        decks[f.name]=hashlib.sha256(f.read_bytes()).hexdigest()
    (OUT/'sf_deck_hashes_postrun.json').write_text(json.dumps(decks,indent=2)+'\n')
    for parent in ('baseline_T24','full_T24','quick_T6'):
        for f in sorted((OUT/parent).glob('*')):
            if f.is_dir() and (f/'metadata.txt').exists() and 'exit_code=0' in (f/'metadata.txt').read_text(encoding='utf-8'):reports.append(validate(f))
    (OUT/'validated_publication.json').write_text(json.dumps(reports,indent=2)+'\n')
    plt.rcParams.update({'font.family':'DejaVu Sans','font.size':7,'pdf.fonttype':42})
    fig,ax=plt.subplots(figsize=(3.45,2.75))
    folder=OUT/'full_T24/secondary3';P,K=sensitivity(folder,3)
    for case,c in [('nobatt','#2866a6'),('cyclic','#d45c32')]:
        es=[e for e in rows(folder/'edges.csv') if e['case']==case]
        pts=[]
        for e in es:
            if e['method']=='unresolved':pts.append([np.nan,np.nan]);continue
            u=np.array([float(e['u1']),float(e['u2'])]);r=float(e['measured_deg']);hi=float(e['upper_deg']);q=u*r
            pts.append(q);ax.plot(*np.array([q,u*hi]).T,color=c,alpha=.7,lw=1.5)
            ax.scatter(*q,color=c,s=25,marker='o' if e['binding']=='backflow' else 's',zorder=3)
        pts=np.array(pts);pts=np.vstack([pts,pts[0]])
        ax.plot(*pts.T,color=c,lw=1,label='Idle batteries' if case=='nobatt' else 'Cyclic batteries')
        pred=[]
        for t in np.linspace(0,2*np.pi,721):
            u=np.array([np.cos(t),np.sin(t)]);A=K@u
            rr=np.min(np.where(A.sum(axis=1)<0,P.sum(axis=1)/np.maximum(-A.sum(axis=1),1e-300),np.inf)) if case=='cyclic' else np.min(np.where(A<0,P/np.maximum(-A,1e-300),np.inf))
            pred.append(u*rr)
        ax.plot(*np.array(pred).T,color=c,ls='--',alpha=.55,lw=.9)
    ax.set(xlim=(-23,23),ylim=(-18,23),xlabel=r'$\delta_2-\delta_1$ (deg)',ylabel=r'$\delta_3-\delta_1$ (deg)')
    ax.axhline(0,color='#ddd',lw=.7);ax.axvline(0,color='#ddd',lw=.7);ax.grid(alpha=.15)
    ax.legend(loc='upper left',frameon=False,fontsize=6.5)
    fig.tight_layout()
    for ext in ('png','pdf','svg'):fig.savefig(OUT/f'three_source_T24.{ext}',dpi=180)
    plt.close(fig)
    sf=OUT/'quick_T6/sf4';validate(sf)
    fig,ax=plt.subplots(figsize=(3.45,2.65))
    es=rows(sf/'edges.csv');labs=[r'$-\delta_2$',r'$+\delta_2$',r'$-\delta_3$',r'$+\delta_3$',r'$-\delta_4$',r'$+\delta_4$']
    for case,c,dy in [('nobatt','#2866a6',-.12),('cyclic','#d45c32',.12)]:
        ee=[e for e in es if e['case']==case]
        ax.scatter([float(e['measured_deg']) for e in ee],np.arange(6)+dy,label='Idle' if case=='nobatt' else 'Cyclic',color=c,s=22)
        ax.scatter([float(e['predicted_deg']) for e in ee],np.arange(6)+dy,marker='|',color=c,s=75)
    ax.set_yticks(range(6),labs);ax.set_xlabel('Reach (deg); dots measured, ticks predicted');ax.set_title('San Francisco P1U Type III\nExploratory T = 6 axis screen',fontsize=8)
    ax.grid(axis='x',alpha=.2);ax.legend(frameon=False,ncol=2,loc='lower right');fig.tight_layout()
    for ext in ('png','pdf','svg'):fig.savefig(OUT/f'sf_axis_T6.{ext}',dpi=180)
    plt.close(fig)
    data=json.loads((ROOT/'results/multiperiod/four_viewer_data_fib78.json').read_text(encoding='utf-8'))
    data.append(append_viewer(sf,'screen_sf4','San Francisco P1U - T6 screen'))
    f=OUT/'quick_T6/secondary4'
    if (f/'metadata.txt').exists() and 'exit_code=0' in (f/'metadata.txt').read_text(encoding='utf-8'):
        validate(f);data.append(append_viewer(f,'screen_8500xf4','IEEE8500 + secondaries - T6 screen'))
    from prediction_polytope import attach_prediction_meshes
    data=attach_prediction_meshes(data,ROOT)
    encoded=json.dumps(data,separators=(',',':'),allow_nan=False)
    (OUT/'four_viewer_data.json').write_text(encoded,encoding='utf-8')
    v=ROOT/'results/multiperiod/viewer'
    html=(v/'viewer_head.html').read_text(encoding='utf-8')+encoded+(v/'viewer_tail.html').read_text(encoding='utf-8')
    (v/'four_regions.html').write_text(html,encoding='utf-8')
    (ROOT.parents[1]/'docs/four_regions.html').write_text(html,encoding='utf-8')
    print(json.dumps({'validated_cases':len(reports),'verified_solves':sum(r['verified_solves'] for r in reports),'viewer_systems':len(data)},indent=2))
if __name__=='__main__':main()
