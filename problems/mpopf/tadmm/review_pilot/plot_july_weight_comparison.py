"""Compare weight changes using identical inputs and each weight's own reference."""
import csv
import hashlib
import json
from pathlib import Path
import sys
import tomllib
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from plot_style import TADMM, THRESHOLD, SYSTEM_COLORS, configure_paper_style

configure_paper_style()
root=Path(__file__).resolve().parent/'results'
out=Path(sys.argv[1]);out.parent.mkdir(parents=True,exist_ok=True)
scenarios=[('ieee2522C_1ph','med2522',
    root/'ddp_transfer_NO/balance2/results/ieee2522C_1ph_T24_rho3800_relative_mu10_NO_31650298',
    root/'accuracy_fix_T24_31645746/results/ieee2522C_1ph_T24_ma57_firstpass_31645746_a0',
    root/'july_weight_comparison/med2522/results/ieee2522C_1ph_T24_july_pair_31650588'),
    ('large10kC_1ph','large10k',
    root/'ddp_transfer_NO/results/large10k_T24_rho48_fixed_NO_31650238',
    root/'accuracy_fix_T24_31645790/results/large10kC_1ph_T24_ma57_firstpass_31645790_a0',
    root/'july_weight_comparison/large10k/results/large10kC_1ph_T24_july_pair_31650589')]
fig,axes=plt.subplots(2,2,figsize=(10,7),sharex='col',layout='constrained')
inputs=[]
for column,(system,name,old,ref,new) in enumerate(scenarios):
    om=tomllib.loads((old/'calibration.toml').read_text())
    rm=tomllib.loads((ref/'metadata.toml').read_text())
    ov=json.loads((old/'calibration_validation.json').read_text())
    nm=tomllib.loads((new/'comparison.toml').read_text())
    nv=json.loads((new/'validation_complete.json').read_text())
    assert (ref/'COMPLETE').is_file() and (new/'ASSESSMENT_VALIDATED').is_file()
    assert ov['stage_status_validation']==ov['source_validation']=='pass'
    assert nv['source_validation']=='pass'
    assert om['reference_summary_sha256']==hashlib.sha256((ref/'summary.toml').read_bytes()).hexdigest()
    assert om['reference_input_sha256']==nm['input_sha256']==rm['data_sha256']
    for key in ('rho_initial','rho_policy','screening'):
        assert om[key]==nm[key]
    assert om.get('rho_balance_ratio',10)==nm['rho_balance_ratio']
    for directory,filename,reference,cb,style,expected in [
        (old,'calibration.csv',om['reference_objective'],rm['C_B'],'--',om['no_iteration']),
        (new/'exact','trace.csv',nm['reference_objective'],nm['C_B'],'-',nv['arms']['exact']['no_iteration'])]:
        file=directory/filename;rows=list(csv.DictReader(file.open()))
        x=[int(r['iteration']) for r in rows]
        gap=[max(abs(float(r['objective'])-reference)/max(1,abs(reference))/.001,1e-12) for r in rows]
        physical=[max(float(r['original_primal'])/1e-5,1e-12) for r in rows]
        first=next((int(r['iteration']) for r,g,p in zip(rows,gap,physical) if g<=1 and p<=1),0)
        assert first==expected
        status=f'NO {first}' if first else f'NO not reached by {len(rows)}'
        label=f'$C_B={cb:.4g}$; {status}'
        axes[0,column].semilogy(x,gap,color=TADMM,linestyle=style,label=label)
        axes[1,column].semilogy(x,physical,color=TADMM,linestyle=style)
        if first:
            for axis in axes[:,column]:
                axis.axvline(first,color=THRESHOLD,linestyle=style,alpha=.35)
        inputs.append(hashlib.sha256(file.read_bytes()).hexdigest()+'  '+file.relative_to(root).as_posix()+'\n')
    axes[0,column].set_title(name,color=SYSTEM_COLORS[system])
    axes[0,column].legend(fontsize=9)
    axes[1,column].set_xlabel('ADMM iteration')
for axis in axes.flat:
    axis.axhline(1,color=THRESHOLD,linestyle=':',linewidth=.8)
    axis.axvline(100,color=THRESHOLD,linestyle=':',linewidth=.8)
    axis.grid(which='major');axis.set_xlim(1,200)
axes[0,0].set_ylabel('Objective gap / NO tolerance')
axes[1,0].set_ylabel('Original violation / NO tolerance')
fig.suptitle('tADMM, $T=24$, MA57: changing only the battery weight\nSame periodic inputs and rho policy per system; objective gaps use their own matched reference',fontsize=12)
for ext in ('png','pdf'):
    fig.savefig(out.with_suffix('.'+ext))
plt.close(fig)
for path in (Path(__file__),Path(__file__).with_name('plot_style.py')):
    inputs.append(hashlib.sha256(path.read_bytes()).hexdigest()+'  '+path.name+'\n')
out.with_suffix('.sha256').write_text(''.join(inputs))
