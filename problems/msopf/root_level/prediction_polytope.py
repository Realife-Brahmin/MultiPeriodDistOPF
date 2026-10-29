"""Exact affine half-space prediction meshes; no radial clipping or artificial sphere."""
import itertools
import numpy as np

def prediction_mesh(P,K,case):
    P=np.asarray(P,dtype=float);K=np.asarray(K,dtype=float)
    if case=='cyclic':A=-K.mean(axis=1);b=P.mean(axis=1);sources=np.arange(P.shape[0])
    else:A=-K.reshape(-1,3);b=P.flatten();sources=np.repeat(np.arange(P.shape[0]),P.shape[1])
    norms=np.linalg.norm(A,axis=1);keep=norms>1e-12
    assert np.all(b[~keep]>=0),'Empty affine prediction'
    A=A[keep]/norms[keep,None];b=b[keep]/norms[keep];sources=sources[keep]
    triples=np.array(list(itertools.combinations(range(len(b)),3)),dtype=int)
    matrices=A[triples];valid=np.abs(np.linalg.det(matrices))>1e-10
    v=np.linalg.solve(matrices[valid],b[triples[valid]][...,None])[...,0]
    chunks=[]
    for part in np.array_split(v,max(1,int(np.ceil(len(v)/4096)))):
        if len(part):chunks.append(part[np.max(A@part.T-b[:,None],axis=0)<=1e-6])
    verts=np.unique(np.round(np.vstack(chunks),8),axis=0)
    assert len(verts)>=4,'Affine prediction is empty/degenerate or unbounded; do not fabricate a cap'
    I=[];J=[];L=[];faces=[];seen=set();plane_errors=[]
    for plane,(n,c,source) in enumerate(zip(A,b,sources)):
        ids=np.flatnonzero(np.abs(verts@n-c)<1e-5)
        key=tuple(ids)
        if len(ids)<3 or key in seen:continue
        seen.add(key)
        ref=np.array([0.,0.,1.]) if abs(n[2])<.9 else np.array([1.,0.,0.])
        e1=np.cross(n,ref);e1/=np.linalg.norm(e1);e2=np.cross(n,e1)
        q=verts[ids]-verts[ids].mean(axis=0);ids=ids[np.argsort(np.arctan2(q@e2,q@e1))]
        for j in range(1,len(ids)-1):
            tri=verts[ids[[0,j,j+1]]]
            if np.linalg.norm(np.cross(tri[1]-tri[0],tri[2]-tri[0]))<1e-8:continue
            I.append(int(ids[0]));J.append(int(ids[j]));L.append(int(ids[j+1]));faces.append(int(source))
            plane_errors.append(float(np.max(np.abs(tri@n-c))))
    assert len(I)>=4
    assert max(plane_errors)<1e-5
    assert np.max(A@verts.T-b[:,None])<1e-5
    return dict(x=verts[:,0].tolist(),y=verts[:,1].tolist(),z=verts[:,2].tolist(),i=I,j=J,k=L,face=faces,
                construction='exact affine half-space intersection; no radius cap',max_plane_residual_deg=max(plane_errors))

def attach_prediction_meshes(data,root):
    """Recover the exact archived sensitivities for every displayed dataset."""
    import csv,hashlib
    from pathlib import Path
    root=Path(root)
    for s in data:
        key=s['id'];names=s['sources'];candidates=[]
        if key.startswith('screen_'):
            label={'screen_sf4':'sf4','screen_8500xf4':'secondary4'}[key]
            candidates=[root/'results/secondary_8500_oct25/quick_T6'/label/f'sensitivity_{label}.csv']
        else:candidates=sorted((root/'results').glob('**/multi_sens_'+key+'_T24.csv'))
        for p in candidates:
            with p.open(encoding='utf-8-sig',newline='') as f:d=list(csv.DictReader(f))
            source='source' if 'source' in d[0] else 'substation'
            prefix='K_';suffix='' if source=='source' else '_kW_per_deg'
            P=np.array([[float(x['import_kW']) for x in d if x[source]==k] for k in names])
            K=np.array([[[float(x[prefix+m+suffix]) for m in names[1:]] for x in d if x[source]==k] for k in names])
            if P.shape[0]!=4 or not np.allclose(-K.sum(axis=1)/P.sum(axis=1)[:,None],np.array(s['normals'])/np.array(s['offsets'])[:,None],rtol=1e-4,atol=1e-7):continue
            for case in ('nobatt','cyclic'):
                s['cases'][case]['mesh']=prediction_mesh(P,K,case)
                s['cases'][case].pop('surface',None)
            s['prediction_source']=dict(path=str(p.relative_to(root)),sha256=hashlib.sha256(p.read_bytes()).hexdigest())
            if key.startswith('screen_'):
                s['sparse']=True
                s['study_note']=s['study_note'].replace('Predicted surfaces are radially clipped for display at 2.5 times each case maximum measured reach; this cap is not a feasibility constraint.','Predictions are exact affine half-space meshes, with no radius cap. Sparse samples do not establish a measured surface.')
            break
        else:raise ValueError('No matching archived sensitivity input for '+key)
    return data
