"""Independently validate a square-cell MSPB endpoint using direct Fourier sums."""
from pathlib import Path
import json
import struct
import sys
import numpy as np
from scipy.spatial import cKDTree


def read_snapshot(path):
    data=Path(path).read_bytes()
    signature,version,dim,count,species,flags=struct.unpack_from('<8sIIQQI',data)
    assert (signature,version,dim,count,species,flags)==(b'MSPB\r\n\x1a\n',1,2,12000,3,3)
    assert len(data)==68+36*count
    basis=np.frombuffer(data,dtype='<f8',count=4,offset=36).reshape(2,2).copy()
    records=np.frombuffer(data,dtype=np.dtype([('id','<u8'),('species','<u4'),('relative','<f8',(2,)),('radius','<f8')]),count=count,offset=68).copy()
    records.sort(order='id')
    assert np.array_equal(records['id'],np.arange(count))
    assert np.array_equal(np.bincount(records['species']),[4000]*3)
    length=np.sqrt(float(count))
    assert np.allclose(basis,np.eye(2)*length,atol=1e-12,rtol=0)
    assert np.isfinite(records['relative']).all() and np.isfinite(records['radius']).all()
    assert np.allclose(records['radius'],np.sqrt(.1/np.pi),rtol=0,atol=1e-15)
    return dict(basis=basis,records=records,length=length)


def validate(folder):
    folder=Path(folder)
    snapshot=read_snapshot(folder/'endpoint.mspb')
    records=snapshot['records']
    length=snapshot['length']
    positions=(records['relative']%1)@snapshot['basis']
    modes=np.loadtxt(folder/'modes.tsv',skiprows=1)
    force=np.zeros_like(positions)
    cc=0.0
    species_max_s=[]
    counts=[]
    for species in range(3):
        selected=records['species']==species
        xyz=positions[selected]
        k=modes[modes[:,0]==species,1:3]
        weights=modes[modes[:,0]==species,3]
        assert np.all(weights==1)
        integer_modes=k*length/(2*np.pi)
        assert np.allclose(integer_modes,np.rint(integer_modes),rtol=0,atol=2e-13)
        assert len(np.unique(np.rint(integer_modes).astype(int),axis=0))==len(k)
        counts.append(len(k))
        local_force=np.zeros_like(xyz)
        max_s=0.0
        for start in range(0,len(k),128):
            block=k[start:start+128]
            phase=block[:,0,None]*xyz[None,:,0]+block[:,1,None]*xyz[None,:,1]
            co=np.cos(phase);si=np.sin(phase)
            re=co.sum(axis=1);im=si.sum(axis=1)
            rho2=re*re+im*im
            cc+=float(rho2.sum())/(3*len(positions))
            max_s=max(max_s,float((rho2/len(xyz)).max()))
            gradient_factor=(2/(3*len(positions)))*(re[:,None]*si-im[:,None]*co)
            for d in range(2):
                local_force[:,d]+=(gradient_factor*block[:,d,None]).sum(axis=0)
        force[selected]=local_force
        species_max_s.append(max_s)
    cc_force=force.copy()
    diameter=2*float(records['radius'][0])
    pairs=cKDTree(positions,boxsize=length).query_pairs(np.nextafter(diameter,np.inf),output_type='ndarray')
    contact=0.0;max_overlap=0.0
    if len(pairs):
        delta=positions[pairs[:,1]]-positions[pairs[:,0]]
        delta-=length*np.rint(delta/length)
        distance=np.linalg.norm(delta,axis=1)
        assert np.all(distance>1e-14)
        overlap=np.maximum(0,1-distance/diameter)
        max_overlap=float(overlap.max())
        contact=1*float(np.sum(overlap*overlap))
        pair_force=(-2/diameter*np.maximum(0,1/distance-1/diameter))[:,None]*delta
        np.add.at(force,pairs[:,0],pair_force)
        np.add.at(force,pairs[:,1],-pair_force)
    energy=cc+contact
    force_max=float(np.max(np.abs(force)))
    criteria=dict(energy_below_1e_minus14=energy<1e-14,
                  max_species_S_below_1e_minus10=max(species_max_s)<1e-10,
                  max_force_component_below_1e_minus7=force_max<1e-7,
                  max_relative_overlap_below_1e_minus6=max_overlap<1e-6)
    energy100=cc+100*contact
    force100_max=float(np.max(np.abs(cc_force+100*(force-cc_force))))
    physical100=dict(independent_energy=energy100,max_force_component=force100_max,
        ground_state=energy100<1e-14 and max(species_max_s)<1e-10 and force100_max<1e-7 and max_overlap<1e-6)
    raw=json.loads((folder/'stdout.log').read_text())
    assert abs(raw['direct_energy']-energy)<=2e-10*max(1,abs(energy))
    assert abs(raw['direct_force_max_abs']-force_max)<=2e-10*max(1,force_max)
    return dict(independent_energy=energy,cc_energy=cc,contact_energy=contact,
                species_max_S=species_max_s,max_S=max(species_max_s),
                mode_counts=counts,chi_actual=[n/(2*(4000-1)) for n in counts],
                max_force_component=force_max,max_relative_overlap=max_overlap,
                criteria=criteria,ground_state=all(criteria.values()),validation_at_contact_strength100=physical100)


if __name__=='__main__':
    result=validate(sys.argv[1])
    path=Path(sys.argv[1])/'physical-validation.json'
    path.write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(result))
