#!/usr/bin/env python3
"""Export numerical sensitivity heatmaps from recorded metrics; no FEFF runs."""
from pathlib import Path
import json
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import Normalize

ROOT=Path(__file__).resolve().parent
report=json.loads((ROOT/'results.json').read_text())
windows=['3-8','8-12','12-16','16-20']
results=[r for r in report['results'] if not r['setting'].startswith('baseline')]
results.sort(key=lambda r:(r['reference'],r['setting']))
labels=[r['reference'].split('-COD')[0].replace('-rutile','')+'  '+r['setting'].replace('_cluster', ' / cluster ').replace('scf','SCF ') for r in results]
panels=[('complex_relative_rms',100.,'Complex NRMS (%)',5.),
        ('amplitude_relative_rms',100.,'Amplitude NRMS (%)',3.),
        ('phase_weighted_rms_rad',1.,'Phase RMS (rad)',.03)]
fig,axes=plt.subplots(1,3,figsize=(13, max(5,len(results)*.35)),sharey=True,constrained_layout=True)
for ax,(metric,scale,title,threshold) in zip(axes,panels):
    data=np.array([[max(p[w][metric] for key,p in r['paths'].items() if key!='complete_O_shell')*scale for w in windows] for r in results])
    vmax=max(float(data.max()),threshold)
    im=ax.imshow(data,cmap='YlOrRd',aspect='auto',norm=Normalize(vmin=0,vmax=vmax))
    ax.set_xticks(range(len(windows)),windows);ax.set_xlabel(r'$k$ window ($\AA^{-1}$)')
    ax.set_title(f'{title}\nReview trigger > {threshold:g}')
    for y in range(len(results)):
        for x in range(len(windows)):
            value=data[y,x]
            ax.text(x,y,f'{value:.3f}',ha='center',va='center',fontsize=8,
                    color='white' if value>.62*vmax else 'black',
                    fontweight='bold' if value>threshold else 'normal')
    fig.colorbar(im,ax=ax,shrink=.8)
axes[0].set_yticks(range(len(results)),labels)
fig.suptitle('Oxide references: radius sensitivity versus SCF 5 / cluster 8 A\nMaximum over individual first-shell paths; native degeneracy, sigma2=0, k^2 weighting. Candidates remain unvalidated.',fontsize=11)
fig.savefig(ROOT/'sensitivity-summary.png',dpi=180)
fig.savefig(ROOT/'sensitivity-summary.svg')
