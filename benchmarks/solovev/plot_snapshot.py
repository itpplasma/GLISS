"""Purpose: compare current production marginality sign against independent DCON.
Dated snapshot: 2026-09-09, GLISS source 3b459021, not future checkout HEAD.
Message: marginal interval disagrees; full-volume accuracy is unqualified.
Axes q0 and negative count, dimensionless; exact raw points, no smoothing.
Okabe-Ito blue and black, circles/squares and solid/dashed. PDF/PNG/grayscale.
"""
import argparse
import csv
from pathlib import Path
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from PIL import Image
P=Path(__file__).resolve().parent
parser=argparse.ArgumentParser(description=__doc__)
parser.add_argument('--output',type=Path,required=True)
args=parser.parse_args()
OUT=args.output
OUT.mkdir(parents=True,exist_ok=True)
G=list(csv.DictReader((P/'qscan.csv').open()))
D=list(csv.DictReader((P/'dcon_reference.csv').open()))
D=sorted([r for r in D if 1.025<=float(r['q0'])<=1.1],key=lambda r:float(r['q0']))
plt.rcParams.update({'font.size':11,'axes.spines.top':False,'axes.spines.right':False})
f,a=plt.subplots(figsize=(8.5,5.8));f.subplots_adjust(left=.10,right=.97,bottom=.40,top=.80)
f.suptitle("Solov'ev n=1: current GLISS marginal bracket disagrees with DCON",fontsize=14)
a.plot([float(r['q0']) for r in G],[int(r['negative_count']) for r in G],'o-',color='#0072B2',label='Fresh GLISS two-component marginality FEEC')
a.plot([float(r['q0']) for r in D],[int(r['zero_crossings']) for r in D],'s--',color='black',label='Archived independent DCON Newcomb zeros')
a.axvspan(1.039062,1.039843,color='#999999',alpha=.35,label='Frozen DCON marginal bracket')
a.axvspan(1.05,1.1,color='#0072B2',alpha=.09,label='Observed GLISS sign-change interval (coarse)')
a.set(xlabel=r'Axis safety factor $q_0$ (dimensionless)',ylabel='Negative directions / Newcomb zeros',ylim=(-.12,2.15),title='M24 equilibrium export; ns64, FEEC degree 2, m=0..8')
a.set_yticks([0,1,2]);a.grid(axis='y',alpha=.2)
a.legend(loc='upper left',bbox_to_anchor=(0,-.22),fontsize=9)
f.text(.03,.02,'GLISS 3b459021 + frozen diagnostic patch; DCON archive f5595c06. Boundary/orientation checks pass.\nFull-volume accuracy is unqualified: implicit-surface errors 7.77e-7, 5.50e-7, 1.73e-7 at ns64/128/256\nexceed the frozen 1e-7 bound. The mismatch does not yet isolate the operator from input interpolation.\nNo physical frequency comparison or marginal-threshold acceptance is claimed.',fontsize=9)
f.savefig(OUT/'solovev_current_comparison.pdf');f.savefig(OUT/'solovev_current_comparison.png',dpi=180)
Image.open(OUT/'solovev_current_comparison.png').convert('L').save(OUT/'solovev_current_comparison_grayscale.png')
