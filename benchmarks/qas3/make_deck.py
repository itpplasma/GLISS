"""Generate a TERPSICHORE deck from the QAS3 benchmark deck with a custom stability table.

usage: make_deck.py TEMPLATE OUT MMS NSMIN NSMAX IVAC AL0

TEMPLATE is QAS3_fx_cur_bench_n1.data from the TERPSICHORE repository.
Stability table: all (m,n) with 0<=m<=MMS, NSMIN<=n<=NSMAX, n % nfp != 0 (N=1 family for nfp=3),
m=0 only for n>0.
"""
import sys, os
base = open(sys.argv[1]).read().split('\n')
out, mms, nsmin, nsmax = sys.argv[2], int(sys.argv[3]), int(sys.argv[4]), int(sys.argv[5])
ivac, al0 = int(sys.argv[6]), float(sys.argv[7])
nfp = 3
lines = list(base)
# header line 4 (index 3): MM NMIN NMAX MMS NSMIN NSMAX NPROCS INSOL
f = lines[3].split()
f[3], f[4], f[5] = str(mms), str(nsmin), str(nsmax)
lines[3] = '     ' + ''.join(f"{int(x):6d}" for x in f)
f = lines[6].split()
f[2] = str(ivac)
lines[6] = '     ' + ''.join(f"{int(x):6d}" for x in f[:5]) + ''.join(f"{int(x):7d}" for x in f[5:])
# locate stability table: after the line containing 'TABLE OF FOURIER COEFFIENTS FOR STABILITY'
k = next(i for i, l in enumerate(lines) if 'STABILITY DISPLACEMENTS' in l)
hdr = k + 2  # 'C M=' line at k+1
start = k + 3
end = next(i for i in range(start, len(lines)) if lines[i].startswith('C'))
modes = []
rows = []
for n in range(nsmin, nsmax + 1):
    flags = []
    for m in range(56):
        on = (m <= mms) and (n % nfp != 0) and (m > 0 or n > 0)
        flags.append(1 if on else 0)
        if on: modes.append((m, n))
    rows.append('     ' + ''.join(f"{v:2d}" for v in flags) + f"{n:3d}")
lines[start:end] = rows
# eigen line: NEV NITMAX AL0 EPSPAM IGREEN MPINIT (last numeric line)
j = next(i for i, l in enumerate(lines) if 'NEV NITMAX' in l) + 1
f = lines[j].split()
lines[j] = f"{int(f[0]):7d}{int(f[1]):7d}{al0:12.3E}{float(f[3]):12.4E}{int(f[4]):7d}{int(f[5]):7d}"
open(out, 'w').write('\n'.join(lines))
print(len(modes), 'modes')
open(out + '.modes', 'w').write('\n'.join(f"{m} {n}" for m, n in modes) + '\n')
# optional equilibrium Boozer-table truncation: EQ_MN="M N"
if os.environ.get('EQ_MN'):
    em, en = map(int, os.environ['EQ_MN'].split())
    lines = open(out).read().split('\n')
    f = lines[3].split(); f[0], f[1], f[2] = str(em), str(-en), str(en)
    lines[3] = '     ' + ''.join(f"{int(x):6d}" for x in f)
    k = next(i for i, l in enumerate(lines) if 'BOOZER COORDINATES' in l)
    start = k + 4
    end = next(i for i in range(start, len(lines)) if lines[i].startswith('C'))
    rows = ['     ' + ''.join(f"{(0 if (m == 0 and n < 0) else 1):2d}" for m in range(37)) + f"{n:3d}" for n in range(-en, en + 1)]
    lines[start:end] = rows
    open(out, 'w').write('\n'.join(lines))
