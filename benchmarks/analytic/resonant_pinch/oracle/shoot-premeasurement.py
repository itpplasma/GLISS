import json
from pathlib import Path
import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import brentq

ROOT = Path(__file__).parent
A, Z, a, pe, rho = .8, 1.2, .7, 17000., 1.7e-7
mu0 = 4*np.pi*1e-7
b2 = (A*a/Z)**2
edge = mu0*pe/Z**2
m = 2
scale = Z**2/(mu0*rho*a*a)

def coefficients(ell, gamma):
    h0 = 1 + gamma*(edge+b2)
    # f=x^(m-1)u, P=-ell*x^m*w/m in pressure-normalized units.
    r1 = m*b2 + 2*b2/h0 - ell/(m*h0)
    r2 = -2*b2/h0 + 4*m*b2*b2/(ell*h0)
    alpha = ((m+2)*r1 + m*r2)/(4*(m+1))
    beta = (m*r1 + (m+2)*r2)/(4*(m+1))
    return alpha,beta

def shoot(ell, gamma, x0=1e-4, rtol=1e-12, atol=1e-14, step=.01, dense=False):
    h0 = 1 + gamma*(edge+b2)
    def rhs(x,y):
        u,w=y
        h=h0+(1-gamma)*b2*x*x
        diff=m*(w-u)/x
        return [diff+x*((m*b2-ell/(m*h))*w+2*b2/h*u),
                -diff+x*(-2*b2/h*w+4*m*b2*b2/(ell*h)*u)]
    alpha,beta=coefficients(ell,gamma)
    result=solve_ivp(rhs,[x0,1],[1+alpha*x0*x0,1+beta*x0*x0],
                     method='DOP853',rtol=rtol,atol=atol,max_step=step,
                     dense_output=dense)
    if not result.success:
        raise RuntimeError(result.message)
    return result if dense else result.y[0,-1]

def phase(ell,gamma,x0=1e-4,rtol=1e-12,atol=1e-14,step=.01):
    h0=1+gamma*(edge+b2)
    def rhs(x,y):
        sine,cosine=np.sin(y[0]),np.cos(y[0])
        h=h0+(1-gamma)*b2*x*x
        value=4*b2*x/h*sine*cosine
        value+=(m*(1+b2*x*x)/x-ell*x/(m*h))*cosine*cosine
        value-=(m/x+4*m*b2*b2*x/(ell*h))*sine*sine
        return [value]
    alpha,beta=coefficients(ell,gamma)
    initial=np.arctan2(1+alpha*x0*x0,1+beta*x0*x0)
    solution=solve_ivp(rhs,[x0,1],[initial],method='DOP853',rtol=rtol,
                       atol=atol,max_step=step)
    if not solution.success:
        raise RuntimeError(solution.message)
    return solution.y[0,-1]

def run(gamma):
    lower=-4*b2*b2/(1+b2+gamma*edge)
    # These are independent 80-digit raw trial-energy/mass integrals
    # already committed in test_resonant_screw_pinch_energy.f90, m=2.
    trial=((-12184.58544802857668475/3.903558731270335941694e-8)
           if gamma==1 else
           (-11268.41463231817697414/3.89965516874636840581e-8))/scale
    lo=lower*(1+1e-9)
    assert phase(lo,gamma)<np.pi<phase(trial,gamma)
    phase_ground=brentq(lambda z:phase(z,gamma)-np.pi,lo,trial,
                        xtol=2e-15,rtol=1e-14)
    # The second root is diagnostic only; identify it by continuous phase
    # 2*pi rather than sign scans near the zero accumulation point.
    overtone_upper=trial/2
    while phase(overtone_upper,gamma)<2*np.pi:
        overtone_upper/=2
    phase_second=brentq(lambda z:phase(z,gamma)-2*np.pi,
                        phase_ground*(1-1e-8),overtone_upper,
                        xtol=2e-15,rtol=1e-14)
    grid=np.linspace(lower*(1-1e-9),phase_second*.95,160)
    roots=[]
    previous=(grid[0],shoot(grid[0],gamma))
    for ell in grid[1:]:
        value=shoot(ell,gamma)
        if previous[1]*value<0:
            root=brentq(lambda z:shoot(z,gamma),previous[0],ell,
                        xtol=2e-15,rtol=1e-14)
            solution=shoot(root,gamma,dense=True)
            values=solution.sol(np.linspace(1e-4,1-1e-7,2001))[0]
            nodes=int(np.count_nonzero(values[:-1]*values[1:]<0))
            roots.append({'ell':root,'lambda':root*scale,'nodes':nodes,
                          'bracket':[previous[0],float(ell)]})
            if len(roots)==2:
                break
        previous=(float(ell),value)
    if len(roots)<2 or roots[0]['nodes']!=0 or roots[1]['nodes']!=1:
        raise RuntimeError(('root ordering failure',roots))
    assert abs(phase_ground/roots[0]['ell']-1)<2e-10
    assert abs(phase_second/roots[1]['ell']-1)<2e-10
    initial=roots[0]['ell']
    studies=[]
    for x0 in [1e-3,1e-4,1e-5,1e-6]:
        for rtol,atol,step in [(1e-10,1e-12,.02),(1e-12,1e-14,.01),(3e-14,1e-15,.005)]:
            lo,hi=roots[0]['bracket']
            root=brentq(lambda z:shoot(z,gamma,x0,rtol,atol,step),lo,hi,
                        xtol=2e-15,rtol=1e-14)
            studies.append({'x0':x0,'rtol':rtol,'atol':atol,'max_step':step,
                            'ell':root,'relative_to_initial':root/initial-1})
    return {'gamma':gamma,'lower_bound_ell':lower,'lambda_scale':scale,
            'negative_trial_upper_bound_ell':trial,
            'continuous_prufer_ground_ell':phase_ground,
            'continuous_prufer_overtone_ell':phase_second,
            'ground_endpoint_phases':[phase(lo,gamma),phase(trial,gamma)],
            'roots':roots,'studies':studies}

if __name__=='__main__':
    result={'parameters':{'a':a,'A':A,'Z':Z,'p_edge':pe,'rho':rho,'m':m},
            'cases':[run(1.),run(5/3)]}
    (ROOT/'shooting.json').write_text(json.dumps(result,indent=2)+'\n')
    for case in result['cases']:
        print('gamma',case['gamma'],'roots',case['roots'],
              'max relative drift',max(abs(s['relative_to_initial']) for s in case['studies']))
