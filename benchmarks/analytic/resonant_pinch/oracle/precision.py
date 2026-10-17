import json
from pathlib import Path
import mpmath as mp

mp.mp.dps=80
root=Path(__file__).parent
data=json.loads((root/'shooting.json').read_text())
A,Z,a,pe,rho=map(mp.mpf,['.8','1.2','.7','17000','1.7e-7'])
mu0=4*mp.pi*mp.mpf('1e-7')
b2=(A*a/Z)**2
edge=mu0*pe/Z**2
m=2
scale=Z**2/(mu0*rho*a*a)
axis_default=mp.mpf('0.0001')

def shoot(ell,gamma,n,x0=axis_default):
    h0=1+gamma*(edge+b2)
    r1=m*b2+2*b2/h0-ell/(m*h0)
    r2=-2*b2/h0+4*m*b2*b2/(ell*h0)
    alpha=((m+2)*r1+m*r2)/(4*(m+1))
    beta=(m*r1+(m+2)*r2)/(4*(m+1))
    u,w=1+alpha*x0*x0,1+beta*x0*x0
    maximum_step=(1-x0)/n
    def rhs(x,u,w):
        h=h0+(1-gamma)*b2*x*x
        diff=m*(w-u)/x
        return (diff+x*((m*b2-ell/(m*h))*w+2*b2/h*u),
                -diff+x*(-2*b2/h*w+4*m*b2*b2/(ell*h)*u))
    x=x0
    while x<1:
        # The irregular mode is -2m/x. Stay inside RK4's real-axis
        # stability interval even for the earliest axis steps.
        step=min(maximum_step,x/(2*m),1-x)
        k1=rhs(x,u,w)
        k2=rhs(x+step/2,u+step*k1[0]/2,w+step*k1[1]/2)
        k3=rhs(x+step/2,u+step*k2[0]/2,w+step*k2[1]/2)
        k4=rhs(x+step,u+step*k3[0],w+step*k3[1])
        u+=step*(k1[0]+2*k2[0]+2*k3[0]+k4[0])/6
        w+=step*(k1[1]+2*k2[1]+2*k3[1]+k4[1])/6
        x+=step
    return u

cases=[]
for gamma,index in [(mp.mpf(1),0),(mp.mpf(5)/3,1)]:
    estimate=mp.mpf(str(data['cases'][index]['roots'][0]['ell']))
    levels=[]
    for n in [1000,2000,4000]:
        value=mp.findroot(lambda ell:shoot(ell,gamma,n),
                          (estimate*(1-mp.mpf('1e-4')),estimate),
                          tol=mp.mpf('1e-50'),verify=True)
        levels.append(value)
        print('gamma',gamma,'steps',n,'ell',mp.nstr(value,40),flush=True)
    extrapolated=levels[-1]+(levels[-1]-levels[-2])/15
    relative=abs(extrapolated/estimate-1)
    assert relative<mp.mpf('2e-10')
    axis_study=[]
    for start in map(mp.mpf,['0.001','0.0001','0.00001']):
        value=mp.findroot(lambda ell:shoot(ell,gamma,4000,start),
                          (estimate*(1-mp.mpf('1e-4')),estimate),
                          tol=mp.mpf('1e-50'),verify=True)
        axis_study.append({'axis_start':str(start),'ell':str(value),
                           'relative_to_default':str(value/levels[-1]-1)})
        assert abs(value/levels[-1]-1)<mp.mpf('2e-10')
    cases.append({'gamma':str(gamma),'axis_study':axis_study,'levels':[
        {'steps':n,'ell':str(v),'lambda':str(v*scale)}
        for n,v in zip([1000,2000,4000],levels)],
        'fourth_order_ratio':str((levels[0]-levels[1])/(levels[1]-levels[2])),
        'richardson_ell':str(extrapolated),'richardson_lambda':str(extrapolated*scale),
        'relative_to_dop853':str(relative)})
(root/'precision.json').write_text(json.dumps({'decimal_precision':80,
   'axis_start':str(axis_default),
   'method':'independent classical RK4 with graded axis steps and fixed maximum step',
   'cases':cases},indent=2)+'\n')
