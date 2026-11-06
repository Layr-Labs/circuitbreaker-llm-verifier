#!/usr/bin/env python3
from fractions import Fraction as F
from math import ceil,floor,isqrt,comb
from pathlib import Path
import json,sys,re,hashlib
EPS=F(1,2**20); NE=F(1,10**6)
SL=dict(zip(('a','s','s0','hc','d','k'),map(F,
 '16829079/524288,32879555/1048576,33062739/1048576,303862521/524288,3210457771/1048576,26132174653767/131072'.split(',')))) # affine-and-standard-wrappers.md:471-476
DX,DY,RATE=120,28,529567 # proof.md:715
QC=1181 # proof.md:687
QS=F(293,2) # proof.md:442
HT,BL=19,14 # linear-and-bilinear-primitives.md:476
RHO=468 # sampler-and-disperser.md:296
VE={
 'I':'301/250,12673/10000,1277/1000,32/25,361/250,1597/1000,8/5,797/500,30037/10000', # linear-and-bilinear-primitives.md:976
 'F':'152/125,32/25,32/25,32/25,1517/1000,38/25,1597/1000,38/25,3321/1000', # linear-and-bilinear-primitives.md:977
 'G':'29819/10000,43249/10000', # linear-and-bilinear-primitives.md:978
 'J':'32951/10000,10251/2000', # linear-and-bilinear-primitives.md:979
 'Ext0':'702/25,702/25,702/25,702/25,702/25,702/25,702/25,702/25,702/25,702/25,702/25,702/25,702/25,702/25,69523/2500,275311/10000,6669/250,256279/10000,14281/625,49089/2500,168499/10000', # linear-and-bilinear-primitives.md:980
 'Ext2':'29859/2500,6919/625,54313/5000,106941/10000,24307/2500,46013/5000', # linear-and-bilinear-primitives.md:981
 'L0':'1237/1000,11/5,11/5,11/5,11/5,11/5,11/5,11/5,11/5,11/5,11/5,11/5,11/5,143/50', # linear-and-bilinear-primitives.md:982
 'Lr':'57027/10000,13243/2500,9603/2000,10877/2500,10877/2500,10877/2500,10877/2500,10877/2500,10877/2500,10877/2500,40389/10000,33109/10000,33109/10000,33109/10000,57027/10000', # linear-and-bilinear-primitives.md:983
 'E':'469/200,39/25,39/25,39/25,39/25,179/50', # linear-and-bilinear-primitives.md:984
}
def lg(x):
 x=F(x);k=max(0,x.numerator.bit_length()-x.denominator.bit_length())
 return k+(2**k<x)
def env(v,ag=False):
 v=list(map(F,v)); h=len(v); A=F(2); B=F(1,2**h)
 for R in v[::-1]:
  la=2*R/(R-1); A=la*(A+1+(R+1)/2); B*=la
 R=v[-1]; kap=(R+1)/(R-1)
 if ag: c=F(2,2**h); pi=6+sum(r+1 for r in v)+F(1,8192)+(h+2)*EPS
 else:
  c=kap/2**h; ch=-F(31,2)*kap+(R+1)/(2*R)+F(3,2)
  pi=2+sum(r+1 for r in v)+kap*(3+(R+1)/2)+(R+1)/2+max(0,h+ch)*EPS
 return dict(c=c,pi=pi,lam=B,g=A)
def inner(v,ag=True):
 d=SL.copy()
 d['env']={n:env(v[n],ag) for n in('I','F','G','J')}
 return d
def nb(v,inn,z=128):
 x,y=env(v['Ext0']),env(v['Ext2']); eta=F(268435456,1975781445633); theta=eta/7
 df=inn['d']+NE; kf=inn['k']+NE
 M0=y['c']*df+y['pi']*(1+NE)+NE; K2=y['lam']*df+y['g']*(1+NE)+NE
 if theta<=7*x['c']: raise ValueError('NB denominator')
 R=(M0+x['pi']*(1+NE)+7+7*NE)/(theta-7*x['c'])+NE
 d0=7*x['c']*R+x['pi']*(1+NE)+NE
 d1=2*(M0+d0+3)+1
 K00=7*x['lam']*R+x['g']*(1+NE); Kother=kf+14*R+z+2+NE
 dc=2*d0+2*d1+df+K2+1
 return dict(dF=df,kF=kf,M0=M0,K2=K2,Rnew=R,d0=d0,d1=d1,K00=K00,Kother=Kother,dc=dc,dC=ceil(dc),kC=ceil(max(K00,Kother)),px0=x,px2=y)
def outer(v,nbb,q,ag=True,causal=True):
 C=4*(q+1); t=4*C-1; j=lg(t); sc=nbb['dC']; kc=nbb['kC']
 E=env(v['E'],ag); R=env(v['Lr'],ag); Z=env(v['L0'],ag)
 dx,dy,r=DX,DY,RATE
 if dx>128: raise ValueError('output cap')
 s0=ceil(Z['c']*sc+Z['pi'])
 TR=(t+1)*(s0+sc+(0 if causal else dx)+2*j*(dy+dx)); TL=(t+1)*(sc+2*j*dy)
 rows=dict(dy_E=E['c']*dx+E['pi'],dx_Lr=R['c']*r+R['pi'],dx_Lm=2*dy+4+NE,s0_L0=Z['c']*sc+Z['pi'],r_3a=(t+1)*dy+dy+3,r_dy=dy,r_Km=dy+2,
  d_3c=(t+1)*(s0+sc)+kc+1,d_3d=TR+E['lam']*dx+E['g']+t*dx+1,d_s0=s0,d_KE=E['lam']*dx+E['g'],d_kC=kc,
  k_3b=TL+R['lam']*r+R['g']+t*r+1,k_K0=Z['lam']*sc+Z['g'])
 d=ceil(max(val for ky,val in rows.items() if ky.startswith('d_')))
 k=ceil(max(d,*[val for ky,val in rows.items() if ky.startswith('k_')]))
 return dict(q=q,C=C,t=t,j=j,dx=dx,dy=dy,r=r,sC=sc,kC=kc,s0=s0,TR=TR,TL=TL,dO=d,kO=k,rows=rows,env=dict(E=E,Lr=R,L0=Z))
def base(ou,Q,h,b,Ro=None):
 Q=F(Q); beta=F(2*(ou['q']+8),7); M=beta+2*ou['q']; p=3*M
 co=(6*p+1)/Q+2+EPS; delta=ceil(ou['dO']*co); kappa=ceil(ou['kO']*co)
 alpha=Q*delta/(4**h*(1-EPS)); R0=2**40
 rc=b+2+(b+lg(b)+4)*EPS; sc=max(rc,alpha/(b-1)+2*EPS)
 i0=F(R0+1,R0)*(Q+1+EPS)+EPS
 if Ro is None: Ro=(beta-i0-sc-rc-h*EPS)/(h*(1+EPS))-1
 Ro=F(Ro)
 if Ro<=3: raise ValueError('base rate')
 spo=(Ro+1)*(1+EPS)+EPS; seed=i0+h*spo+sc+rc
 if seed>beta: raise ValueError('base seed')
 lo=4*Ro/(Ro-3); g=lo*(1+4*EPS)+3*spo/4; kt=b*(sc-1)
 kfac=lo**h*(kt+g/(lo-1)); baseH=ceil(kfac/Q)+7
 T=ceil((12*p+2)/Q+EPS); rawH=4*ou['C']*delta+kappa+T
 return dict(Q=Q,beta=beta,M=M,p=p,h=h,b=b,Ro=Ro,alpha=alpha,R0=R0,rc=rc,sc=sc,i0=i0,spo=spo,seed=seed,lo=lo,g=g,kt=kt,kfac=kfac,ceilk=ceil(kfac),baseH=baseH,rawH=rawH,H=max(baseH,rawH),K=4*max(baseH,rawH),co=co,delta=delta,kappa=kappa,T=T)
P='proof.md'; W='affine-and-standard-wrappers.md'; V='linear-and-bilinear-primitives.md'
S='sampler-and-disperser.md'; R='restricted-nonmalleability-and-breaking.md'
TEXT={};ROOT=None
checks=0;ac=0;fail=[];groups={}
def ck(nm,a,b=True,op='==',text=False):
 global checks,ac
 ok={'==':a==b,'<=':a<=b,'>=':a>=b,'<':a<b,'>':a>b}[op]
 checks+=1;ac+=not text;g=nm.split('.')[0];groups[g]=groups.get(g,0)+1
 if not ok:fail.append(dict(check=nm,left=str(a),right=str(b),op=op))
def batch(prefix,rows):
 for nm,*args in rows:ck(prefix+nm,*args)
def tc(nm,a,b=True,op='=='):ck(nm,a,b,op,True)
def num(f,rx):
 m=re.search(rx,TEXT[f],re.M);tc('parse.'+rx,m is not None)
 return F(m[1]) if m else None
def pt(f,token):
 if TEXT:tc('text.'+token,token in TEXT[f])
def et(e):return ','.join(k+'='+str(e[k]) for k in('c','pi','lam','g'))
def run():
 v={n:list(map(F,s.split(','))) for n,s in VE.items()}
 i=inner(v); n=nb(v,i); q=QC;o=outer(v,n,q)
 Q=QS;h=HT;b=BL
 z=base(o,Q,h,b); beta=z['beta'];p=z['p'];C=o['C'];t=o['t'];j=o['j'];H=z['H'];K=z['K'];A=z['alpha'];Ro=z['Ro'];U=2**40
 for nm,vec in v.items():
  ag=nm not in('Ext0','Ext2');en=env(vec,ag);prod=F(1);ss=F(0)
  for ix,rr in enumerate(vec):
   ck('PXT.rate',rr,1,'>');prod*=2*rr/(rr-1);ss+=prod*(rr+3)/2
  ck('PXT.lambda',en['lam'],prod/2**len(vec));ck('PXT.g',en['g'],2*prod+ss)
  if ag:
   ck('PXG.terminal-min',vec[-1],F(5,3),'>=');ck('PXG.terminal-max',vec[-1],64,'<=')
   ck('PXG.error-sum',3*2**len(vec)-1,2**(len(vec)+4),'<')
   output={'I':i['hc']+EPS,'F':i['hc']+EPS,'G':i['a']+EPS,'J':i['s']+EPS,'L0':o['sC'],'Lr':o['r'],'E':o['dx']}[nm]
   ck('PXG.output',output,256*2**len(vec),'<=')
  if nm in i['env']:pt(W,f' {nm}: '+et(en))
  elif nm in o['env']:pt(P,f' {nm}: '+et(en))
  else:pt(P,f'{nm}=PXT(vNB'+('0' if nm=='Ext0' else '2')+'): '+et(en))
 I,FF,G,J=[i['env'][a] for a in('I','F','G','J')]
 rows=[('a',FF['c']*(i['hc']+EPS)+FF['pi']),('a',J['c']*(i['s']+EPS)+J['pi']),('s',G['c']*(i['a']+EPS)+G['pi']),('s0',I['c']*(i['hc']+EPS)+I['pi']),('hc',(15+J['lam'])*(i['s']+EPS)+J['g']+1),('hc',F(128)),('hc',i['s']+EPS),('d',15*(i['s0']+EPS)+(77+G['lam'])*(i['a']+EPS)+G['g']+1),('d',i['s0']+EPS),('k',(89+FF['lam'])*(i['hc']+EPS)+FF['g']+1),('k',I['lam']*(i['hc']+EPS)+I['g']),('k',i['d']+EPS)]
 for ky,val in rows:
  ck('F.slot-'+ky,i[ky],val,'>=');pt(W,f' {ky}: [{val}]w <= ({i[ky]})w.')
 for ky in('a','s','s0','hc','d','k'):pt(W,f' {ky}=ceil(({i[ky]})w).')
 for ky,val in o['rows'].items():
  slot=o[dict(dy='dy',dx='dx',s0='s0',r='r',d='dO',k='kO')[ky.split('_')[0]]]
  ck('A.slot-'+ky,slot,val,'>=');pt(P,f' {ky}: [{val}]Om <= {slot}Om')
 for ky in('dy','dx','r','sC','s0','kC'):pt(P,f'{ky}={o[ky]}Om')
 for ky,k in [('d_O','dO'),('k_O','kO')]:pt(P,f'{ky}={o[k]}Om')
 pt(P,f'TR={o["TR"]},TL={o["TL"]}')
 pt(W,'TR=TR0+H(t+1)(2dy+2dx)-(t+1)dx');pt(W,'T=TR+TL+(t+1)dx')
 for nr in(0,1,2,j):
  rt=o['s0']+o['sC']+o['dx'];lf=o['sC'];mx=0
  for rr in range(nr):
   rt+=o['dy'];lf+=o['dy'];mx=max(mx,rt)
   rt+=o['dx']+o['dy'];lf+=o['dy'];mx=max(mx,rt);rt+=o['dx']
  if nr:ck('TR.causal',mx,o['s0']+o['sC']+2*nr*(o['dy']+o['dx']))
  ck('TR.full',rt+lf,o['s0']+2*o['sC']+o['dx']+nr*(4*o['dy']+2*o['dx']))
 for ky,label in [('dF','dF'),('kF','kF'),('M0','M0'),('K2','K2bar'),('Rnew','Rnew'),('d0','d0bar')]:
  if ky in('dF','kF'):pt(P,label+'='+str(n[ky]))
  elif TEXT:tc('NB.print-'+ky,str(n[ky]) in TEXT[P])
 if TEXT:tc('NB.dc-text',str(n['dc']) in TEXT[P]);tc('NB.d1-text',str(n['d1']) in TEXT[P])
 ck('NB.branch0',n['K00'],n['kC'],'<');ck('NB.branch1',n['Kother'],n['kC'],'<');ck('NB.seed',n['dc'],n['dC'],'<')
 pt(W,'c0='+str(n['px0']['c'])+',pi0='+str(n['px0']['pi']))
 pt(W,'7c0='+str(7*n['px0']['c'])+'<theta=268435456/13830470119431')
 eta=F(268435456,1975781445633);theta=eta/7;gamma=F(256,20481)
 ck('NB.denominator',theta,7*n['px0']['c'],'>');ck('NB.pi0',n['px0']['pi'],16,'>=')
 batch('NB.',[
 ('unit',max(len(v['Ext0']),len(v['Ext2'])),38,'<='),('z',o['dx'],128,'<='),
 ('advice-numerator',F(38,7)+1+F(53,U),7,'<='),('log',8*41+240,NE*U,'<='),
 ('ceilings',12,NE*U,'<='),('logconstant',lg(ceil(n['Rnew']*n['M0'])+1)+1,240,'<'),
 ('quarter',2*41+240,F(U,4),'<'),('prime',7*2**20/eta,3200000,'>'),
 ('capacity',F(8,7)*eta*(1+eta),1,'<'),('floor',53-F(17*4+8,7)-1,F(288,7),'=='),
 ('advice',F(288,7),40,'>'),('error',14,16,'<'),('final-error',3+2*(1+1+2+2),16,'<'),
 ('deficit',112,64,'>'),('absorption',112*2**20,64*2**20,'>')])
 for f in(P,V,S):pt(f,'Q_S='+str(Q))
 for f in(P,V,S):
  pt(f,'alpha='+str(A));pt(f,'Ro='+str(Ro));pt(f,'ceilkfac='+str(z['ceilk']))
 for ky in('h','b'):pt(V,ky+'='+str(z[ky]))
 for ky in('B0','M','p','q','C','t'):pt(P,ky+'='+str({'B0':beta,'M':z['M'],'p':p,'q':q,'C':C,'t':t}[ky]))
 pt(P,'c_O=(6p+1)/Q_S+2+2^-20='+str(z['co']))
 pt(P,f'delta_O=ceil({o["dO"]}c_O)={z["delta"]}')
 pt(P,f'kappa_O=ceil({o["kO"]}c_O)={z["kappa"]}')
 pt(P,'T=ceil((12p+2)/Q_S+2^-20)='+str(z['T']))
 pt(P,f'Rawbranch={z["rawH"]},basebranch={z["baseH"]},H={H}.');pt(P,'K=4H='+str(K))
 ck('result.H',H,7412884525208893);ck('result.K',K,29651538100835572) # bound.json:1; proof.md:930
 ck('result.stronger',H,10342869631519417,'<')
 batch('base.',[
 ('Ro',Ro,3,'>'),('q',F(7,2)*beta-8,q,'=='),('C',C,4*(q+1),'=='),('M',z['M'],beta+2*q,'=='),
 ('beta',beta,511,'<'),('Q',Q,2,'>='),('seed',z['seed'],beta,'=='),('output',4**h*A*(1-EPS),Q*z['delta'],'=='),
 ('ceil',z['kfac'],z['ceilk'],'<='),('ceil-min',z['ceilk']-1,z['kfac'],'<'),('bits',z['ceilk'],2**84,'<'),
 ('padding',(z['R0']+1)*z['ceilk']+z['i0'],2**128,'<'),('error',3*4**h,2**64,'<'),
 ('hash',4*(b-1),2**(lg(b)+2),'<'),('error-sum',1+F(4*(4**h-1),3)+4**h+F(4**h,4),3*4**h,'<'),
 ('p',p,3*Q,'>'),('source',H,2*z['delta'],'>='),('walk-support',z['M']*U,134217728,'>'),
 ('native',13*2**j-8,26*t,'<'),('Lm-log',o['r']*o['dy']+1,2**33,'<='),
 ('Lm-overhead',66+4*41,NE*U,'<='),('sampler-failure',F(1,128),F(1,100),'<')])
 e0=33+lg(C)+lg(26*t);allow=max(e0+42+36,74+2*lg(C)+72)
 rb=RHO;ck('dyadic.rho-constant',rb,468)
 for nm,ic,sl,bd in [('a',144,8,EPS),('rho',rb,10,EPS),('rho0',Q+234,9,EPS),
  ('outer',20*(288*q+8)+allow,2*(288*q+8),EPS*Q),('runtime',480*q+6,48*q,p),
  ('NB-log',248,8,NE),('NB-quarter',242,2,F(1,4)),('Lm',70,4,NE)]:
  an=ic+sl*40;ck('dyadic.'+nm,an,bd*U,'<');ck('dyadic.slope-'+nm,sl,an,'<=')
 ck('PXG.log-anchor',22+2*20,2**6,'<=');ck('PXG.log-slope',2,22+2*20,'<=')
 ck('PXG.dimension',1199*2**20+9,2048*2**20,'<=');ck('PXG.product',2048*256,2**19)
 ends={s for k in range(40,81) for s in(2**k,2**(k+1)-1)}
 ck('dyadic.endpoints',len(ends),82)
 for u in sorted(ends|{2**k+d for k in(40,41,63,80) for d in range(-1,8) if 2**k+d>=U}):
  k=u.bit_length()-1;u1=2**(128+lg(u))
  sv=2*ceil(beta*u/2);dv=sv-beta*u;lc=lg(sv+1);rf=u+8*lc;aa=rf+64
  rho=aa+2*lg(u1)+2;ell=ceil(A*u);rt=b*rho+2*aa+lg(b)+2;st=max(rt,ceil(F(ell,b-1)))
  so=ceil((Ro+1)*rho);ki=b*(st-rho);walk=sv+2*q*rf;tag=3*walk
  batch('BT.',[('log-dimension',sv+1,2**(k+10),'<='),('lgS',lg(u),k+1,'<='),
   ('rho-definition',rho,rf+322+2*lg(u),'=='),('rho-dyadic',rho-u,10*k+rb,'<='),
   ('rho-conservative',aa+322+2*lg(u)-u,10*k+rb,'<='),('a-dyadic',aa-u,8*k+144,'<='),
   ('rho-envelope',rho,(1+EPS)*u,'<='),('a-envelope',aa,(1+EPS)*u,'<='),
   ('rho-lower',rho,u,'>='),('So',so,z['spo']*u,'<='),('So-denominator',so,4*rho,'>'),
   ('ratio',F(4*(so-rho),so-4*rho),z['lo'],'<='),('a-rounding',aa+3,(1+4*EPS)*u,'<=')])
  batch('BT.',[('dv0',dv,0,'>='),('dv2',dv,2,'<'),('walk-even',walk%2,0,'=='),('ellc',lc,u.bit_length()-1+10,'<='),
   ('walk-margin',F(7,16)*dv+9-8*lc,0,'<'),('r',rt,z['rc']*u,'<='),('s',st,z['sc']*u,'<='),('Kh',ki,z['kt']*u,'<='),
   ('subset',rt,st,'<='),('leaf',ell,(b-1)*st,'<='),('terminal',b,ki,'<=')])
  for level in range(h,0,-1):
   bi=ceil(F(4*(ki+aa+3),so-4*rho));mi=ceil(F(bi*so,4));kp=3*mi+ki+aa
   batch('BT.',[('ceil',bi*(F(so,4)-rho),ki+aa+3,'>='),('capacity',kp,bi*(so-rho),'<'),('count',bi,kp,'<='),
    ('child',mi,ki,'>='),('pad',mi,kp,'<='),('envelope',kp,z['lo']*ki+z['g']*u,'<')])
   ki=kp
  for L in(ceil(Q*u),ceil(Q*(u+1))-1):
   for shift in(2,3):
    rho0=aa+L+shift+lg(ki)+2;sinit=ceil(F(z['R0']+1,z['R0'])*rho0)
    for nm,lf,rhs in [('lgK0',lg(ki),84+lg(u)),('rho0-dyadic',rho0-(Q+1)*u,Q+234+9*k),
     ('rho0-envelope',rho0,(Q+1+EPS)*u),('s0',sinit,z['i0']*u),('source',lg(H*L),L),
     ('source-K0',ki,H*L),('seed-pad',sinit+h*so+st+rt,sv)]:ck('BT.'+nm,lf,rhs,'<=')
    u0=sinit-rho0
    batch('BT.',[('floor',floor(L/Q),u,'=='),('K0',ki,z['kfac']*u,'<='),('R0',z['R0'],rho0,'<='),
     ('u0',u0,F(rho0,z['R0']),'>='),('ratio',F(sinit,u0),z['R0']+1,'<='),('padding',(z['R0']+1)*ki+sinit,2**(128+lg(u)),'<='),
     ('seed',sinit+h*so+st+rt,beta*u,'<='),('output',4**h*ell,z['delta']*L,'>='),('entropy',ki+7,H*L,'<='),
     ('tag',tag,2*L,'>='),('runtime',tag,p*L,'<='),('outer-overhead',2*(288*q+8)*lc+allow,EPS*L,'<=')])
 batch('',[('SP.gamma',F(1,4)-20*gamma,F(1,81924),'=='),('SP.guard',37*81924,3200000,'<'),
 ('SP.subset',(F(1,16)-gamma)*3200000,3,'>'),('SP.dimension',18*gamma+F(33,3200000),F(1,4),'<'),
 ('SP.power',F(7,3)-F(1,16)+F(1,4),F(121,48),'=='),('SP.eta',(92+EPS)*eta,gamma,'=='),('SP.cap',eta,F(1,7000),'<='),
 ('SP.cardinal',5+3+1,9,'=='),('SP.sixteen',9,16,'<'),('SP.endpoint',F(15+4)*3/8,F(57,8),'=='),('SP.endpoint-loss',F(57,8),8,'<'),
 ('CZ.units',1-F(1021,7000),F(5979,7000),'=='),('CZ.plane',F(31,16)-F(1020,7000)-1,F(2217,2800),'=='),
 ('CZ.projection',F(31,16)-F(1041,7000)-F(5,4),F(7543,14000),'=='),('CZ.energy',2**63+2**57,2**64,'<'),
 ('CZ.smoothing',6+37,45,'<='),('CZ.hit',35,36,'<'),('CZ.root16',16,2**(20*64),'<'),
 ('prime.anchor',F(32,5),F(99,16),'>'),('prime.square',631,26**2,'<'),
 ('native.union',64*2**16*2*1024,2**33,'=='),('native.collision',F(2,2048)+F(1,4),F(257,1024),'=='),
 ('native.r',F(257,1020)+F(3,1024*2**256),F(13,50),'<'),('native.mono',6-8*F(13,50),0,'>'),
 ('native.quartic',6*F(13,50)-4*F(13,50)**2-2,F(-444,625),'=='),('native.margin',F(-444,625)+F(1,1000)+F(12,1024),F(-1,2),'<='),
 ('native.g',F(4,2**4095),F(1,1000),'<'),('native.delete',F(1,2**27),F(1,10000**2),'<'),('native.perturb',F(432,10000),F(1,8),'<'),
 ('native.moment',3+F(3,1024),4,'<='),('native.success',1-F(2,100),F(98,100),'=='),
 ('F.entropy',5*15+14,89,'=='),('F.seed',5*15+2,77,'=='),('F.error',76,128,'<'),('F.round1',1+12,13,'=='),
 ('F.round2',2*13+8,34,'=='),('F.round3',2*34+8,76,'==')])
 ck('SP.constants',(8**4,2**15,2**17*2**20,2**15*2**18),(4096,32768,2**37,2**33))
 ck('native.ratio',F(257,1024)*F(256,255),F(257,1020))
 ck('native.bias',2*(F(1,1024*2**256)+F(1,2**33*C*2**768)),F(3,1024*2**256),'<=')
 ck('native.deletion',F(1,2**28)*F(256,255),F(1,2**27),'<')
 primes=[2,3,5,7,13,23,43,83,163,317,631]
 for ix,pp in enumerate(primes):
  ck('prime.trial',all(pp%d for d in range(2,isqrt(pp)+1)))
  if ix:ck('prime.chain',pp,2*primes[ix-1],'<')
 codes=['11000','10100','10010','01100','01010','00110','10001','01001']
 pt(W,','.join(codes));ck('F.distinct',len(set(codes)),8)
 for aa in codes:
  ck('F.weight',aa.count('1'),2)
  for bb in codes:
   if aa!=bb:ck('F.antichain',any(x=='1' and y=='0' for x,y in zip(aa,bb)))
 for d in range(1,20):
  for nn in range(1,40):
   vv=2*nn//(d+2);ck('SP.dimension-count',(vv+1)*comb(d+2,2),nn*(d+1),'>')
 for d0 in range(1,12):
  for d1 in range(d0):
   TT=d0+d1+7
   ck('SP.codimension',comb(TT+2,2)-comb(TT-d0+2,2)-comb(TT-d1+2,2)+comb(TT-d0-d1+2,2),d0*d1)
 def px(nm,nn,zz,ee):
  vec=v[nm];hh=len(vec);aa=ee+hh+4;rho=aa+2*lg(4*nn)+2
  slots=[ceil((rr+1)*rho) for rr in vec];kk=ceil(F(zz,2**hh))+2*aa
  for ix in range(hh-1,-1,-1):
   nxt=kk;kk=ceil(2*vec[ix]/(vec[ix]-1)*(kk+aa+F(slots[ix],2)+1))
   mm=ceil(F(ceil(F(kk,slots[ix]-rho))*slots[ix],2))
   ck('PX.split',kk,mm+nxt+aa,'>=');ck('PX.child',mm,nxt,'>=')
   if ix==hh-1:leaf=mm
  ww=max(2**20,ee+hh+6+2*lg(4*nn));ag=nm not in('Ext0','Ext2')
  if ag: tail=2*(lg(leaf*ceil(F(zz,2**hh))+1)+ceil(F(zz,2**hh))+2*aa)
  else:tail=leaf
  dd=2*rho+sum(slots)+tail;en=env(vec,ag)
  ck('PX.threshold',kk,en['lam']*zz+en['g']*ww,'<=')
  if not ag or zz<=256*2**hh*ww:ck('PX.seed',dd,en['c']*zz+en['pi']*ww,'<=')
  ck('PX.padding',2*kk+2*rho,4*nn,'<=');return dd,kk
 L=ceil(Q*U);sv=2*ceil(beta*U/2);rf=U+8*lg(sv+1);tag=3*(sv+2*q*rf)
 ee0=6*tag+rf+33+lg(C)+lg(26*t);om=ee0+2*L+42
 for nm in v:px(nm,2**256,1,1)
 for nn,zz,adv,ex in [(2**256,1,1,1),(o['dO']*om,o['dx']*om,tag+2+lg(C),ee0)]:
  ee=ex+4;wn=max(2**20,ee+44+2*lg(4*nn));wf=max(2**20,ee+2*lg(4*nn)+32,ceil(F(zz,128)))
  m1,kf=ceil(i['d']*wf),ceil(i['k']*wf);m0,k2=px('Ext2',nn,m1,ee)
  sig=ee+4;rr=ceil((m0+n['px0']['pi']*wn+F(38,7)*ee+adv+53)/(theta-7*n['px0']['c']))
  ss=floor((eta*rr-17*sig-8)/7);d0,k00=px('Ext0',nn,7*rr,ee)
  el=ceil(F(d0+ee-1,2))+ee;d1=2*(lg(rr*m0+1)+m0+2*el);d2=d0+d1+m1+k2+ee
  dc=d0+d1+d2;kc=max(k00,kf+14*rr+zz+2*ee)
  batch('NBexec.',[('I1',ss,m0+2*el,'>='),('advice',ss,adv+40,'>='),('density',F(adv+sig,rr+adv),eta/2,'<='),
   ('source',rr-sig,(1-eta)*(rr+adv+1)+1,'>='),('capacity',floor(F(8,7)*eta*(rr+adv))-8,8*ss+17*sig,'>='),
   ('deficit',eta*(rr+adv),64*sig,'>='),('prime',rr+adv,3200000,'>'),('I4',d2-d0-d1-m1-ee,k2,'=='),
   ('I5',kc-14*rr-zz-2*ee,kf,'>='),('I6',F(d0+ee-1,2)-el,-ee,'<='),('guard-source',kc,nn,'<='),('guard-seed',d2,nn,'<=')])
  if ex==ee0:
   ck('NBexec.B9seed',dc,o['sC']*om,'<');ck('NBexec.B9source',kc,o['kC']*om,'<')
 def mul(x,y):
  z=0
  while y:
   if y&1:z^=x
   x<<=1;y>>=1
   if x&8:x^=11
  return z
 for bb in(2,3):
  for code in range(1,8**bb):
   coeff=[(code>>(3*ix))&7 for ix in range(bb)];roots=0
   for yy in range(8):
    val=0
    for c in coeff[::-1]:val=mul(val,yy)^c
    roots+=val==0
   ck('terminal.roots',roots,bb-1,'<=')
  source=list(range(8**bb//2));source[-1]=8**bb-1;tv=F(0)
  for yy in range(8):
   counts={}
   for code in source:
    z0=code&7;yp=yy;out=0
    for ix in range(1,bb):out|=(((code>>(3*ix))&7)^mul(yp,z0))<<(3*(ix-1));yp=mul(yp,yy)
    counts[out]=counts.get(out,0)+1
   tv+=sum(abs(F(counts.get(out,0),len(source))-F(1,8**(bb-1))) for out in range(8**(bb-1)))/16
  ck('terminal.TV',tv*tv,F(bb-1,16),'<=')
 if TEXT:
  for nm,tag in [('I','MI'),('F','MF'),('G','MG'),('J','MJ'),('Ext0','NB0'),('Ext2','NB2'),('L0','L0'),('Lr','Lr'),('E','E')]:
   pt(V,' v'+tag+'=('+VE[nm]+')')
  for f,rx,val in [(P,r'p=[0-9/]+,q=(\d+),C=',q),(P,r'Q_S=([0-9/]+)',Q),
   (V,r' h=(\d+), B=4\^h',h),(V,r', b=(\d+), R0=',b),(S,r'rho-L<=10k\+(\d+)',rb)]:
   tc('literal.'+rx,num(f,rx)==val)
  for f in TEXT:
   tc('refs.removed-name','sum-product-and-ordinary-condensation.md' not in TEXT[f])
   for m in re.finditer(r'Section (SP\.\d+(?:\.\d+)?)',TEXT[f]):
    tc('refs.'+m[1],bool(re.search(r'^#+ '+re.escape(m[1])+r'(?:[. ]|$)',TEXT[S],re.M)))
   for ky,value in [('Q_S',Q),('ceilkfac',z['ceilk'])]:
    for m in re.finditer(r'(?<![A-Za-z0-9_])'+ky+r'=([0-9]+(?:/[0-9]+)?)(?![0-9/])',TEXT[f]):tc('cross.'+ky,F(m[1]),value)
  for f in(W,P):
   pt(f,'Ext0=PXT(vNB0)');pt(f,'Ext2=PXT(vNB2)')
  for f in(V,S):pt(f,'2^(-r_F)')
  spheads=re.findall(r'^#+ (SP\.[0-9.]+)',TEXT[S],re.M);tc('refs.unique-SP',len(spheads),len(set(spheads)))
  sp=TEXT[S].split('# SP. Finite rings')[1]
  labels=re.findall(r'^ +.*\((SP-[\w]+)\)\s*$',sp,re.M)
  tc('refs.SP-labels',set(labels),set('SP-1 SP-2 SP-3 SP-4 SP-5 SP-5s SP-Mi0 SP-Mi1 SP-Mi2 SP-Mi3 SP-Mi4 SP-17 SP-18 SP-19m SP-21m SP-S2 SP-M SP-I'.split()))
  tc('refs.SP-definitions',len(labels),len(set(labels)))
  for label in re.findall(r'\((SP-[\w]+)\)',sp):tc('refs.'+label,label in labels)
  tc('refs.stale-SP',any(s in sp for s in('theoremI','theoremM','4(SP-2)')),False)
  files=[x for x in ROOT.iterdir() if x.is_file()];tc('limits.files',len(files),8,'<=')
  for f in files:tc('limits.'+f.name,f.stat().st_size,65536,'<=')
  tc('limits.total',sum(f.stat().st_size for f in files),262144,'<=')
  bj=json.loads((ROOT/'bound.json').read_text());sj=json.loads((ROOT/'submission.json').read_text())
  kd=[int(x['coefficient']['numerator']) for x in bj['guarantee']['remainder']['terms']];tc('meta.K',kd,[K])
  tc('meta.experiment',[e['source'] for e in sj.get('experiments',[])],['verify.py']) # submission.json:1
  for fn,digest in JSON_HASHES.items():tc('frozen.'+fn,hashlib.sha256((ROOT/fn).read_bytes()).hexdigest(),digest)
  tc('audited.proofs',hashlib.sha256(b''.join((ROOT/f).read_bytes() for f in TEXT)).hexdigest(),PROOF_HASH)
 return dict(all_passed=not fail,checks=checks,arithmetic_checks=ac,text_checks_run=bool(TEXT),
  skipped_reason=None if TEXT else 'Package files absent or incomplete',failures=fail,groups=groups,H=str(H),K=str(K),
  scope='Finite arithmetic witnesses; optional package integrity. Universal proofs require Markdown.')
JSON_HASHES={'bound.json': '7b4a769d770daceaacc3b13026dbcbfbb6e74df14dcbc376d9ccb4d0f2663df2', 'submission.json': 'fc87333917a761e3b9f763f436b0df2519b8d07c732a2c9d5f83098da00c3ecf'}
PROOF_HASH='38c2396defb5925044d0a1ec154246980762e4bf2ebfe8647a0679dc9363f8ff'
if __name__=='__main__':
 if json.load(sys.stdin)!={}:raise ValueError('verify.py takes {}')
 home=Path(__file__).resolve().parent
 ROOT=next((d for d in(home,Path('submission'),Path('/snapshot/files/submission')) if all((d/f).is_file() for f in(P,W,V,S,R,'bound.json','submission.json'))),None)
 if ROOT:TEXT={f:(ROOT/f).read_text() for f in(P,W,V,S,R)}
 result=run();out=json.dumps(result,sort_keys=True);assert len(out.encode())<65536
 print(out);raise SystemExit(0 if result['all_passed'] else 1)
