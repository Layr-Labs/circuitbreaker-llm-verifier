#!/usr/bin/env python3
"""Exact finite witness; universal proofs are in Markdown."""
from fractions import Fraction as F
from math import ceil,floor,isqrt,comb
from pathlib import Path
import json,sys,re,hashlib
# Budget literals: proof.md:821; linear-and-bilinear-primitives.md:509
EPS=F(1,2**20); EPS2=F(1,2**23); NE=F(1,10**6); EB=F(157,2**38);ER=F(1,2**40);N0=2**41;AB=3
SL={'a': '450483/16384', 's': '13273387/524288', 's0': '7145233/262144', 'hc': '274466439/524288', 'd': '2763214913/1048576', 'k': '147686758667065/1048576'} # affine-and-standard-wrappers.md:470
VE={'I': '33/25,2341/2000,2359/2000,1537/1250,6867/5000,12783/10000,379/250,5249/2500', 'F': '5853/5000,6093/5000,6093/5000,3039/2500,301/200,1597/1000,38/25,19599/10000', 'G': '3257/1250,7409/2500', 'J': '783/200,12801/2500,10251/2000', 'Ext0': '275543/10000,275819/10000,275267/10000,275267/10000,275267/10000,275267/10000,275267/10000,275267/10000,275267/10000,275267/10000,275267/10000,275267/10000,275267/10000,17187/625,136443/5000,272041/10000,267027/10000,129421/5000,122127/5000,219589/10000,99463/5000', 'Ext2': '5061/1250,5061/1250,5061/1250,20107/5000,1609/400,7921/2000,4823/1250,80199/10000', 'L0': '2807/2000,11/5,11/5,11/5,11/5,11/5,11/5,11/5,11/5,11/5,11/5,11/5,11/5', 'Lr': '44143/10000,10003/2500,10877/2500,10877/2500,10877/2500,10877/2500,10877/2500,10877/2500,10877/2500,40389/10000,33109/10000,33109/10000,33109/10000,23761/5000', 'E': '469/200,39/25,39/25,39/25,39/25,179/50'} # linear-and-bilinear-primitives.md:1042
QS=F('1498513008899/500000000000');QC=20;DX=77;DY=25;RATE=8428 # proof.md:442,742
P='proof.md';W='affine-and-standard-wrappers.md';V='linear-and-bilinear-primitives.md'
S='sampler-and-disperser.md';R='restricted-nonmalleability-and-breaking.md'
TEXT={};ROOT=None;fail=[];meta=[];ct=0;ac=0;gs={}
def ck(nm,a,b=True,op='==',text=False):
 global ct,ac
 ok=(a==b if op=='==' else a<b if op=='<' else a>b if op=='>' else a<=b if op=='<=' else a>=b)
 ct+=1;ac+=not text;key=nm.split('.')[0];gs[key]=gs.get(key,0)+1
 if not ok:(meta if nm.startswith('meta.') else fail).append(dict(check=nm,left=str(a),right=str(b),op=op))
def batch(prefix,rows):
 for nm,*a in rows:ck(prefix+nm,*a)
def tc(nm,*a):ck(nm,*a,text=True)
def pt(f,tok):
 if TEXT:tc('text.'+tok[:70],tok in TEXT[f])
def lg(x):
 x=F(x);k=max(0,x.numerator.bit_length()-x.denominator.bit_length());return k+(2**k<x)
def env(vec,ag=False):
 p=F(1);g=F(0)
 for r in vec:
  lam=2*r/(r-1);g+=p*(lam+(r+1)/2);p*=lam
 h=len(vec);r=vec[-1];kap=(r+1)/(r-1)
 return dict(lam=p/2**h,g=g+2*p,c=F(2,2**h) if ag else kap/2**h,
  pi=4+sum(r+1 for r in vec)+F(1,8192)+(h+2)*EPS if ag else sum(r+1 for r in vec)+3*kap+(r+1)/2+max(0,h+F(3,2)-16*kap)*EPS)
def et(e):return ','.join(k+'='+str(e[k]) for k in ('c','pi','lam','g'))
def run():
 v={k:list(map(F,s.split(','))) for k,s in VE.items()};i={k:F(s) for k,s in SL.items()}
 ev={k:env(a,k not in ('Ext0','Ext2')) for k,a in v.items()}
 I,FF,G,J=[ev[k] for k in ('I','F','G','J')];x,y=ev['Ext0'],ev['Ext2']
 eta=F(268435456,1975781445633);theta=eta/7;gamma=F(256,20481);U=2**40
 df=i['d']+NE;kf=i['k']+NE
 M0=y['c']*df+y['pi']*(1+NE)+NE;K2=y['lam']*df+y['g']*(1+NE)+NE
 rr=(M0+x['pi']*(1+NE)+7+7*NE)/(theta-7*x['c'])+NE
 d0=7*x['c']*rr+x['pi']*(1+NE)+NE;d1=2*(M0+d0+3)+1
 k00=7*x['lam']*rr+x['g']*(1+NE);kother=kf+14*rr+DX+2+NE
 dc=2*d0+2*d1+df+K2+1;sc=ceil(dc);kc=ceil(max(k00,kother))
 q=QC;C=4*(q+1);t=4*C-1;j=lg(t);Q=QS;beta=F(2*(q+8),7);p=3*(beta+2*q)
 EE,RR,ZZ=[ev[k] for k in ('E','Lr','L0')]
 s0=ceil(ZZ['c']*sc+ZZ['pi']);tr=(t+1)*(s0+sc+2*j*(DY+DX));tl=(t+1)*(sc+2*j*DY)
 rows=dict(dy_E=EE['c']*DX+EE['pi'],dx_Lr=RR['c']*RATE+RR['pi'],dx_Lm=2*DY+4+NE,s0_L0=ZZ['c']*sc+ZZ['pi'],r_3a=(t+2)*DY+3,r_dy=DY,r_Km=DY+2,
 d_3c=(t+1)*(s0+sc)+kc+1,d_3d=tr+EE['lam']*DX+EE['g']+t*DX+1,d_s0=s0,d_KE=EE['lam']*DX+EE['g'],d_kC=kc,
 k_3b=tl+RR['lam']*RATE+RR['g']+t*RATE+1,k_K0=ZZ['lam']*sc+ZZ['g'])
 do=ceil(max(val for key,val in rows.items() if key.startswith('d_')));ko=ceil(max(do,*[val for key,val in rows.items() if key.startswith('k_')]))
 co=(6*p+1)/Q+2+EPS2;delta=ceil(do*co);kappa=ceil(ko*co);T=ceil((12*p+2)/Q+EPS)
 a0=F(N0+1,N0);gb=(beta-a0*(Q+1+EB)-3*ER)/(2*(1+EB));xi=gb/(gb-1);gg=1-(xi-1)**2;aa=1/gg;bb=ceil((xi*gb+xi/U)/gg)
 Bbase=bb;raw=4*C*delta+kappa+T;bh=ceil(aa*delta+bb*(1+EB)/Q+ER);H=max(raw,bh);K=4*H
 kb=aa*delta*Q*(1+F(1,U))+bb*(1+EB)+ER;b0c=a0*(Q+1+EB)+ER
 ir=[('a',FF['c']*(i['hc']+EPS)+FF['pi']),('a',J['c']*(i['s']+EPS)+J['pi']),('s',G['c']*(i['a']+EPS)+G['pi']),('s0',I['c']*(i['hc']+EPS)+I['pi']),('hc',(15+J['lam'])*(i['s']+EPS)+J['g']+1),('hc',F(128)),('hc',i['s']+EPS),('d',15*(i['s0']+EPS)+(77+G['lam'])*(i['a']+EPS)+G['g']+1),('d',i['s0']+EPS),('k',(89+FF['lam'])*(i['hc']+EPS)+FF['g']+1),('k',I['lam']*(i['hc']+EPS)+I['g']),('k',i['d']+EPS)]
 for ky,val in ir:
  ck('inner.'+ky,i[ky],val,'>=');pt(W,f'{ky}: [{val}]w <= ({i[ky]})w.')
 for ky,val in i.items():pt(W,f'{ky}=ceil(({val})w).')
 for ky,val in rows.items():
  lim=dict(dy=DY,dx=DX,r=RATE,s0=s0,d=do,k=ko)[ky.split('_')[0]]
  ck('outer.'+ky,lim,val,'>=');pt(P,f'[{val}] <= {lim}.')
 for nm,vec in v.items():
  en=ev[nm];h=len(vec)
  for rate in vec:
   ck('PT.rate',rate,1,'>')
   sh=ceil(7*(rate+1));nn=ceil(F(3*(sh-14),2))-2
   bc=ceil(F(2*(nn+2),sh-14));mm=ceil(F(bc*sh,2));kk=mm+nn+1
   ck('PT.boundary',kk,bc*(sh-7),'<')
  f=W if nm in ('I','F','G','J') else P
  lab=nm if nm not in ('Ext0','Ext2') else nm+'=PXT(vNB'+('0' if nm=='Ext0' else '2')+')'
  pt(f,lab+': '+et(en))
  if nm not in ('Ext0','Ext2'):
   out=dict(I=i['hc']+EPS,F=i['hc']+EPS,G=i['a']+EPS,J=i['s']+EPS,L0=sc,Lr=RATE,E=DX)[nm]
   batch('PXG.',[('min',vec[-1],F(5,3),'>='),('max',vec[-1],64,'<='),('output',out,256*2**h,'<='),('error',3*2**h-2,2**(h+4),'<')])
 for nm,val in [('dF',df),('kF',kf),('M0',M0),('K2bar',K2),('Rnew',rr),('d0bar',d0)]:pt(P,str(val))
 for val in (d1,dc,k00,kother,co,gb,xi,gg,kb):pt(V if val in (gb,xi,gg,kb) else P,str(val))
 for val in (k00,kother):ck('NB.source',val,kc,'<=')
 ck('NB.seed',dc,sc,'<=');ck('NB.source-dimension',kc,do,'<=');ck('NB.seed-dimension',sc,do,'<=')
 batch('NB.',[('denominator',theta,7*x['c'],'>'),('pi0',x['pi'],16,'>='),('height',max(len(v['Ext0']),len(v['Ext2'])),38,'<='),('cap',DX,128,'<='),('advice',F(38,7)+1+F(53,U),7,'<='),('constant',lg(ceil(rr*M0)+1)+1,240,'<'),('prime',7*2**20/eta,3200000,'>'),('capacity',F(8,7)*eta*(1+eta),1,'<'),('floor',53-F(76,7)-1,F(288,7),'=='),('advice40',F(288,7),40,'>'),('error',14,16,'<'),('final-error',15,16,'<')])
 for hh in (0,1,2,j):
  rt=s0+sc+DX;lf=sc;mx=0
  for _ in range(hh):
   rt+=DY;lf+=DY;mx=max(mx,rt);rt+=DX+DY;lf+=DY;mx=max(mx,rt);rt+=DX
  if hh:ck('causal.E',mx,s0+sc+2*hh*(DY+DX))
  ck('causal.full',lf+rt,s0+2*sc+DX+hh*(4*DY+2*DX))
 batch('AF.',[('gamma',gb,2,'>'),('rank',gg*aa,1,'=='),('rounding',gg*bb-xi*gb,xi/U,'>='),('B-minimal',gg*(bb-1)-xi*gb,xi/U,'<'),('seed',a0*(Q+1+EB)+2*gb*(1+EB)+3*ER,beta,'=='),('Kbits',kb,2**48,'<'),('pad0',(N0+1)*2**48+b0c,2**90,'<'),('pad1',xi*2**48+gb*(1+EB)+EB,2**90,'<'),('error',6,2**AB,'<'),('q',F(7,2)*beta-8,q,'=='),('p',p,3*Q,'>'),('H',H,2*delta,'>='),('runtime-power',ceil(2*p/Q),97,'=='),('incidence',48*U,134217728,'>')])
 # Affine-in-k induction bases and slopes; finite endpoints supplement the proofs.
 for nm,c,b,budget in [('rho',227,10,EB),('rho0',Q+97,9,EB),('outer',10*(288*q+8)+160,2*(288*q+8),EPS*Q),('outer23',5*(288*q+8)+96,288*q+8,EPS2*Q),('runtime',48*q*5+6,48*q,p),('NB-log',248,8,NE),('NB-quarter',242,2,F(1,4)),('Lm',70,4,NE)]:
  ck('induction.'+nm,c+40*b,budget*U,'<');ck('induction.slope',b,c+40*b,'<=')
 batch('PXG.',[('log-base',62,64,'<='),('log-slope',2,62,'<='),('mh',4*(256+3)+F(65,2),1199,'<'),('dimension',1199*2**20+9,2048*2**20,'<='),('product',2048*256,2**19,'==')])
 for ss in sorted({s for k in range(40,81) for s in (2**k,2**(k+1)-1)}|{U+j for j in range(8)}):
  k=ss.bit_length()-1;s=8*ss;ellc=lg(s+255);rf=ss+8*ellc;ab=rf+AB;u1=2**(90+lg(ss));rho=ab+2*lg(u1)+2
  walk=s+40*rf;tag=3*walk
  batch('AFexec.',[('rho',rho,rf+185+2*lg(ss),'=='),('rho-dyadic',rho-ss,10*k+227,'<='),('rho-bound',rho,(1+EB)*ss,'<='),('log',ellc,k+5,'<='),('walk',16-8*ellc,0,'<'),('parity',walk%2,0,'==')])
  for L in (ceil(Q*ss),ceil(Q*(ss+1))-1):
   d=delta*L;kk=ceil(aa*d+bb*rho)
   for shift in (2,3):
    rho0=ab+L+shift+lg(kk)+2;b0=ceil(a0*rho0);u0=b0-rho0;b=ceil(gb*rho);uu=b-rho
    v0=b0*ceil(F(kk,u0));v1=b*ceil(F(kk,uu));cc=v1-kk;v2=b*ceil(F(cc,uu));D=b0+2*b
    ee=6*tag+rf+33+lg(C);om=ee+lg(26*t)+2*L+42
    batch('AFexec.',[('floor',floor(L/Q),ss,'=='),('rho0',rho0-(Q+1)*ss,9*k+Q+97,'<='),('rho0-bound',rho0,(Q+1+EB)*ss,'<='),('N0',rho0,N0,'>='),('u0',u0,1,'>='),('ratio0',F(b0,u0),N0+1,'<='),('u',uu,1,'>='),('ratio',F(b,uu),xi,'<='),('Kbound',kk,kb*ss,'<='),('Klog',lg(kk),48+lg(ss),'<='),('v0',v0,u1,'<='),('v1',v1,u1,'<='),('c',cc,1,'>='),('out',v1-v2,d,'>='),('rank-bound',v1-v2,gg*kk-xi*b,'>='),('seed',D,s,'<='),('source',kk,bh*L,'<='),('tag',tag,2*L,'>='),('runtime',tag,2*p*L/Q,'<='),('unit',om,co*L,'<='),('tail',8+2*ee,T*L,'<='),('combined',4*C*d+ko*om+8+2*ee,H*L,'<=')])
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
  ww=max(2**20,rho);ell=ceil(F(zz,2**hh))
  slots=[ceil((rr+1)*rho) for rr in vec];kk=ell+2*aa
  ck('PT.leaf-threshold',kk,F(zz,2**hh)+2*ww-11,'<=')
  ck('PT.expansion-guard',4*(4*nn)**2,2**(rho-aa),'<=')
  ck('PT.a-unit',aa,ww-6,'<=')
  for ix in range(hh-1,-1,-1):
   nxt=kk;bb=ceil(F(2*(nxt+aa+1),slots[ix]-2*rho));mm=ceil(F(bb*slots[ix],2));kk=mm+nxt+aa
   ck('PX.capacity',kk,bb*(slots[ix]-rho),'<');ck('PX.count',bb,kk,'<=')
   ck('PX.split',kk,mm+nxt+aa,'==');ck('PX.child',mm,nxt,'>')
   if ix==hh-1:
    leaf=mm;kap=(vec[-1]+1)/(vec[-1]-1)
    ck('PT.terminal',leaf,kap*(ceil(F(zz,2**hh))+3*aa+1)+F(slots[-1],2)+1,'<')
  ag=nm not in ('Ext0','Ext2')
  if ag: tail=2*(lg(leaf*ceil(F(zz,2**hh))+1)+ceil(F(zz,2**hh))+2*aa)
  else:tail=leaf
  dd=sum(slots)+tail;en=env(vec,ag)
  kap=(vec[-1]+1)/(vec[-1]-1)
  ck('PT.terminal-unit',kap*(ell+3*aa+1)+F(slots[-1],2)+1,
   kap*F(zz,2**hh)+(3*kap+(vec[-1]+1)/2)*ww+F(3,2)-16*kap,'<=')
  ck('PT.seed-slices',sum(slots),sum(rr+1 for rr in vec)*ww+hh,'<=')
  if ag and zz<=256*2**hh*ww:
   batch('PXG.',[('ell',ell,256*ww,'<='),('kap',kap,4,'<='),('tau',(vec[-1]+1)/(2*vec[-1]),F(4,5),'<='),
    ('Sh',slots[-1],65*ww+1,'<='),('mh',leaf,1199*ww+9,'<='),('mh2',1199*ww+9,2048*ww,'<='),
    ('lg',lg(leaf*ell+1),20+2*lg(ww),'<='),('lg-budget',20+2*lg(ww),F(ww,2**14),'<=')])
  ck('PX.threshold',kk,en['lam']*zz+en['g']*ww,'<=')
  if not ag or zz<=256*2**hh*ww:ck('PX.seed',dd,en['c']*zz+en['pi']*ww,'<=')
  ck('PX.padding',kk,nn,'<=');return dd,kk
 for nm in v:
  for nn,zz,ee in ((2**256,1,1),(2**512,100000,100),(2**2048,1,2**22)):
   px(nm,nn,zz,ee)
 # NB's literal actual algorithm at the first permitted outer unit.
 L=ceil(Q*U);rf=U+8*lg(8*U+255);tag=3*(8*U+40*rf)
 ee0=6*tag+rf+33+lg(C)+lg(26*t);om=ee0+2*L+42
 nn=do*om;zz=DX*om;adv=tag+2+lg(C);ee=ee0+4
 wn=max(2**20,ee+44+2*lg(4*nn));wf=max(2**20,ee+2*lg(4*nn)+32,ceil(F(zz,128)))
 m1,kf1=ceil(i['d']*wf),ceil(i['k']*wf);m0,k2=px('Ext2',nn,m1,ee)
 sig=ee+4;rr1=ceil((m0+x['pi']*wn+F(38,7)*ee+adv+53)/(theta-7*x['c']))
 ss=floor((eta*rr1-17*sig-8)/7);d01,k01=px('Ext0',nn,7*rr1,ee)
 el=ceil(F(d01+ee-1,2))+ee;d11=2*(lg(rr1*m0+1)+m0+2*el);d2=d01+d11+m1+k2+ee
 batch('NBexec.',[('I1',ss,m0+2*el,'>='),('advice',ss,adv+40,'>='),('density',F(adv+sig,rr1+adv),eta/2,'<='),('source',rr1-sig,(1-eta)*(rr1+adv+1)+1,'>='),('capacity',floor(F(8,7)*eta*(rr1+adv))-8,8*ss+17*sig,'>='),('deficit',eta*(rr1+adv),64*sig,'>='),('prime',rr1+adv,3200000,'>'),('dC',d01+d11+d2,sc*om,'<='),('kC',max(k01,kf1+14*rr1+zz+2*ee),kc*om,'<=')])
 # Kernel-duality examples include non-surjective matrices and dependent rows.
 def rank(rows):
  piv={}
  for row in rows:
   while row:
    k=row.bit_length()
    if k not in piv:piv[k]=row;break
    row^=piv[k]
  return len(piv)
 for dim in range(1,4):
  spcs={frozenset([0])}
  for _ in range(dim):
   spcs|={frozenset(set(w)|{x^b for x in w}) for w in list(spcs) for b in range(2**dim)}
  for code in range(2**(dim*dim)):
   rows=[(code>>(dim*j))%(2**dim) for j in range(dim)]
   ker={x for x in range(2**dim) if all((x&r).bit_count()%2==0 for r in rows)}
   for spc in spcs:
    perp={x for x in range(2**dim) if all((x&w).bit_count()%2==0 for w in spc)}
    if ker&perp!={0}:continue
    bas=[]
    for r in sorted(ker):
     if rank(bas+[r])>len(bas):bas.append(r)
    for d in range(len(bas)+1):
     im={sum(((x&r).bit_count()%2)<<j for j,r in enumerate(bas[:d])) for x in spc}
     ck('duality.rank',len(im),2**d)
 for dim in range(1,20):
  for rr in range(dim):ck('duality.rank-loss',F(2**rr,2**dim),F(1,2),'<=')
 ck('result.H',H,73355899284943);ck('result.K',K,293423597139772)
 ck('result.raw',raw,73355899284943);ck('result.base',bh,73355899234929)
 ck('result.improvement',K,355336740966308,'<')
 if TEXT:
  for nm,tag in [('I','MI'),('F','MF'),('G','MG'),('J','MJ'),('Ext0','NB0'),('Ext2','NB2'),('L0','L0'),('Lr','Lr'),('E','E')]:pt(V,' v'+tag+'=('+VE[nm]+')')
  for f in (P,V,S):pt(f,'Q_S='+str(Q))
  for tok in ('TR=TR0+H(t+1)(2dy+2dx)-(t+1)dx','T=TR+TL+(t+1)dx'):pt(W,tok)
  for f in (P,W):pt(f,'Ext0=PXT(vNB0)');pt(f,'Ext2=PXT(vNB2)')
  for f,tok in [(P,f'dy={DY}Om, dx={DX}Om, r={RATE}Om'),(P,f's_C={sc}Om'),(P,f'k_Cslot={kc}Om'),(P,f'd_O={do}Om, k_O={ko}Om'),(P,f'={tr}Om'),(P,f'={tl}Om'),(P,f'K=4H={K}'),(V,f'B_B=ceil((xi_B*gamma_B+xi_B/2^40)/G_B)={Bbase}')]:pt(f,tok)
  sp=TEXT[S].split('# SP. Finite rings')[1]
  heads=re.findall(r'^#+ (SP\.[0-9.]+)',TEXT[S],re.M);tc('refs.headings',len(heads),len(set(heads)))
  labs=re.findall(r'^ +.*\((SP-[\w]+)\)\s*$',sp,re.M)
  tc('refs.labels',set(labs),set('SP-1 SP-2 SP-3 SP-4 SP-5 SP-5s SP-Mi0 SP-Mi1 SP-Mi2 SP-Mi3 SP-Mi4 SP-17 SP-18 SP-19m SP-21m SP-S2 SP-M SP-I'.split()))
  tc('refs.unique',len(labs),len(set(labs)))
  for lab in re.findall(r'\((SP-[\w]+)\)',sp):tc('refs.'+lab,lab in labs)
  for f in TEXT:
   tc('refs.removed','sum-product-and-ordinary-condensation.md' not in TEXT[f])
   for sec in re.findall(r'Section (SP\.\d+(?:\.\d+)?)',TEXT[f]):tc('refs.section',bool(re.search(r'^#+ '+re.escape(sec)+r'(?:[. ]|$)',TEXT[S],re.M)))
  tc('refs.journal','Mathematics 4(2)' in sp);tc('refs.stale','theoremI' in sp,False)
  files=[p for p in ROOT.iterdir() if p.is_file()];tc('limits.count',len(files),8,'<=')
  for p in files:tc('limits.file',p.stat().st_size,65536,'<=')
  tc('limits.total',sum(p.stat().st_size for p in files),262144,'<=')
  for f,h in HASHES.items():tc('hash.'+f,hashlib.sha256((ROOT/f).read_bytes()).hexdigest(),h)
  bj=json.loads((ROOT/'bound.json').read_text());sj=json.loads((ROOT/'submission.json').read_text())
  decl=int(bj['guarantee']['remainder']['terms'][0]['coefficient']['numerator']);tc('meta.K',decl,K)
  tc('meta.experiment',[e.get('source') for e in sj['experiments']],['verify.py'])
  tc('meta.evidence',all((ROOT/f).is_file() for claim in sj['claims'] for f in claim['evidence']))
  for claim in sj['claims']:tc('frozen.declaration',json.loads(claim['statement'].split('Declaration: ')[1]),bj)
 return dict(all_passed=not fail,all_checks_passed=not fail and not meta,checks=ct,arithmetic_checks=ac,text_checks_run=bool(TEXT),skipped_reason=None if TEXT else 'Package absent or incomplete; arithmetic still executed',H=str(H),K=str(K),H_raw=str(raw),H_base=str(bh),failures=fail,metadata_consistent=(not meta) if TEXT else None,metadata_failures=meta,groups=gs,scope='Exact arithmetic/integrity. Markdown supplies universal proofs; metadata is separate.')
HASHES={'proof.md': '57b1ec450744fd4699b114c59ce55ce9b6609e33edf96c41b672274041fd3ee2', 'bound.json': '8b0eee9ab6d36d13b7322d0bf857218e2a39922520dd02fa0bd51ceecd3d71a7', 'sampler-and-disperser.md': 'ac3f554c44e04325f428a361ab1566ee7a5181fc90a4bcd392ff272a2f687f97', 'restricted-nonmalleability-and-breaking.md': 'a4f77f941438492ae58003ef3dca568341d9f62a52b23606f6c508395e6a8a6f', 'linear-and-bilinear-primitives.md': 'b8c2b2986bd8069ccbca7d9e6102b9a23eb18207c6942edff2c44e197073f28e', 'affine-and-standard-wrappers.md': '7998a24fffb71f6924f511b3e5c494beebfd8488289df0c84f6efcf7c7ba2bb1', 'submission.json': 'f0c17d1db6aeced2d7556b4539f428c30b0eebb60672c4b382a6cf0e117f4c8b'}
if __name__=='__main__':
 if json.load(sys.stdin)!={}:raise ValueError('Expected {}')
 roots=(Path(__file__).resolve().parent,Path('submission'),Path('/snapshot/files/submission'))
 ROOT=next((p for p in roots if all((p/f).is_file() for f in (P,W,V,S,R,'bound.json','submission.json'))),None)
 if ROOT:TEXT={f:(ROOT/f).read_text() for f in (P,W,V,S,R)}
 res=run();out=json.dumps(res,sort_keys=True);assert len(out.encode())<65536
 print(out);raise SystemExit(0 if res['all_passed'] else 1)
