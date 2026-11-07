from fractions import Fraction as F
from math import ceil,floor,isqrt,comb
from pathlib import Path
import json,sys,re,hashlib
EPS=F(1,2**20);NU=F(1,2**80); NE=F(1,10**18); EB=F(1,10**18);ER=F(1,2**80);N0=2**60;AB=3
SL={'a': '4866429720000000001/200000000000000000', 's': '4545765360000000001/200000000000000000', 's0': '6035509597479093741/250000000000000000', 'hc': '71267993438661999391/250000000000000000', 'd': '2022585474044559308479/1000000000000000000', 'k': '78206013596976137519594059/1000000000000000000'}
VN={'I': '11888889,12000000,12000000,12000000,12988999,13000000,13000000', 'F': '10788990,12990000,13000000,13000000,13000000,13000000,13000000', 'G': '22986109,22641416', 'J': '35987358,33689994', 'Ext0': '254207875,237000000,227000000,225000000,225000000,225000000,225000000,225000000,225000000,225000000,225000000,225000000,225000000,225000000,225000000,224000000,215000000,204999990,175000000,157792025', 'Ext2': '42900000,43889999,43900000,44100000,43800000,42600000,64587881', 'L0': '38043035,39043128,29053087,29053963,29063718,29063678,30060594,30051514,30024294,29943548,29704763,29004651,27099079', 'Lr': '45000000,45000000,45000000,45000000,45000000,45000000,45000000,45000000,45000000', 'E': '12400000,12400000,12400000,12400000,12400000,12400000'}
VE={k:','.join(str(F(int(x),10**7)) for x in s.split(',')) for k,s in VN.items()}
QS=F('29970175275836701176998397/10000000000000000000000000');QC=20;DX=80 # linear-and-bilinear-primitives.md:3; proof.md:745
P='proof.md';W='affine-and-standard-wrappers.md';V='linear-and-bilinear-primitives.md'
S='sampler-and-disperser.md';R='restricted-nonmalleability-and-breaking.md'
TAG=dict(I='MI',F='MF',G='MG',J='MJ',Ext0='NB0',Ext2='NB2',L0='L0',Lr='Lr',E='E')
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
  pi=4+sum(r+1 for r in vec) if ag else sum(r+1 for r in vec)+3*kap+(r+1)/2+max(0,h+F(3,2)-16*kap)*EPS)
def et(e):return ','.join(k+'='+str(e[k]) for k in ('c','pi','lam','g'))
def walkextra(s):return ((2**16*25*(s+4)**40).bit_length()-1)//7+1
def slot(x):return F(ceil(x*10**18)+2,10**18)
def run():
 v={k:list(map(F,s.split(','))) for k,s in VE.items()};i={k:F(s) for k,s in SL.items()}
 ev={k:env(a,k not in ('Ext0','Ext2')) for k,a in v.items()}
 I,FF,G,J=[ev[k] for k in ('I','F','G','J')];x,y=ev['Ext0'],ev['Ext2']
 DR=max(0,len(v['Ext0'])+F(3,2)-16*x['c']*2**len(v['Ext0']));PB=x['pi']-DR*EPS
 eta=F(1208925819614629174706175,4274761698157328761761038896);theta=eta/6;gamma=F(241785163922925834941235,19342813113834066795298816);U=2**40
 alpha=F(433)/(433+2*QS)+NE
 df=i['d']*alpha+NE;kf=i['k']*alpha+NE
 M0=y['c']*df+y['pi']*alpha+NE;K2=y['lam']*df+y['g']*alpha+NE
 rr=(M0+PB*alpha+F(11,2)*alpha+(53+DR)*NE)/(theta-6*x['c'])+NE
 d0=6*x['c']*rr+PB*alpha+(DR+1)*NE;d1=2*(M0+d0+3*alpha)+NE
 k00=6*x['lam']*rr+x['g']*alpha;kother=kf+12*rr+DX+2*alpha
 dc=2*d0+2*d1+df+K2+alpha;sc=slot(dc);kc=slot(max(k00,kother))
 q=QC;C=4*(q+1);t=4*C-1;j=lg(t);Q=QS;beta=F(2*(q+8),7);p=3*(beta+2*q)
 EE,RR,ZZ=[ev[k] for k in ('E','Lr','L0')]
 DY=ceil(EE['c']*DX+EE['pi']);RATE=(t+2)*DY+3
 s0=slot(ZZ['c']*sc+ZZ['pi']);tr=(t+1)*(s0+sc+2*j*(DY+DX));tl=(t+1)*(sc+2*j*DY)
 rows=dict(dy_E=EE['c']*DX+EE['pi'],dx_Lr=RR['c']*RATE+RR['pi'],dx_Lm=2*DY+4+NE,s0_L0=ZZ['c']*sc+ZZ['pi'],r_3a=(t+2)*DY+3,r_dy=DY,r_Km=DY+2,
 d_3c=(t+1)*(s0+sc)+kc+alpha,d_3d=tr+EE['lam']*DX+EE['g']+t*DX+alpha,d_s0=s0,d_KE=EE['lam']*DX+EE['g'],d_kC=kc,
 k_3b=tl+RR['lam']*RATE+RR['g']+t*RATE+alpha,k_K0=ZZ['lam']*sc+ZZ['g'])
 do=slot(max(val for key,val in rows.items() if key.startswith('d_')));ko=slot(max(val for key,val in rows.items() if key.startswith('k_')))
 co=(3*p+1)/Q+2+F(1,10**15);delta=ceil(do*co);kappa=ceil(ko*co);T=ceil((6*p+2)/Q+EPS)
 a0=F(N0+1,N0);gb=(beta-a0*(Q+1+EB)-3*ER)/(2*(1+EB));xi=gb/(gb-1);gg=1-(xi-1)**2;aa=1/gg;bb=ceil((xi*gb+xi/U)/gg)
 Bbase=bb;raw=4*C*delta+kappa+T;bh=ceil(aa*delta+bb*(1+EB)/Q+ER);H=max(raw,bh);K=4*H
 kb=aa*delta*Q*(1+F(1,U))+bb*(1+EB)+ER;b0c=a0*(Q+1+EB)+ER
 IR=lambda I,FF,G,J:[('aJ','a',J['c']*(i['s']+NE)+J['pi']),('aF','a',FF['c']*(i['hc']+NE)+FF['pi']),('sG','s',G['c']*(i['a']+NE)+G['pi']),('s0I','s0',I['c']*(i['hc']+NE)+I['pi']),('hJ','hc',(9+J['lam'])*(i['s']+NE)+J['g']+1),('hcap','hc',F(128)),('h>=s','hc',i['s']+NE),('dG','d',13*(i['s0']+NE)+(65+G['lam'])*(i['a']+NE)+G['g']+1),('d>=s0','d',i['s0']+NE),('kF','k',(38+FF['lam'])*(i['hc']+NE)+39*(i['s']+NE)+FF['g']+1),('kI','k',I['lam']*(i['hc']+NE)+I['g']),('k>=d','k',i['d']+NE)];ir=IR(I,FF,G,J)
 for nm,ky,val in ir:
  ck('inner.'+nm,i[ky],val,'>=');pt(W,f'{nm}: [{ceil(val/NE)*NE}]w <= ({i[ky]})w.')
 for ky,val in i.items():pt(W,f'{ky}=ceil(({val})w).')
 for ky,val in rows.items():
  lim=dict(dy=DY,dx=DX,r=RATE,s0=s0,d=do,k=ko)[ky.split('_')[0]]
  ck('outer.'+ky,lim,(F(ceil(val*10**18),10**18)+NE if ky.split('_')[0] in ('s0','d','k') else val),'>=');pt(P,f'{ky}: '+('U18(rhs)+eps_N' if ky.split('_')[0] in ('s0','d','k') else 'rhs')+f' <= {lim}.')
 for nm,vec in v.items():
  en=ev[nm];h=len(vec)
  for rate in vec:
   ck('PT.rate',rate,1,'>')
   sh=ceil(7*(rate+1));nn=ceil(F(3*(sh-14),2))-2
   bc=ceil(F(2*(nn+2),sh-14));mm=ceil(F(bc*sh,2));kk=mm+nn+1
   ck('PT.boundary',kk,bc*(sh-7),'<')
  f=W if nm in ('I','F','G','J') else P
  if nm in ('I','F','G','J'):pt(W,'Section7g recurrences at vMI,vMF,vMG,vMJ')


  if nm not in ('Ext0','Ext2'):
   out=dict(I=i['hc']+EPS,F=i['hc']+EPS,G=i['a']+EPS,J=i['s']+EPS,L0=sc,Lr=RATE,E=DX)[nm]
   batch('PXG.',[('min',vec[-1],F(16,15),'>='),('max',vec[-1],64,'<='),('output',out,256*2**h,'<='),('error',3*2**h-2,2**(h+4),'<')])
 pt(P,'bases 1067<')
 pt(W,'p0='+str(PB));pt(W,'D0='+str(DR))
 pt(W,'c0='+str(x['c']));pt(W,'pi0='+str(x['pi']))
 pt(W,'6c0='+str(6*x['c'])+'<theta='+str(theta));pt(P,f'delta_O=ceil({do}c_O)={delta}');pt(P,f'kappa_O=ceil({ko}c_O)={kappa}');pt(P,'D_F=alpha*d_inner+eps_N='+str(df));pt(P,'K_F=alpha*k_inner+eps_N='+str(kf));pt(P,str(co))
 for val in (gb,xi,gg,kb):pt(V,str(val))
 cl=[ceil(z) for z in (M0,K2,rr,d0,d1,dc)];pt(P,'ceil(M0)=%d, ceil(K2bar)=%d,\n ceil(Rnew)=%d, ceil(d0bar)=%d,\n ceil(d1bar)=%d, ceil(DNB)=%d.'%tuple(cl))
 pt(P,'KC='+str(kc));pt(P,'SC=Slot(DNB)='+str(sc))
 # AE: one-sided 10^-30 enclosures reproduce every ceiling and slot.
 eu={k:{c:F(ceil(e[c]*10**40),10**40) for c in e} for k,e in ev.items()}
 for c in ('c','pi','lam','g'):
  for k in ev:ck('AE.enclosure',0<=eu[k][c]-ev[k][c]<F(1,10**40),True)
 xu,yu=eu['Ext0'],eu['Ext2'];M0u=yu['c']*df+yu['pi']*alpha+NE;K2u=yu['lam']*df+yu['g']*alpha+NE
 ck('AE.denominator',theta,6*xu['c'],'>');ru=(M0u+(xu['pi']-DR*EPS)*alpha+F(11,2)*alpha+(53+DR)*NE)/(theta-6*xu['c'])+NE
 d0u=6*xu['c']*ru+(xu['pi']-DR*EPS)*alpha+(DR+1)*NE;d1u=2*(M0u+d0u+3*alpha)+NE;dcu=2*d0u+2*d1u+df+K2u+alpha
 aec=[ceil(z) for z in (M0u,K2u,ru,d0u,d1u,dcu,6*xu['lam']*ru+xu['g']*alpha,kf+12*ru+DX+2*alpha)]
 ck('AE.seed',dcu,sc-NE,'<=');ck('AE.source',max(6*xu['lam']*ru+xu['g']*alpha,kf+12*ru+DX+2*alpha),kc-NE,'<=');
 Iu,Fu,Gu,Ju=[eu[k] for k in ('I','F','G','J')]
 for nm,ky,val in IR(Iu,Fu,Gu,Ju):ck('AE.inner.'+nm,i[ky],val,'>=')
 Eu,Ru,Zu=eu['E'],eu['Lr'],eu['L0']
 for nm,lim,val in [('dyE',DY,Eu['c']*DX+Eu['pi']),('dxLr',DX,Ru['c']*RATE+Ru['pi']),('s0L0',s0,Zu['c']*sc+Zu['pi']),('d3d',do,tr+Eu['lam']*DX+Eu['g']+t*DX+alpha),('dKE',do,Eu['lam']*DX+Eu['g']),('k3b',ko,tl+Ru['lam']*RATE+Ru['g']+t*RATE+alpha),('kK0',ko,Zu['lam']*sc+Zu['g'])]:ck('AE.outer.'+nm,lim,val+(NE if nm in ('s0L0','d3d','dKE','k3b','kK0') else 0),'>=')
 for val in (k00,kother):ck('NB.source',val,kc,'<=')
 ck('NB.seed',dc,sc,'<=');ck('NB.source-dimension',kc,do,'<=');ck('NB.seed-dimension',sc,do,'<=')
 batch('NB.',[('denominator',theta,6*x['c'],'>'),('p0',PB,16,'>='),('height',max(len(v['Ext0']),len(v['Ext2'])),38,'<='),('cap',DX,128,'<='),('advice',F(11,2)+F(53,U),28,'<='),('AG2-p0',6*PB,2826,'>='),('AG2-slope',160+F(600,2**80),2826,'<'),('AG2-const',6*16*44+318,600,'>'),('constant',ceil(2*(rr+2)*M0)+1,2**240,'<'),('Lm',RATE*DY+1,2**33,'<'),('heights',max(len(v[k]) for k in ('I','F','G','J'))+13,21,'<='),('outer-heights',max(len(v[k]) for k in ('L0','Lr','E'))+12,25,'<='),('NB-heights',max(len(v['Ext0']),len(v['Ext2']))+6,27,'<='),('prime',6*2**80/eta,2**90,'>'),('capacity',F(7,6)*6,7,'=='),('floor',53-F(34,3)-1,F(122,3),'=='),('advice40',F(122,3),40,'>'),('error',14,16,'<'),('final-error',15,16,'<')])
 for hh in (0,1,2,j):
  rt=s0+sc+DX;lf=sc;mx=0
  for _ in range(hh):
   rt+=DY;lf+=DY;mx=max(mx,rt);rt+=DX+DY;lf+=DY;mx=max(mx,rt);rt+=DX
  if hh:ck('causal.E',mx,s0+sc+2*hh*(DY+DX))
  ck('causal.full',lf+rt,s0+2*sc+DX+hh*(4*DY+2*DX))
 batch('AF.',[('ceil-round',1,ER*2**80,'<='),('gamma',gb,2,'>'),('rank',gg*aa,1,'=='),('rounding',gg*bb-xi*gb,xi/U,'>='),('B-minimal',gg*(bb-1)-xi*gb,xi/U,'<'),('seed',a0*(Q+1+EB)+2*gb*(1+EB)+3*ER,beta,'=='),('Kbits',kb,2**48,'<'),('pad0',(N0+1)*2**48+b0c,2**110,'<'),('pad1',xi*2**48+gb*(1+EB)+EB,2**110,'<'),('error',6,2**AB,'<'),('q',F(7,2)*beta-8,q,'=='),('p',p,3*Q,'>'),('H',H,2*delta,'>='),('runtime-power',ceil(2*p/Q),97,'=='),('incidence',48*U,134217728,'>')])
 for nm,c,b,budget in [('rho',267,10,EB),('rho0',Q+97,9,EB),('outer',10*(288*q+8)+160,2*(288*q+8),EPS*Q),('runtime',48*q*5+6,48*q,p),('NB-log',248,8,NE),('NB-quarter',242,2,F(1,4)),('Lm',70,4,NE)]:
  kk=80 if nm in ('rho','rho0','NB-log','Lm') else 40;ck('induction.'+nm,c+kk*b,budget*2**kk,'<');ck('induction.slope',b,c+kk*b,'<=')
 ck('induction.Om-base',5768*85+96,F(1,10**15)*Q*2**80,'<');ck('induction.Om-slope',5768,5768*85+96,'<=')
 ck('PXG.legal-small',F(1,8192)+15*EPS,2,'<')
 batch('PXG.',[('log-base',64,64,'<='),('log-slope',2,64,'<='),('mh',31*(256+3)+F(65,2),8062,'<'),('dimension',8062*2**20+9,8192*2**20,'<='),('product',8192*256,2**21,'==')])
 batch('local.',[('alpha',alpha,1,'<'),('output',F(DX,128)+NE,alpha,'<='),('source',do,2**30,'<'),('Omega',co*Q*(1+F(1,2**80)),1024,'<'),('unit-base',490628,NE*(433+2*Q)*2**80,'<'),('unit-slope',5770,490628,'<='),('hash-base',804,NE*2**80,'<'),('hash-slope',4,804,'<='),('round',1,NE*2**80,'<='),('advice-charge',F(11,2),F(11,2)),('AG1-coef',9*F(1,3),1,'>='),('AG2-coef',120*F(1,3)+120,160,'==')])
 pt(P,'5770*80+29028=490628');pt(V,'16/15<=Rh<=64')
 for ss in sorted({s for k in range(80,121) for s in (2**k,2**(k+1)-1)}|{2**80+j for j in range(8)}):
  k=ss.bit_length()-1;s=8*ss;ellc=lg(s+4);rf=ss+walkextra(s);ab=rf+AB;u1=2**(110+lg(ss));rho=ab+2*lg(u1)+2
  walk=s+40*rf;tag=3*walk
  batch('AFexec.',[('rho',rho,rf+225+2*lg(ss),'=='),('rho-dyadic',rho-ss,10*k+267,'<='),('rho-bound',rho,(1+EB)*ss,'<='),('log',ellc,k+5,'<='),('walk',2**(7*(rf-ss)),2**16*25*(s+4)**40,'>'),('parity',walk%2,0,'==')])
  for L in (ceil(Q*ss),ceil(Q*(ss+1))-1):
   ee=3*tag+rf+42+lg(C);om=ee+lg(26*t)+2*L+42;d=floor(do*om);kk=ceil(aa*d+bb*rho)
   eeN=ee+lg(26*t)+4;np=d;adv=tag+2+lg(C)
   wn=max(2**80,eeN+44+2*lg(4*np));rlo=ceil((PB*max(2**80,eeN+44+2*lg(4*np))+DR+F(11,2)*eeN+53)/(theta-6*x['c']));wf=max(2**80,eeN+2*lg(4*np)+32,ceil(F(DX*om,128)))
   batch('localexec.',[('Omega',om,1024*ss,'<='),('nlog',lg(4*np),42+lg(ss),'<='),('wN',wn,alpha*om,'<='),('wF',wf,alpha*om,'<='),('advice',3*adv,eeN,'<='),('advice-exact',3*adv,eeN-(rf+40),'=='),('AG1-term',PB*wn+3*eeN+F(122,3),adv+40,'>='),('AG2-term',6*PB*wn+33*eeN+318,120*(adv+eeN+4+1),'>='),('AG2-wN',2826*wn,160*wn+600,'>'),('AG1-exec',floor((eta*rlo-15*(eeN+4)-8)/6),adv+40,'>='),('AG2-exec',eta*rlo,120*(adv+eeN+4+1),'>='),('hash',480+4*lg(om),NE*om,'<=')])
   for shift in (2,3):
    rho0=ab+L+shift+lg(kk)+2;b0=ceil(a0*rho0);u0=b0-rho0;b=ceil(gb*rho);uu=b-rho
    v0=b0*ceil(F(kk,u0));v1=b*ceil(F(kk,uu));cc=v1-kk;v2=b*ceil(F(cc,uu));D=b0+2*b
    ee=3*tag+rf+42+lg(C);om=ee+lg(26*t)+2*L+42
    batch('AFexec.',[('floor',floor(L/Q),ss,'=='),('rho0',rho0-(Q+1)*ss,9*k+Q+97,'<='),('rho0-bound',rho0,(Q+1+EB)*ss,'<='),('N0',rho0,N0,'>='),('u0',u0,1,'>='),('ratio0',F(b0,u0),N0+1,'<='),('u',uu,1,'>='),('ratio',F(b,uu),xi,'<='),('Kbound',kk,kb*ss,'<='),('Klog',lg(kk),48+lg(ss),'<='),('v0',v0,u1,'<='),('v1',v1,u1,'<='),('c',cc,1,'>='),('out',v1-v2,d,'>='),('rank-bound',v1-v2,gg*kk-xi*b,'>='),('seed',D,s,'<='),('source',kk,bh*L,'<='),('tag',tag,2*L,'>='),('runtime',tag,2*p*L/Q,'<='),('unit',om,co*L,'<='),('tail',8+2*ee,T*L,'<='),('combined',4*C*d+ko*om+8+2*ee,H*L,'<=')])
 batch('',[('SP.gamma',F(1,4)-20*gamma,F(1,2**82),'=='),('SP.guard',37*2**82,2**90,'<'),
 ('SP.subset',(F(1,16)-gamma)*2**90,3,'>'),('SP.dimension',18*gamma+F(33,2**90),F(1,4),'<'),
 ('SP.power',F(7,3)-F(1,16)+F(1,4),F(121,48),'=='),('SP.eta',(F(221,5)+7*NU)*eta,gamma,'=='),('SP.cap',eta,F(1,3000),'<='),('SP.theta',6*theta,eta,'=='),('CZ.exponent',6*F(4,3)+6*(F(91,30)+NU)+18,F(221,5)+6*NU,'=='),('CZ.hered',(1000-(F(91,30)+NU)-3)*64/NU,6,'>'),('CZ.zero',1-27*eta,0,'>'),('CZ.T',2+F(7,3)*eta,4,'<='),('CZ.smoothing-exponent',F(7,6),1,'>='),
 ('SP.cardinal',5+3+1,9,'=='),('SP.sixteen',9,16,'<'),('SP.endpoint',F(15+4)*3/8,F(57,8),'=='),('SP.endpoint-loss',F(57,8),8,'<'),
 ('CZ.units',1-F(1021,3000),F(1979,3000),'=='),('CZ.plane',F(31,16)-F(1020,3000)-1,F(239,400),'=='),
 ('CZ.projection',F(31,16)-F(1041,3000)-F(5,4),F(681,2000),'=='),('CZ.energy',2**63+2**57,2**64,'<'),
 ('CZ.smoothing',30+3+2+10,45,'<='),('CZ.hit',35,36,'<'),('CZ.root16',16,2**(20*64),'<'),
 ('prime.anchor',F(32,5),F(99,16),'>'),('prime.square',631,26**2,'<'),
 ('native.union',32*2**16*2*2**20,2**42,'=='),('native.collision',F(2,2048)+F(1,4),F(257,1024),'=='),
 ('native.r',F(257,1020)+F(1,2**17),F(13,50),'<'),('native.mono',6-8*F(13,50),0,'>'),
 ('native.quartic',6*F(13,50)-4*F(13,50)**2-2,F(-444,625),'=='),('native.margin',F(-444,625)+F(1,1000)+F(4,1024),F(-1,2),'<='),
 ('native.g',F(4,2**4095),F(1,1000),'<'),('native.delete',F(1,2**27),F(1,10000**2),'<'),('native.perturb',F(432,10000),F(1,8),'<'),
 ('native.moment',3+F(3,1024),4,'<='),('native.success',1-2*F(2,5),F(1,5),'=='),
 ('F.entropy',2*13+12,38,'=='),('F.seed3',3*13,39,'=='),('F.seed',5*13,65,'=='),('F.labels',1+2*4,9,'=='),('F.error',76,128,'<'),('F.round',64,76,'<='),('native.H128',H*128,2**128,'<')])
 ck('SP.constants',(8**4,2**15,2**17*2**20,2**15*2**18),(4096,32768,2**37,2**33))
 ck('native.ratio',F(257,1024)*F(256,255),F(257,1020))
 ck('native.bias',F(1,2**18)+F(1,2**37),F(1,2**17),'<')
 ck('native.weight-count',sum(2**j for j in range(1,5)),32,'<=');ck('native.moment',32*F(4,2**17),F(1,1024),'==')
 ck('native.deletion',F(1,2**28)*F(256,255),F(1,2**27),'<')
 primes=[2,3,5,7,13,23,43,83,163,317,631]
 for ix,pp in enumerate(primes):
  ck('prime.trial',all(pp%d for d in range(2,isqrt(pp)+1)))
  if ix:ck('prime.chain',pp,2*primes[ix-1],'<')
 codes=['11000','10100','10010','01100','01010','00110','10001']
 pt(W,','.join(codes));ck('F.distinct',len(set(codes)),7)
 for aa in codes:
  first=aa.index('1');ck('F.prefix9',1+2*sum(bb[first]=='0' for bb in codes if bb!=aa),9,'<=')
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
 def px(nm,nn,zz,ee,unit=0):
  vec=v[nm];hh=len(vec);aa=ee+hh+4;rho=aa+2*lg(4*nn)+2
  ww=max(unit,2**20,rho);ell=ceil(F(zz,2**hh))
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
   batch('PXG.',[('ell',ell,256*ww,'<='),('kap',kap,31,'<='),('tau',(vec[-1]+1)/(2*vec[-1]),F(31,32),'<='),
    ('Sh',slots[-1],65*ww+1,'<='),('mh',leaf,8062*ww+9,'<='),('mh2',8062*ww+9,8192*ww,'<='),
    ('lg',lg(leaf*ell+1),22+2*lg(ww),'<='),('lg-budget',22+2*lg(ww),F(ww,2**14),'<=')])
  ck('PX.threshold',kk,en['lam']*zz+en['g']*ww,'<=')
  if not ag or zz<=256*2**hh*ww:ck('PX.seed',dd,en['c']*zz+en['pi']*ww,'<=')
  if nm=='Ext0':ck('PX.finite-seed',dd,en['c']*zz+PB*ww+DR,'<=')
  ck('PX.padding',kk,nn,'<=');return dd,kk
 for nm in v:
  for nn,zz,ee in ((2**256,1,1),(2**512,100000,100),(2**2048,1,2**22)):
   px(nm,nn,zz,ee)
 for nm in ('I','F','G','J','L0','Lr','E'):
  hh=len(v[nm]);px(nm,2**48,1,1,2**80);px(nm,2**256,256*2**hh*2**20,2**20-hh-522)
 # NB's literal algorithm at S=2^80, inside the proved local-unit range.
 UU=2**80;L=ceil(Q*UU);rf=UU+walkextra(8*UU);tag=3*(8*UU+40*rf)
 ee0=3*tag+rf+42+lg(C)+lg(26*t);om=ee0+2*L+42
 nn=floor(do*om);zz=DX*om;adv=tag+2+lg(C);ee=ee0+4
 wn=max(2**80,ee+44+2*lg(4*nn));wf=max(2**80,ee+2*lg(4*nn)+32,ceil(F(zz,128)))
 m1,kf1=ceil(i['d']*wf),ceil(i['k']*wf);m0,k2=px('Ext2',nn,m1,ee)
 sig=ee+4;rr1=ceil((m0+PB*wn+DR+F(11,2)*ee+53)/(theta-6*x['c']))
 ss=floor((eta*rr1-15*sig-8)/6);d01,k01=px('Ext0',nn,6*rr1,ee)
 el=ceil(F(d01+ee-1,2))+ee;d11=2*(lg(2*(rr1+adv+1)*m0+1)+m0+2*el);d2=d01+d11+m1+k2+ee
 batch('NBexec.',[('I1',ss,m0+2*el,'>='),('advice',ss,adv+40,'>='),('density',F(adv+sig,rr1+adv),eta/120,'<='),('AG1',ss,adv+40,'>='),('AG2',eta*rr1,120*(adv+sig+1),'>='),('AG1-margin',ss,1884*(adv+40),'>'),('AG2-margin',10*eta*(rr1+adv),236*120*(adv+sig),'>'),('source',rr1-sig,(1-eta/60)*(rr1+adv+1)+1,'>='),('capacity',F(7,6)*eta*(rr1+adv),7*ss+15*sig,'>='),('density-error',eta*(rr1+adv),6*ss+15*sig+8,'>='),('prime',rr1+adv,2**90,'>'),('dC',d01+d11+d2,floor(sc*om),'<='),('kC',max(k01,kf1+12*rr1+zz+2*ee),floor(kc*om),'<=')])
 ck('refine.log-slope',361*F(40,7),F(14440,7),'==');ck('refine.Econst',361*F(253,21)+49,F(92362,21),'==')
 ck('refine.Oconst',361*F(253,21)+105,F(93538,21),'==');ck('refine.tailconst',2*F(92362,21)+8,F(184892,21),'==')
 UO=4*C*do+ko;LA=433/Q+2;A=UO*LA+866/Q;B=max(F(14440,7)*(UO+2),do*F(14440,7)/gg+F(54,7)*Bbase);CC=ceil(max(F(93538,21)*UO+F(184892,21),do*F(93538,21)/gg+F(9907,42)*Bbase+1))
 ck('refine.Qpower',Q**12,2**19,'>');ck('refine.logmargin',(8/Q+F(4,2**80))**12,2**17,'<');
 batch('refine.',[('source',A,do*LA/gg+Bbase/Q,'>='),('log',B,do*F(14440,7)/gg+F(54,7)*Bbase,'>='),('constant',CC,do*F(93538,21)/gg+F(9907,42)*Bbase+1,'>='),('leading',A,H-1,'<='),('guard-base',B*81+CC+1,2**80,'<'),('guard-slope',B,B*81+CC+1,'<=')])
 ck('result.H',H,11972535894429);ck('result.K',K,47890143577716)
 ck('result.raw',raw,11972535894429);ck('result.base',bh,11972535894429)
 ck('result.improvement',4*A,F('30174155329677423964304027159190216193690196821309827327/624378650877472663098703562500000000000000'),'<');ck('result.loss',[4*A,4*B,4*(CC+1)],[F('59803166541457854438806728560915898445486487090621074199/1248757303159862549041599875000000000000000'),F('9365125743275850011523800440169425204763237471148364360573127285194332394893541644304633832958258860865141378159087/13885696126015665580455349976853313907393050325362870163359438155420808470158040320158750000000000000'),1456282357241572])
 if TEXT:
  for key,val in [('c1',4*A),('c2',4*B),('c0',4*(CC+1))]:pt(P,key+'='+str(val))
  for nm,tag in [('I','MI'),('F','MF'),('G','MG'),('J','MJ'),('Ext0','NB0'),('Ext2','NB2'),('L0','L0'),('Lr','Lr'),('E','E')]:pt(V,' v'+tag+'=('+VE[nm]+')')
  for f in (P,V,S):pt(f,'Q_S='+str(Q))
  for tok in ('TR=TR0+H(t+1)(2dy+2dx)-(t+1)dx','T=TR+TL+(t+1)dx'):pt(W,tok)
  for f,tok in [(P,f'dy={DY}Om, dx={DX}Om, r={RATE}Om'),(P,f's_C=floor(({sc})Om)'),(P,f'k_Cslot=floor(({kc})Om)'),(P,f'd_O=floor(({do})Om), k_O=floor(({ko})Om)'),(P,f'={tr}Om'),(P,f'={tl}Om'),(P,f'K=4H={K}'),(V,f'B_B=ceil((xi_B*gamma_B+xi_B/2^40)/G_B)={Bbase}')]:pt(f,tok)
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
  terms=bj['guarantee']['remainder']['terms'];declared=[F(t['coefficient']['numerator'])/F(t['coefficient']['denominator']) for t in terms];new=[4*A,4*B,4*(CC+1)]
  tc('meta.declared',declared==new)
  tc('meta.experiment',[e.get('source') for e in sj['experiments']],['verify.py'])
  tc('meta.evidence',all((ROOT/f).is_file() for claim in sj['claims'] for f in claim['evidence']))
  for claim in sj['claims']:tc('frozen.declaration',json.loads(claim['statement'].split('Declaration: ')[1]),bj)
 return dict(refined_loss=[str(4*A),str(4*B),str(4*(CC+1))],all_passed=not fail,all_checks_passed=not fail and not meta,checks=ct,arithmetic_checks=ac,text_checks_run=bool(TEXT),skipped_reason=None if TEXT else 'no package; arithmetic ran',H=str(H),K=str(K),H_raw=str(raw),H_base=str(bh),failures=fail,metadata_consistent=(not meta) if TEXT else None,metadata_failures=meta,groups=gs,scope='exact arithmetic+integrity')
HASHES={'proof.md': '14b503ce82de35a2595b10a370ae2f8bb6cdbb14beb6ca54e2afffe35ffbaad0', 'bound.json': '4d59046a6a9489d435f26781d637cf65b5011038e2ae9d7eb6852fa339dc245a', 'sampler-and-disperser.md': 'b31ed2048a2049bfb26a6d6c9cd07d1913652575a8127a8740df7b38de344ecc', 'restricted-nonmalleability-and-breaking.md': '059ea11be13a948e9e26777be3fb80fdc04ace967ee9939e991a8d081dfa89ec', 'linear-and-bilinear-primitives.md': '66bca252fad1c314850dcdd20b34ad19e038aae4653548c3acd3f24e9d2af09e', 'affine-and-standard-wrappers.md': 'ea7a968adbf948fff2705aaf3aca32ec8e14275aa30294b1fc50f8faaeeeee5b', 'submission.json': '1dbff83eeaf4d4ed2dc1d3185b2b19a273d931b7c6776506315e5d2c44deb6c7'}
if __name__=='__main__':
 if json.load(sys.stdin)!={}:raise ValueError('Expected {}')
 roots=(Path(__file__).resolve().parent,Path('submission'),Path('/snapshot/files/submission'))
 ROOT=next((p for p in roots if all((p/f).is_file() for f in (P,W,V,S,R,'bound.json','submission.json'))),None)
 if ROOT:TEXT={f:(ROOT/f).read_text() for f in (P,W,V,S,R)}
 res=run();out=json.dumps(res,sort_keys=True);assert len(out.encode())<65536
 print(out);raise SystemExit(0 if res['all_passed'] else 1)
