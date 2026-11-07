from fractions import Fraction as F
from math import ceil,floor,isqrt,comb,prod
from pathlib import Path
import json,sys,re,hashlib
EPS=F(1,2**20);NU=F(1,2**80); NE=F(1,10**18); EB=F(1,10**18);ER=F(1,2**80);N0=2**60;AB=3
SL={'a': '25892596493064607817/1000000000000000000', 's': '24504334046532303911/1000000000000000000', 's0': '25889184193064607817/1000000000000000000', 'hc': '299996959556134900129/1000000000000000000', 'd': '1061923791119220918633/500000000000000000', 'k': '613584029589046517201433/250000000000000000'}
VN={'I': '14620761,14620737,14620576,14619745,14613953,14577362,14344183', 'F': '14625759,14625733,14625585,14624666,14618871,14582204,14348622', 'G': '28567875,27012483', 'J': '40212127,36192167', 'Ext0': '423643930,423626488,423529329,422999506,420107819,405083279,415401918', 'Ext2': '157471115,154226162,144036645', 'L0': '32736084,32735887,32736659,32736084,32744880,32742017,32734657,32714367,32656862,32492665,32025143,30714605,27189496', 'Lr': '46193304,46183418,46163938,46122171,46028828,45806071,45256935,43909426,40730226', 'E': '12272017,12309611,12339099,12351961,12323100,12145019'}
VE={k:','.join(str(F(int(x),10**7)) for x in s.split(',')) for k,s in VN.items()}
QS=F('29970195661516549467101813/10000000000000000000000000');QC=20;DX=80 # linear-and-bilinear-primitives.md:3; proof.md:745
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
JS=(5,5,5,5,5,5,6);JY=(5,5,4);BM=prod(JS);BY=prod(JY)
def envm(vec,js=JS):
 assert len(vec)==len(js)
 p=F(1);g=F(0)
 for j,r in zip(js,vec):
  lam=j*r/(r-j+1);g+=p*(lam+F(j-1,j)*(r+1));p*=lam
 h=len(vec);B=prod(js);j=js[-1];kap=(vec[-1]+1)/(vec[-1]-j+1);dr=max(0,h+1+F(1,j)+(j-18)*kap)
 return dict(lam=p/B,g=g+2*p,c=kap/B,pi=sum(r+1 for r in vec)+3*kap+(vec[-1]+1)/j+dr*EPS)

def et(e):return ','.join(k+'='+str(e[k]) for k in ('c','pi','lam','g'))
def walkextra(s):return ((2**16*25*(s+4)**40).bit_length()-1)//7+1
def slot(x):return F(ceil(x*10**18)+2,10**18)
def run():
 v={k:list(map(F,s.split(','))) for k,s in VE.items()};i={k:F(s) for k,s in SL.items()}
 ev={k:(envm(a,JS if k=='Ext0' else JY) if k in ('Ext0','Ext2') else env(a,True)) for k,a in v.items()}
 I,FF,G,J=[ev[k] for k in ('I','F','G','J')];x,y=ev['Ext0'],ev['Ext2']
 DR=max(0,len(v['Ext0'])+1+F(1,JS[-1])+(JS[-1]-18)*x['c']*BM);PB=x['pi']-DR*EPS
 eta=F(91153006798943039772845595,41819161952109252411436045792);theta=eta/6;gamma=F(18230601359788607954569119,193428131138340667952988160);U=2**40
 alpha=F(433)/(433+2*QS)+NE
 df=i['d']*alpha+NE;kf=i['k']*alpha+NE
 M0=y['c']*df+y['pi']*alpha+NE;K2=y['lam']*df+y['g']*alpha+NE
 rr=(M0+PB*alpha+F(149,18)*alpha+(53+DR)*NE)/(theta-6*x['c'])+NE
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
 pt(P,f'T=ceil((6p+2)/Q_S+eps)={T}')
 pt(P,f'T=ceil((6p+2)/Q_S+eps)={T}')
 pt(P,f'T=ceil((6p+2)/Q_S+eps)={T}')
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
  for ridx,rate in enumerate(vec):
   ck('PT.rate',rate,JS[ridx]-1 if nm=='Ext0' else JY[ridx]-1 if nm=='Ext2' else 1,'>')
   sh=ceil(7*(rate+1));nn=ceil(F(3*(sh-14),2))-2
   bc=ceil(F(2*(nn+2),sh-14));mm=ceil(F(bc*sh,2));kk=mm+nn+1
   ck('PT.boundary',kk,bc*(sh-7),'<')
  f=W if nm in ('I','F','G','J') else P
  if nm in ('I','F','G','J'):pt(W,'Section7g recurrences at vMI,vMF,vMG,vMJ')


  if nm not in ('Ext0','Ext2'):
   out=dict(I=i['hc']+EPS,F=i['hc']+EPS,G=i['a']+EPS,J=i['s']+EPS,L0=sc,Lr=RATE,E=DX)[nm]
   batch('PXG.',[('min',vec[-1],F(16,15),'>='),('max',vec[-1],64,'<='),('output',out,256*2**h,'<='),('error',3*2**h-2,2**(h+4),'<')])
 pt(V,'## 7h. Mixed-arity exact primitive PXTM');pt(W,'Ext0 is PXTM with7 mixed-arity levels');pt(W,'Ext2=PXTM(JY,vNB2)');pt(P,'bases 1067<')
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
 ck('AE.denominator',theta,6*xu['c'],'>');ru=(M0u+(xu['pi']-DR*EPS)*alpha+F(149,18)*alpha+(53+DR)*NE)/(theta-6*xu['c'])+NE
 d0u=6*xu['c']*ru+(xu['pi']-DR*EPS)*alpha+(DR+1)*NE;d1u=2*(M0u+d0u+3*alpha)+NE;dcu=2*d0u+2*d1u+df+K2u+alpha
 aec=[ceil(z) for z in (M0u,K2u,ru,d0u,d1u,dcu,6*xu['lam']*ru+xu['g']*alpha,kf+12*ru+DX+2*alpha)]
 ck('AE.seed',dcu,sc-NE,'<=');ck('AE.source',max(6*xu['lam']*ru+xu['g']*alpha,kf+12*ru+DX+2*alpha),kc-NE,'<=');
 Iu,Fu,Gu,Ju=[eu[k] for k in ('I','F','G','J')]
 for nm,ky,val in IR(Iu,Fu,Gu,Ju):ck('AE.inner.'+nm,i[ky],val,'>=')
 Eu,Ru,Zu=eu['E'],eu['Lr'],eu['L0']
 for nm,lim,val in [('dyE',DY,Eu['c']*DX+Eu['pi']),('dxLr',DX,Ru['c']*RATE+Ru['pi']),('s0L0',s0,Zu['c']*sc+Zu['pi']),('d3d',do,tr+Eu['lam']*DX+Eu['g']+t*DX+alpha),('dKE',do,Eu['lam']*DX+Eu['g']),('k3b',ko,tl+Ru['lam']*RATE+Ru['g']+t*RATE+alpha),('kK0',ko,Zu['lam']*sc+Zu['g'])]:ck('AE.outer.'+nm,lim,val+(NE if nm in ('s0L0','d3d','dKE','k3b','kK0') else 0),'>=')
 for val in (k00,kother):ck('NB.source',val,kc,'<=')
 ck('NB.seed',dc,sc,'<=');ck('NB.source-dimension',kc,do,'<=');ck('NB.seed-dimension',sc,do,'<=')
 batch('NB.',[('denominator',theta,6*x['c'],'>'),('p0',PB,16,'>='),('height',max(lg(BM),lg(BY)),38,'<='),('cap',DX,128,'<='),('advice',F(11,2)+F(25,3)+F(53,U),28,'<='),('constant',ceil(2*(rr+2)*M0)+1,2**240,'<'),('Lm',RATE*DY+1,2**33,'<'),('heights',max(len(v[k]) for k in ('I','F','G','J'))+13,21,'<='),('outer-heights',max(len(v[k]) for k in ('L0','Lr','E'))+12,25,'<='),('NB-heights',max(lg(BM),lg(BY))+6,27,'<='),('prime',6*2**80/eta,2**90,'>'),('capacity',F(7,6)*6,7,'=='),('floor',53-F(34,3)-1,F(122,3),'=='),('advice40',F(122,3),40,'>'),('error',14,16,'<'),('final-error',15,16,'<')])
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
 batch('local.',[('alpha',alpha,1,'<'),('output',F(DX,128)+NE,alpha,'<='),('source',do,2**30,'<'),('Omega',co*Q*(1+F(1,2**80)),1024,'<'),('unit-base',490628,NE*(433+2*Q)*2**80,'<'),('unit-slope',5770,490628,'<='),('hash-base',804,NE*2**80,'<'),('hash-slope',4,804,'<='),('round',1,NE*2**80,'<='),('advice-charge',F(11,2)+F(25,3)*F(1,3),F(149,18))])
 pt(P,'5770*80+29028=490628');pt(V,'16/15<=Rh<=64')
 for ss in sorted({s for k in range(80,121) for s in (2**k,2**(k+1)-1)}|{2**80+j for j in range(8)}):
  k=ss.bit_length()-1;s=8*ss;ellc=lg(s+4);rf=ss+walkextra(s);ab=rf+AB;u1=2**(110+lg(ss));rho=ab+2*lg(u1)+2
  walk=s+40*rf;tag=3*walk
  batch('AFexec.',[('rho',rho,rf+225+2*lg(ss),'=='),('rho-dyadic',rho-ss,10*k+267,'<='),('rho-bound',rho,(1+EB)*ss,'<='),('log',ellc,k+5,'<='),('walk',2**(7*(rf-ss)),2**16*25*(s+4)**40,'>'),('parity',walk%2,0,'==')])
  for L in (ceil(Q*ss),ceil(Q*(ss+1))-1):
   ee=3*tag+rf+28+lg(C);om=ee+lg(26*t)+2*L+42;d=floor(do*om);kk=ceil(aa*d+bb*rho)
   eeN=ee+lg(26*t)+4;np=d;adv=tag+2+lg(C)
   wn=max(2**80,eeN+44+2*lg(4*np));wf=max(2**80,eeN+2*lg(4*np)+32,ceil(F(DX*om,128)))
   batch('localexec.',[('Omega',om,1024*ss,'<='),('nlog',lg(4*np),42+lg(ss),'<='),('wN',wn,alpha*om,'<='),('wF',wf,alpha*om,'<='),('advice',3*adv,eeN,'<='),('hash',480+4*lg(om),NE*om,'<=')])
   for shift in (2,3):
    rho0=ab+L+shift+lg(kk)+2;b0=ceil(a0*rho0);u0=b0-rho0;b=ceil(gb*rho);uu=b-rho
    v0=b0*ceil(F(kk,u0));v1=b*ceil(F(kk,uu));cc=v1-kk;v2=b*ceil(F(cc,uu));D=b0+2*b
    ee=3*tag+rf+28+lg(C);om=ee+lg(26*t)+2*L+42
    batch('AFexec.',[('floor',floor(L/Q),ss,'=='),('rho0',rho0-(Q+1)*ss,9*k+Q+97,'<='),('rho0-bound',rho0,(Q+1+EB)*ss,'<='),('N0',rho0,N0,'>='),('u0',u0,1,'>='),('ratio0',F(b0,u0),N0+1,'<='),('u',uu,1,'>='),('ratio',F(b,uu),xi,'<='),('Kbound',kk,kb*ss,'<='),('Klog',lg(kk),48+lg(ss),'<='),('v0',v0,u1,'<='),('v1',v1,u1,'<='),('c',cc,1,'>='),('out',v1-v2,d,'>='),('rank-bound',v1-v2,gg*kk-xi*b,'>='),('seed',D,s,'<='),('source',kk,bh*L,'<='),('tag',tag,2*L,'>='),('runtime',tag,2*p*L/Q,'<='),('unit',om,co*L,'<='),('tail',8+2*ee,T*L,'<='),('combined',4*C*d+ko*om+8+2*ee,H*L,'<=')])
 batch('',[('SP.gamma',F(377,200)-20*gamma,F(377,200)*NU,'=='),('SP.guard',37,F(377,200)*NU*2**90,'<'),
 ('SP.subset',(F(3771,40000)-gamma)*2**90,3,'>'),('SP.dimension',18*gamma+F(33,2**90),F(2),'<'),
 ('SP.power',2+18*gamma+F(33,2**90),4,'<'),('SP.eta',(F(1081,25)+7*NU)*eta,gamma,'=='),('SP.cap',eta,F(1,458),'<='),('SP.theta',6*theta,eta,'=='),('CZ.exponent',6*F(7,6)+6*(F(76,25)+NU)+18,F(1081,25)+6*NU,'=='),('CZ.hered',(F(121,20)-(F(76,25)+NU)-3)*64/NU,6,'>'),('CZ.zero',1-9*eta,0,'>'),('CZ.T',2+F(13,6)*eta,4,'<='),('CZ.smoothing-exponent',F(7,6),1,'>='),
 ('SP.cardinal',5+3+1,9,'=='),('SP.sixteen',9,16,'<'),('SP.endpoint',F(15+4)*3/8,F(57,8),'=='),('SP.endpoint-loss',F(57,8),8,'<'),
 ('CZ.units',1-F(493,27480),F(26987,27480),'=='),('CZ.plane',F(76229,40000)-F(433,27480)-1,F(24456323,27480000),'=='),
 ('CZ.projection',F(76229,40000)-F(563,27480)-F(377,200),F(6523,27480000),'=='),('CZ.energy',2**63+2**57,2**64,'<'),
 ('CZ.smoothing',30+3+2+10,45,'<='),('CZ.hit',35,36,'<'),('CZ.root16',16**6,2**(7*64),'<'),
 ('prime.anchor',F(32,5),F(99,16),'>'),('prime.square',631,26**2,'<'),
 ('native.union',32*2**16*2*2**6,2**28,'=='),('native.collision',F(2,2048)+F(1,4),F(257,1024),'=='),
 ('native.r',F(257,1020)+F(1,15*2**256),F(13,50),'<'),('native.mono',6-8*F(13,50),0,'>'),
 ('native.quartic',6*F(13,50)-4*F(13,50)**2-2,F(-444,625),'=='),('native.margin',F(-444,625)+F(1,1000)+F(2048,3825),F(-1,8),'<'),
 ('native.g',F(4,2**4095),F(1,1000),'<'),('native.delete',F(1,2**27),F(1,10000**2),'<'),('native.perturb',F(432,10000),F(1,8),'<'),
 ('native.moment',3+F(3,1024),4,'<='),('native.success',1-2*F(2,5),F(1,5),'=='),
 ('F.entropy',2*13+12,38,'=='),('F.seed3',3*13,39,'=='),('F.seed',5*13,65,'=='),('F.labels',1+2*4,9,'=='),('F.error',76,128,'<'),('F.round',64,76,'<='),('native.H128',H*128,2**128,'<')])
 ck('SP.constants',(8**4,2**15,2**17*2**20,2**15*2**18),(4096,32768,2**37,2**33))
 ck('native.ratio',F(257,1024)*F(256,255),F(257,1020))
 ck('native.bias',F(1,16)+F(1,2**23),F(1,15),'<')
 ck('native.marginal',3+F(24*256,15*255),16,'<');ck('native.tail',F(432,10000),F(1,16),'<')
 ck('native.weight-count',sum(2**j for j in range(1,5)),32,'<=');ck('native.moment',F(8*256,15*255),F(11,20),'<')
 ck('native.deletion',F(1,2**28)*F(256,255),F(1,2**27),'<')
 for kp_p in (5,7,11):
  kp_d=next(x for x in range(1,kp_p) if all(t*t%kp_p!=x for t in range(kp_p)))
  kp_mul=lambda x,y:((x[0]*y[0]+kp_d*x[1]*y[1])%kp_p,(x[0]*y[1]+x[1]*y[0])%kp_p)
  kp_conj=lambda x:(x[0],(-x[1])%kp_p)
  kp_vals=[(a,b) for a in range(kp_p) for b in range(kp_p)]
  ck('quadratic.basis',kp_mul((0,1),(0,1)),(kp_d,0))
  ck('quadratic.roots',[x for x in kp_vals if kp_mul(x,x)==(kp_d,0)],[(0,1),(0,kp_p-1)])
  for kp_x in kp_vals:
   if kp_x!=(0,0):
    kp_n=(kp_x[0]**2-kp_d*kp_x[1]**2)%kp_p;ck('quadratic.norm',kp_n,0,'>')
    kp_ni=pow(kp_n,kp_p-2,kp_p);kp_y=(kp_x[0]*kp_ni%kp_p,-kp_x[1]*kp_ni%kp_p)
    ck('quadratic.inverse',kp_mul(kp_x,kp_y),(1,0))
   for kp_y in kp_vals:ck('quadratic.conjugation',kp_mul(kp_conj(kp_x),kp_conj(kp_y)),kp_conj(kp_mul(kp_x,kp_y)))
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
  vec=v[nm];hh=len(vec);js=JS if nm=='Ext0' else JY if nm=='Ext2' else (2,)*hh;BB=prod(js);aa=ee+lg(BB)+4;rho=aa+2*lg(4*nn)+2
  ww=max(unit,2**20,rho);ell=ceil(F(zz,BB))
  slots=[ceil((rr+1)*rho) for rr in vec];kk=ell+2*aa
  ck('PT.leaf-threshold',kk,F(zz,BB)+2*ww-11,'<=')
  ck('PT.expansion-guard',4*(4*nn)**2,2**(rho-aa),'<=')
  ck('PT.a-unit',aa,ww-6,'<=')
  for ix in range(hh-1,-1,-1):
   nxt=kk;jj=js[ix];bb=ceil(F(jj*(nxt+aa+jj-1),slots[ix]-jj*rho));mm=ceil(F(bb*slots[ix],jj));kk=(jj-1)*mm+nxt+aa
   ck('PX.capacity',kk,bb*(slots[ix]-rho),'<');ck('PX.count',bb,kk,'<=')
   ck('PX.split',kk,(jj-1)*mm+nxt+aa,'==');ck('PX.child',mm,nxt,'>')
   if ix==hh-1:
    leaf=mm;kap=(vec[-1]+1)/(vec[-1]-jj+1)
    ck('PT.terminal',leaf,kap*(ceil(F(zz,BB))+3*aa+jj-1)+F(slots[-1],jj)+1,'<')
  jj=js[-1];ag=nm not in ('Ext0','Ext2')
  if ag: tail=2*(lg(leaf*ceil(F(zz,BB))+1)+ceil(F(zz,BB))+2*aa)
  else:tail=leaf
  dd=sum(slots)+tail;en=envm(vec,JS if nm=='Ext0' else JY) if nm in ('Ext0','Ext2') else env(vec,ag)
  kap=(vec[-1]+1)/(vec[-1]-jj+1)
  ck('PT.terminal-unit',kap*(ell+3*aa+jj-1)+F(slots[-1],jj)+1,
   kap*F(zz,BB)+(3*kap+(vec[-1]+1)/jj)*ww+1+F(1,jj)+(jj-18)*kap,'<=')
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
 ee0=3*tag+rf+28+lg(C)+lg(26*t);om=ee0+2*L+42
 nn=floor(do*om);zz=DX*om;adv=tag+2+lg(C);ee=ee0+4
 wn=max(2**80,ee+44+2*lg(4*nn));wf=max(2**80,ee+2*lg(4*nn)+32,ceil(F(zz,128)))
 m1,kf1=ceil(i['d']*wf),ceil(i['k']*wf);m0,k2=px('Ext2',nn,m1,ee)
 sig=ee+4;rr1=ceil((m0+PB*wn+DR+F(11,2)*ee+F(25,3)*adv+53)/(theta-6*x['c']))
 ss=floor((eta*rr1-15*sig-8)/6);d01,k01=px('Ext0',nn,6*rr1,ee)
 el=ceil(F(d01+ee-1,2))+ee;d11=2*(lg(2*(rr1+adv+1)*m0+1)+m0+2*el);d2=d01+d11+m1+k2+ee
 batch('NBexec.',[('I1',ss,m0+2*el,'>='),('advice',ss,adv+40,'>='),('density',F(adv+sig,rr1+adv),eta/50,'<='),('source',rr1-sig,(1-eta/50)*(rr1+adv+1)+1,'>='),('capacity',F(7,6)*eta*(rr1+adv),7*ss+15*sig,'>='),('density-error',eta*(rr1+adv),6*ss+15*sig+8,'>='),('prime',rr1+adv,2**90,'>'),('dC',d01+d11+d2,floor(sc*om),'<='),('kC',max(k01,kf1+12*rr1+zz+2*ee),floor(kc*om),'<=')])
 ck('refine.log-slope',361*F(40,7),F(14440,7),'==');ck('refine.Econst',361*F(14576685580,1210361901)+35,F(5304546160915,1210361901),'==')
 ck('refine.Oconst',361*F(14576685580,1210361901)+91,F(5372326427371,1210361901),'==');ck('refine.tailconst',2*F(5304546160915,1210361901)+8,F(10618775217038,1210361901),'==')
 UO=4*C*do+ko;LA=433/Q+2;A=UO*LA+866/Q;B=max(F(14440,7)*(UO+2),do*F(14440,7)/gg+F(54,7)*Bbase);CC=ceil(max(F(5372326427371,1210361901)*UO+F(10618775217038,1210361901),do*F(5372326427371,1210361901)/gg+F(570991104275,2420723802)*Bbase+1))
 ck('refine.Qpower',Q**64924,2**102809,'>');ck('refine.logmargin',(8/Q+F(4,2**80))**64924,2**91963,'<');ck('refine.pow25',25**10653,2**49471,'<');
 batch('refine.',[('source',A,do*LA/gg+Bbase/Q,'>='),('log',B,do*F(14440,7)/gg+F(54,7)*Bbase,'>='),('constant',CC,do*F(5372326427371,1210361901)/gg+F(570991104275,2420723802)*Bbase+1,'>='),('leading',A,H-1,'<='),('guard-base',B*81+CC+1,2**80,'<'),('guard-slope',B,B*81+CC+1,'<=')])
 ck('result.H',H,1119447578339);ck('result.K',K,4477790313356)
 ck('result.raw',raw,1119447578339);pt(P,f'H_raw=4C*delta_O+kappa_O+T={raw}');ck('result.base',bh,1119447578339)
 ck('result.improvement',4*A,F('375352142468710314199053223984/1498508685085561'),'<');ck('result.loss',[4*A,4*B,4*(CC+1)],[F('4193757868301077995136034350912430016438245474423786143/936568614422392170846931656250000000000000'),F('1458421632047435346858531996366559157395802213565226688486875351650704021577175687432506561823908154429171261351207/23126996603663247089766053616958354107773090821825892022574324329273248504494221250431250000000000000'),135688130693588])
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
HASHES={'proof.md': 'cae06639504e944cfb9421ea5a41ae395cf58240a4df9f5920a41823a3dfc339', 'bound.json': 'cb483ae34813c3a3f3737b32ce94653b49281aa37444d795d46d7aab69e63f02', 'sampler-and-disperser.md': 'eb80a44363137f72bf1511901a3a20429111bfa2055007774283e70237efea77', 'restricted-nonmalleability-and-breaking.md': 'a0b6a75f052c7d8f261e420ae69a4be6280e3b25d59c979e49172ae8b778af23', 'linear-and-bilinear-primitives.md': '41fd9fed4253de53b25768249ad293f6aa152c9658f87a7bc5c0cf968a710312', 'affine-and-standard-wrappers.md': '0575aab061bfe2ea4b32d6becd2d34b9efea1353c4191560cd5dc4e3dff79401', 'submission.json': '8654309ef77b0a38ad76fcf1b70bba7a7ae324930c361857113060b37e9ae45f'}
if __name__=='__main__':
 if json.load(sys.stdin)!={}:raise ValueError('Expected {}')
 roots=(Path(__file__).resolve().parent,Path('submission'),Path('/snapshot/files/submission'))
 ROOT=next((p for p in roots if all((p/f).is_file() for f in (P,W,V,S,R,'bound.json','submission.json'))),None)
 if ROOT:TEXT={f:(ROOT/f).read_text() for f in (P,W,V,S,R)}
 res=run();out=json.dumps(res,sort_keys=True);assert len(out.encode())<65536
 print(out);raise SystemExit(0 if res['all_passed'] else 1)
