from fractions import Fraction as F
from math import ceil,floor,isqrt,comb,prod
from pathlib import Path
import json,sys,re,hashlib
EPS=F(1,2**20);NU=F(1,2**80); NE=F(1,10**18); EB=F(1,10**18);ER=F(1,2**80);N0=2**60;AB=3
SL={'a': '8865705637228501403/500000000000000000', 's': '8188722818614250703/500000000000000000', 's0': '8865099637228501403/500000000000000000', 'hc': '114352533880452100939/1000000000000000000', 'd': '420043978928429477387/500000000000000000', 'k': '542651172076908700115091/1000000000000000000'}
VN={'I': '2183488000,2183441000,3186811000', 'F': '2183935000,2183880000,3187137000', 'G': '2511740000', 'J': '4542688000', 'Ext0': '32805689000,32799192000,32762327000,32554520000,31436509000,26611077000', 'Ext2': '12326277000,14868047000,12694543000', 'L0': '2282932000,2282928000,2282916000,2282870000,2282709000,2282134000,2280087000,2272808000,2246950000,2155452000', 'Lr': '6601218373,6600100000,6597355000,6590621000,6574211000,6534146000,6437152000,6207296000,5691885000', 'E': '1187570000,1187581000,1187542000,1187225000,1182928000'}
VE={k:','.join(str(F(int(x),10**9)) for x in s.split(',')) for k,s in VN.items()}
QS=F('10696926570874890811620649/5000000000000000000000000');QC=17;DX=81 # linear-and-bilinear-primitives.md:3; proof.md:747
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
 if TEXT:tc('text.'+tok[:70],' '.join(tok.split()) in ' '.join(TEXT[f].split()))
def lg(x):
 x=F(x);k=max(0,x.numerator.bit_length()-x.denominator.bit_length());return k+(2**k<x)
def env(vec,ag=False):
 p=F(1);g=F(0)
 for r in vec:
  lam=2*r/(r-1);g+=p*(lam+(r+1)/2);p*=lam
 h=len(vec);r=vec[-1];kap=(r+1)/(r-1)
 return dict(lam=p/2**h,g=g+2*p,c=F(1,2**h),pi=4+sum(r+1 for r in vec))
JS=(5,5,5,5,5,5);JY=(3, 4, 4);JI=(3,3,4);BM=prod(JS);BY=prod(JY)
def envm(vec,js=JS):
 assert len(vec)==len(js)
 p=F(1);g=F(0)
 for j,r in zip(js,vec):
  lam=j*r/(r-j+1);g+=p*(lam+F(j-1,j)*(r+1)+max(0,(j-7)*lam+F(j-1,j))*EPS);p*=lam
 h=len(vec);B=prod(js);j=js[-1];kap=(vec[-1]+1)/(vec[-1]-j+1);dr=max(0,h+1+F(1,j)+(j-18)*kap)
 return dict(lam=p/B,g=g+2*p,c=F(1,B),pi=4+sum(r+1 for r in vec))
def ar(nm,h):return JS if nm=='Ext0' else JY if nm=='Ext2' else JI if nm in ('I','F') else (2,)*h
def et(e):return ','.join(k+'='+str(e[k]) for k in ('c','pi','lam','g'))
def walkextra(s):return 0
def slot(x):return F(ceil(x*10**18)+2,10**18)
def canonical_leaf(N,m,a):
 D=1<<(2*a+1);B=D+1;d=m+4*a+lg(4*(N+1));T=1<<d
 xs=list(range(1,1<<N));weights=[1]*len(xs);table=[]
 def potential(rows):
  total=0
  for x,w in zip(xs,weights):
   free=0;bad=False
   for row in rows:
    fixed=0;undetermined=False
    for j,v in enumerate(row):
     if x>>j&1:
      if v is None:undetermined=True
      else:fixed^=v
    if undetermined:free+=1
    elif fixed:bad=True;break
   total+=w*((D<<m)+(0 if bad else 1<<(m-free)))
  return total
 for s in range(T):
  rows=[[None]*N for _ in range(m)]
  for i in range(m):
   for j in range(N):
    rows[i][j]=0;p0=potential(rows);rows[i][j]=1;p1=potential(rows);rows[i][j]=int(p1<p0)
  masks=[sum(v<<j for j,v in enumerate(row)) for row in rows];table.append(masks)
  for t,x in enumerate(xs):weights[t]*=B if all((r&x).bit_count()%2==0 for r in masks) else D
 counts=[sum(all((r&x).bit_count()%2==0 for r in masks) for masks in table) for x in xs]
 assert max(counts)<T*F(1,1<<m)*(1+F(1,1<<(2*a)))
 return table,counts

def run():
 table,counts=canonical_leaf(5,1,1)
 ck('CLtable.size',len(table),1024)
 for nzero in counts:ck('CLtable.collision',nzero,640,'<')
 for source in (list(range(8)),[0,1,2,4,8,16,3,31],[1,3,5,7,9,11,13,15]):
  excess=sum(abs(2*sum((row[0]&xx).bit_count()%2 for xx in source)-len(source)) for row in table)
  ck('CLtable.strong',F(excess,2*len(source)*len(table)),F(1,2),'<')
 v={k:list(map(F,s.split(','))) for k,s in VE.items()};i={k:F(s) for k,s in SL.items()}
 ev={k:(envm(a,ar(k,len(a))) if k in ('Ext0','Ext2','I','F') else env(a,True)) for k,a in v.items()}
 I,FF,G,J=[ev[k] for k in ('I','F','G','J')];x,y=ev['Ext0'],ev['Ext2']
 SD=F(999,1000)*(sum(rr+1 for rr in v['Ext2'])+i['d']/BY+4);DR=F(0);PB=x['pi']-DR*EPS
 eta=F(1172658045026190299464990235,142840025753656702822820907112);theta=eta/4;h_R=F(1456,485)+2*NU;gamma=F(24111854012685236344524388261014376388159010267723399,172682995208014168848897321782865682988787662388723712);delta_R=gamma+NU;q_R=14*gamma-eta+NU;U=2**40
 eb=F(2599,7)-9*F(1,2**100);alpha=eb/(eb+2*QS)+NE
 df=i['d']*alpha+NE;kf=i['k']*alpha+NE
 M0=y['c']*df+y['pi']*alpha+NE;K2=y['lam']*df+y['g']*alpha+NE
 rr=(M0+PB*alpha+F(23,4)*alpha+(56+DR)*NE)/(theta-(6-eta/2)*x['c'])+NE
 d0=(6-eta/2)*x['c']*rr+PB*alpha+(DR+2)*NE;d1=2*(M0+d0+3*alpha)+NE
 k00=(6-eta/2)*x['lam']*rr+x['g']*alpha+NE;kother=kf+(12-eta)*rr+DX+2*alpha+NE
 dc=2*d0+2*d1+df+K2+alpha;sc=slot(dc);kc=slot(max(k00,kother))
 q=QC;C=4*(q+1);t=4*C-1;j=lg(t);Q=QS;beta=F(2*(q+8),7);p=144
 EE,RR,ZZ=[ev[k] for k in ('E','Lr','L0')]
 DY=ceil(EE['c']*DX+EE['pi']);RATE=(t+2)*DY+3
 s0=slot(ZZ['c']*sc+ZZ['pi']);tr=(t+1)*(s0+sc+2*j*(DY+DX));tl=(t+1)*(sc+2*j*DY)
 rows=dict(dy_E=EE['c']*DX+EE['pi'],dx_Lr=RR['c']*RATE+RR['pi'],dx_Lm=2*DY+4+NE,s0_L0=ZZ['c']*sc+ZZ['pi'],r_3a=(t+2)*DY+3,r_dy=DY,r_Km=DY+2,
 d_3c=(t+1)*(s0+sc)+kc+alpha,d_3d=tr+EE['lam']*DX+EE['g']+t*DX+alpha,d_s0=s0,d_KE=EE['lam']*DX+EE['g'],d_kC=kc,
 k_3b=tl+RR['lam']*RATE+RR['g']+t*RATE+alpha,k_K0=ZZ['lam']*sc+ZZ['g'])
 do=slot(max(val for key,val in rows.items() if key.startswith('d_')));ko=slot(max(val for key,val in rows.items() if key.startswith('k_')))
 co=(3*p+1)/Q+2+F(1,10**15);delta=ceil(do*co);kappa=ceil(ko*co);T=ceil((6*p+2)/Q+EPS)
 pt(P,f'T=ceil((6p+2)/Q_S+eps)={T}')
 a0=F(N0+1,N0);gb=(beta-a0*(Q+1+EB)-3*ER)/(2*(1+EB));xi=gb/(gb-1);gg=1-(xi-1)**2;aa=1/gg;bb=ceil((xi*gb+xi/U)/gg)
 Bbase=bb;raw=4*C*delta+kappa+T;bh=ceil(aa*delta+bb*(1+EB)/Q+ER);H=max(raw,bh);K=4*H
 kb=aa*delta*Q*(1+F(1,U))+bb*(1+EB)+ER;b0c=a0*(Q+1+EB)+ER
 IR=lambda I,FF,G,J:[('aJ','a',J['c']*(i['s']+NE)+J['pi']),('aF','a',FF['c']*(i['hc']+NE)+FF['pi']),('sG','s',G['c']*(i['a']+NE)+G['pi']),('s0I','s0',I['c']*(i['hc']+NE)+I['pi']),('hJ','hc',(5+J['lam'])*(i['s']+NE)+J['g']+1),('hcap','hc',F(100)),('h>=s','hc',i['s']+NE),('dG','d',9*(i['s0']+NE)+(36+G['lam'])*(i['a']+NE)+G['g']+1),('d>=s0','d',i['s0']+NE),('kF','k',(26+FF['lam'])*(i['hc']+NE)+18*(i['s']+NE)+FF['g']+1),('kI','k',I['lam']*(i['hc']+NE)+I['g']),('k>=d','k',i['d']+NE)];ir=IR(I,FF,G,J)
 for nm,ky,val in ir:
  ck('inner.'+nm,i[ky],val,'>=')
 for ky,val in i.items():pt(W,f'{ky}=ceil(({val})w).')
 for ky,val in rows.items():
  lim=dict(dy=DY,dx=DX,r=RATE,s0=s0,d=do,k=ko)[ky.split('_')[0]]
  ck('outer.'+ky,lim,(F(ceil(val*10**18),10**18)+NE if ky.split('_')[0] in ('s0','d','k') else val),'>=');pt(P,f'{ky}: '+('U18(rhs)+eps_N' if ky.split('_')[0] in ('s0','d','k') else 'rhs')+f' <= {lim}.')
 den=theta-(6-eta/2)*x['c'];runit=(i['d']/BY+y['pi']+PB+6)/den
 batch('CLguard.',[('Om',co*Q*(1+F(1,2**80)),440,'<='),('wF',DX*440,128*370,'<='),('Fguard',(100*370-DX*440)*2**160,1200,'>'),('Ext2cap',i['d']+NU,256*BY,'<='),('Ext0cap',(6-eta/2)*runit+NU,256*BM,'<='),('rround',56+den+F(1,BY),F(2**80,4),'<=')])
 for nm,vec in v.items():
  en=ev[nm];h=len(vec)
  for ridx,rate in enumerate(vec):
   ck('PT.rate',rate,ar(nm,len(vec))[ridx]-1,'>')
   sh=ceil(7*(rate+1));nn=ceil(F(3*(sh-14),2))-2
   bc=ceil(F(2*(nn+2),sh-14));mm=ceil(F(bc*sh,2));kk=mm+nn+1
   ck('PT.boundary',kk,bc*(sh-7),'<')
  f=W if nm in ('I','F','G','J') else P
  if nm in ('I','F','G','J'):pt(W,'Sections7g/7h recurrences at vMI,vMF,vMG,vMJ')
  kap=(vec[-1]+1)/(vec[-1]-ar(nm,len(vec))[-1]+1)
  batch('CLguard.',[('kap',kap,31,'<='),('last',vec[-1],64,'<='),('height',h,13,'<=')])
  if nm not in ('Ext0','Ext2'):
   out=dict(I=i['hc']+EPS,F=i['hc']+EPS,G=i['a']+EPS,J=i['s']+EPS,L0=sc,Lr=RATE,E=DX)[nm]
   ck('CLguard.output',out,256*prod(ar(nm,h)),'<=')
 pt(V,'## 7h. Mixed-arity exact primitive PXTM');pt(W,'Ext0 is PXTM with6 mixed-arity levels');pt(W,'Ext2=PXTM(JY,vNB2)');pt(P,'bases 1067<')
 pt(W,'p0='+str(PB));pt(W,'D0='+str(DR))
 pt(W,'c0='+str(x['c']));pt(W,'pi0='+str(x['pi']))
 pt(W,'bc0='+str((6-eta/2)*x['c'])+'<theta='+str(theta));pt(P,f'delta_O=ceil({do}c_O)={delta}');pt(P,f'kappa_O=ceil({ko}c_O)={kappa}');pt(P,'D_F=alpha*d_inner+eps_N='+str(df));pt(P,'K_F=alpha*k_inner+eps_N='+str(kf));pt(P,str(co))
 for val in (gb,xi,gg,kb):pt(V,str(val))
 cl=[ceil(z) for z in (M0,K2,rr,d0,d1,dc)];pt(P,'ceil(M0)=%d, ceil(K2bar)=%d,\n ceil(Rnew)=%d, ceil(d0bar)=%d,\n ceil(d1bar)=%d, ceil(DNB)=%d.'%tuple(cl))
 pt(P,'KC='+str(kc));pt(P,'SC=Slot(DNB)='+str(sc))
 eu={k:{c:F(ceil(e[c]*10**40),10**40) for c in e} for k,e in ev.items()}
 for c in ('c','pi','lam','g'):
  for k in ev:ck('AE.enclosure',0<=eu[k][c]-ev[k][c]<F(1,10**40),True)
 xu,yu=eu['Ext0'],eu['Ext2'];M0u=yu['c']*df+yu['pi']*alpha+NE;K2u=yu['lam']*df+yu['g']*alpha+NE
 ck('AE.denominator',theta,(6-eta/2)*xu['c'],'>');ru=(M0u+(xu['pi']-DR*EPS)*alpha+F(23,4)*alpha+(56+DR)*NE)/(theta-(6-eta/2)*xu['c'])+NE
 d0u=(6-eta/2)*xu['c']*ru+(xu['pi']-DR*EPS)*alpha+(DR+2)*NE;d1u=2*(M0u+d0u+3*alpha)+NE;dcu=2*d0u+2*d1u+df+K2u+alpha
 aec=[ceil(z) for z in (M0u,K2u,ru,d0u,d1u,dcu,(6-eta/2)*xu['lam']*ru+xu['g']*alpha+NE,kf+(12-eta)*ru+DX+2*alpha+NE)]
 ck('AE.seed',dcu,sc-NE,'<=');ck('AE.source',max((6-eta/2)*xu['lam']*ru+xu['g']*alpha+NE,kf+(12-eta)*ru+DX+2*alpha+NE),kc-NE,'<=');
 Iu,Fu,Gu,Ju=[eu[k] for k in ('I','F','G','J')]
 for nm,ky,val in IR(Iu,Fu,Gu,Ju):ck('AE.inner.'+nm,i[ky],val,'>=')
 Eu,Ru,Zu=eu['E'],eu['Lr'],eu['L0']
 for nm,lim,val in [('dyE',DY,Eu['c']*DX+Eu['pi']),('dxLr',DX,Ru['c']*RATE+Ru['pi']),('s0L0',s0,Zu['c']*sc+Zu['pi']),('d3d',do,tr+Eu['lam']*DX+Eu['g']+t*DX+alpha),('dKE',do,Eu['lam']*DX+Eu['g']),('k3b',ko,tl+Ru['lam']*RATE+Ru['g']+t*RATE+alpha),('kK0',ko,Zu['lam']*sc+Zu['g'])]:ck('AE.outer.'+nm,lim,val+(NE if nm in ('s0L0','d3d','dKE','k3b','kK0') else 0),'>=')
 for val in (k00,kother):ck('NB.source',val,kc,'<=')
 ck('NB.seed',dc,sc,'<=');ck('NB.source-dimension',kc,do,'<=');ck('NB.seed-dimension',sc,do,'<=')
 batch('NB.',[('denominator',theta,(6-eta/2)*x['c'],'>'),('p0',PB,16,'>='),('height',max(lg(BM),lg(BY)),38,'<='),('cap',DX,128,'<='),('advice',F(23,4)+F(56,U),28,'<='),('unit-p0',4*PB,64,'>='),('AG2-slope',F(3880,3),F(3883,3),'<'),('AG2-const',2**80,5820,'>'),('AG2-r-lower',eta*(PB+SD)/(theta-(6-eta/2)*x['c']),F(3883,3),'>='),('AG2-seed-round',32,F(2**80,1000),'<='),('seed-aa-base',440,F(2**160,1000),'<'),('seed-aa-slope',2,440,'<='),('constant',ceil(2*(rr+2)*M0)+1,2**240,'<'),('Lm',RATE*DY+1,2**33,'<'),('heights',max(lg(prod(ar(k,len(v[k])))) for k in ('I','F','G','J'))+13,21,'<='),('outer-heights',max(len(v[k]) for k in ('L0','Lr','E'))+12,25,'<='),('NB-heights',max(lg(BM),lg(BY))+6,27,'<='),('prime',4*16*2**80/eta,2**90,'>'),('capacity',F(5,4)*4,5,'=='),('raw-blocks',1+1+(2-eta/2)+2,6-eta/2,'=='),('raw-round',max(2,x['lam']),NE*2**80,'<='),('c0',x['c'],1,'<'),('raw-rows',4+1,5,'=='),('matching-charge',2*4+2,10,'=='),('sig-capacity',F(5,4)*11,11,'>'),('floor',56-x['c']-13-1,F(41),'>='),('advice40',F(41),40,'>'),('error',14,16,'<'),('final-error',15,16,'<')])
 for hh in (0,1,2,j):
  rt=s0+sc+DX;lf=sc;mx=0
  for _ in range(hh):
   rt+=DY;lf+=DY;mx=max(mx,rt);rt+=DX+DY;lf+=DY;mx=max(mx,rt);rt+=DX
  if hh:ck('causal.E',mx,s0+sc+2*hh*(DY+DX))
  ck('causal.full',lf+rt,s0+2*sc+DX+hh*(4*DY+2*DX))
 batch('AF.',[('ceil-round',1,ER*2**80,'<='),('gamma',gb,2,'>'),('rank',gg*aa,1,'=='),('rounding',gg*bb-xi*gb,xi/U,'>='),('B-minimal',gg*(bb-1)-xi*gb,xi/U,'<'),('seed',a0*(Q+1+EB)+2*gb*(1+EB)+3*ER,beta,'=='),('Kbits',kb,2**48,'<'),('pad0',(N0+1)*2**48+b0c,2**110,'<'),('pad1',xi*2**48+gb*(1+EB)+EB,2**110,'<'),('error',6,2**AB,'<'),('q',F(7,2)*beta-8,q,'=='),('p',p,3*Q,'>'),('H',H,2*delta,'>='),('runtime-power',ceil(2*p/Q),135,'=='),('incidence',48*U,134217728,'>')])
 for nm,c,b,budget in [('rho',267,10,EB),('rho0',Q+97,9,EB),('outer',10*5768+160,2*5768,EPS*Q),('runtime',48*q*5+6,48*q,p),('NB-log',248,8,NE),('NB-quarter',242,2,F(1,4)),('Lm',70,4,NE)]:
  kk=80 if nm in ('rho','rho0','NB-log','Lm') else 40;ck('induction.'+nm,c+kk*b,budget*2**kk,'<');ck('induction.slope',b,c+kk*b,'<=')
 ck('induction.Om-base',5768*85+96,F(1,10**15)*Q*2**80,'<');ck('induction.Om-slope',5768,5768*85+96,'<=')
 ck('PXG.legal-small',F(1,8192)+15*EPS,2,'<')
 batch('PXG.',[('log-base',64,64,'<='),('log-slope',2,64,'<='),('mh',31*(256+3)+F(65,2),8062,'<'),('dimension',8062*2**20+9,8192*2**20,'<='),('product',8192*256,2**21,'==')])
 batch('local.',[('alpha',alpha,1,'<'),('output',F(DX,128)+NE,alpha,'<='),('source',do,2**30,'<'),('Omega',co*Q*(1+F(1,2**80)),1024,'<'),('unit-base',490628,NE*(eb+2*Q)*2**80,'<'),('unit-slope',5770,490628,'<='),('hash-base',804,NE*2**80,'<'),('hash-slope',4,804,'<='),('round',1,NE*2**80,'<='),('advice-charge',F(23,4),F(23,4)),('AG1-coef',9*F(1,3),1,'>='),('AG2-coef',970*F(1,3)+970,F(3880,3),'==')])
 pt(P,f'z={DX}Om');pt(P,'5770*80+29028=490628');pt(P,'3a=3s_tag+27=Ee-(r_F+19) EXACTLY');pt(W,'(AG2) eta*r>=970(a+sig+2)');pt(P,'Rnew=(M0+p0*alpha+(23/4)*alpha+');pt(V,'16/15<=Rh<=64')
 for ss in sorted({s for k in range(160,201) for s in (2**k,2**(k+1)-1)}|{2**160+j for j in range(8)}):
  k=ss.bit_length()-1;s=2*ceil(((F(50,7)-F(1,2**100))*ss-lg(ss+4))/2);ellc=lg(s+4);rf=ss+walkextra(s);ab=rf+AB;u1=2**(110+lg(ss));rho=ab+2*lg(u1)+2
  walk=s+34*rf;tag=3*walk
  batch('AFexec.',[('rho',rho,rf+225+2*lg(ss),'=='),('rho-dyadic',rho-ss,10*k+267,'<='),('rho-bound',rho,(1+EB)*ss,'<='),('log',ellc,k+5,'<='),('walk',F(7,16)*F(ss,2**100),100+20*(k+1),'>'),('parity',walk%2,0,'==')])
  for L in (ceil(Q*ss),ceil(Q*(ss+1))-1):
   ee=3*tag+rf+22+lg(C);om=ee+13+2*L+42;d=floor(do*om);kk=ceil(aa*d+bb*rho)
   eeN=ee+13+4;np=d;adv=tag+2+lg(C)
   wn=max(2**80,eeN+44+2*lg(4*np));rlo=ceil(((PB+SD)*wn+DR+F(23,4)*eeN+56)/(theta-(6-eta/2)*x['c']));wf=max(2**80,eeN+2*lg(4*np)+32,ceil(F(DX*om,128)))
   batch('localexec.',[('Omega',om,1024*ss,'<='),('nlog',lg(4*np),42+lg(ss),'<='),('wN',wn,alpha*om,'<='),('wF',wf,alpha*om,'<='),('advice',3*adv,eeN,'<='),('advice-exact',3*adv,eeN-(rf+19),'=='),('AG1-term',PB*wn+3*eeN+F(41),adv+40,'>='),('AG2-term',eta*(PB+SD)*wn/(theta-(6-eta/2)*x['c']),970*(adv+eeN+4+2),'>='),('AG2-wN',eta*(PB+SD)*wn/(theta-(6-eta/2)*x['c']),F(3880,3)*wn+5820,'>'),('AG1-exec',floor((eta*rlo-11*(eeN+4)-8)/4),adv+40,'>='),('AG2-exec',eta*rlo,970*(adv+eeN+4+2),'>='),('hash',480+4*lg(om),NE*om,'<=')])
   for shift in (2,3):
    rho0=ab+L+shift+lg(kk)+2;b0=ceil(a0*rho0);u0=b0-rho0;b=ceil(gb*rho);uu=b-rho
    v0=b0*ceil(F(kk,u0));v1=b*ceil(F(kk,uu));cc=v1-kk;v2=b*ceil(F(cc,uu));D=b0+2*b
    ee=3*tag+rf+22+lg(C);om=ee+13+2*L+42
    batch('AFexec.',[('floor',floor(L/Q),ss,'=='),('rho0',rho0-(Q+1)*ss,9*k+Q+97,'<='),('rho0-bound',rho0,(Q+1+EB)*ss,'<='),('N0',rho0,N0,'>='),('u0',u0,1,'>='),('ratio0',F(b0,u0),N0+1,'<='),('u',uu,1,'>='),('ratio',F(b,uu),xi,'<='),('Kbound',kk,kb*ss,'<='),('Klog',lg(kk),48+lg(ss),'<='),('v0',v0,u1,'<='),('v1',v1,u1,'<='),('c',cc,1,'>='),('out',v1-v2,d,'>='),('rank-bound',v1-v2,gg*kk-xi*b,'>='),('seed',D,s,'<='),('source',kk,bh*L,'<='),('tag',tag,2*L,'>='),('runtime',tag,2*p*L/Q,'<='),('unit',om,co*L,'<='),('tail',8+2*ee,T*L,'<='),('combined',4*C*d+ko*om+8+2*ee,H*L,'<=')])
 batch('',[('SP.gamma',q_R-14*gamma+eta,NU,'=='),('SP.guard',223,NU*2**90,'<'),('SP.log-slope',2,223,'<='),('SP.gamma-fourteen',gamma,F(1,7),'<'),('SP.ring-log',29,F(2*2**90,7),'<'),
 ('SP.subset',(delta_R-gamma)*2**90,3,'>'),('SP.delta-cap',delta_R,F(1,7),'<'),('SP.projection-degree',q_R,1,'>'),('SP.hereditary-cap',h_R,7,'<'),('SP.dimension',12*gamma+F(29,2**90),F(2),'<'),
 ('SP.power',2+12*gamma+F(29,2**90),4,'<'),('SP.fiber',14*gamma+(h_R+F(5,2))*eta,2-NU,'=='),('SP.eta',(F(8249,485)+5*NU)*eta,gamma,'<='),('SP.cap',eta,F(1,120),'<='),('SP.theta',4*theta,eta,'=='),('CZ.quotient',4*F(5,4)+2*(F(486,485)+NU)+2*(F(1457,970)+NU)+7,F(8249,485)+4*NU,'=='),('CZ.exponent',4*F(5,4)+4*(F(486,485)+NU)+8,F(8249,485)+4*NU,'=='),('CZ.hered',(h_R-(F(486,485)+NU)-2)*64/NU,6,'>'),('CZ.zero',1-9*eta,0,'>'),('CZ.T',2+F(9,4)*eta,4,'<='),('CZ.smoothing-exponent',F(5,4),1,'>='),
 ('SP.cardinal',5+3+1,9,'=='),('SP.sixteen',9,16,'<'),('SP.endpoint',F(17+4)*3/8,F(63,8),'=='),('SP.endpoint-loss',F(63,8),8,'<'),
 ('CZ.units',1-(h_R+F(9,4))*eta,F(763,800),'>='),('CZ.plane',1-delta_R-(h_R+F(5,4))*eta,F(4569,5600),'>='),
 ('CZ.projection',2-(h_R+F(7,2))*eta-q_R,0,'>='),('CZ.subprojection',q_R-delta_R,1,'>'),('CZ.energy',2**28+2**26,2**29,'<'),('CZ.path',F(1,4)*F(1,64)**2,F(1,2**14),'=='),('CZ.path-bad',F(1,64)/F(1,4),F(1,16),'=='),('CZ.path-labels',2**14*4**4,2**22,'=='),
 ('CZ.smoothing',36+2,45,'<='),('CZ.hit',35,36,'<'),('CZ.root16',16**6,2**(7*64),'<'),
 ('prime.anchor',F(32,5),F(99,16),'>'),('prime.square',631,26**2,'<'),
 ('native.union',47*2*(F(2,3)+F(3,256))*2**16,2**22,'<'),('native.collision',F(2,2048)+F(1,4),F(257,1024),'=='),
 ('native.r',F(257,1020)+F(3,35*2**256),F(13,50),'<'),('native.mono',6-8*F(13,50),0,'>'),
 ('native.quartic',6*F(13,50)-4*F(13,50)**2-2,F(-444,625),'=='),('native.margin',F(-444,625)+F(1,1000)+F(6144,8925),F(-1,50),'<'),
 ('native.g',F(4,2**4095),F(1,1000),'<'),('native.delete',F(1,2**31),F(1,40000**2),'<'),('native.perturb',F(432,40000),F(1,80),'<'),
 ('native.moment',3+F(3,1024),4,'<='),('native.success',1-2*F(2,5),F(1,5),'=='),
 ('F.entropy',2*9+8,26,'=='),('F.seed2',2*9,18,'=='),('F.seed',4*9,36,'=='),('F.labels',1+2*2,5,'=='),('F.error',76,128,'<'),('F.round',52,76,'<='),('native.H128',H*128,2**128,'<')])
 ck('SP.constants',(4**4,1024*16*64,512*16*16,2**21*2**18,2**17*2**16),(256,2**20,2**17,2**39,2**33))
 ck('native.ratio',F(257,1024)*F(256,255),F(257,1020))
 ck('native.bias',F(4,47)+F(1,2**17),F(3,35),'<')
 ck('native.marginal',3+F(72*256,35*255),16,'<');ck('native.tail',F(432,40000),F(1,80),'<')
 ck('native.weight-count',F(2,3)+F(3,256),1,'<=');ck('native.moment',F(24*256,35*255),F(7,10),'<')
 ck('native.deletion',F(1,2**32)*F(256,255),F(1,2**31),'<')
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
 codes=['1100','1010','1001','0110','0101']
 pt(W,','.join(codes));ck('F.distinct',len(set(codes)),5)
 for aa in codes:
  first=aa.index('1');ck('F.prefix5',1+2*sum(bb[first]=='0' for bb in codes if bb!=aa),5,'<=')
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
  vec=v[nm];hh=len(vec);js=ar(nm,hh);BB=prod(js);aa=ee+lg(BB)+4;rho=aa+2*lg(4*nn)+2
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
  jj=js[-1];ag=nm not in ('Ext0','Ext2');tail=ell+4*aa+lg(4*(leaf+1))
  dd=sum(slots)+tail;en=envm(vec,ar(nm,len(vec))) if nm in ('Ext0','Ext2','I','F') else env(vec,ag)
  kap=(vec[-1]+1)/(vec[-1]-jj+1)
  ck('PT.seed-slices',sum(slots),sum(rr+1 for rr in vec)*ww+hh,'<=')
  if zz<=256*BB*ww:
   batch('CL.',[('ell',ell,256*ww,'<='),('kap',kap,31,'<='),('rate',vec[-1],64,'<='),('height',hh,13,'<='),('arity',jj,18,'<='),('mh',leaf,8192*ww,'<='),('log',lg(4*(leaf+1)),16+lg(ww),'<=')])
  ck('PX.threshold',kk,en['lam']*zz+en['g']*ww,'<=')
  if zz<=256*BB*ww:ck('PX.seed',dd,en['c']*zz+en['pi']*ww,'<=')
  if nm=='Ext0' and zz<=256*BB*ww:ck('PX.finite-seed',dd,en['c']*zz+PB*ww+DR,'<=')
  ck('PX.padding',kk,nn,'<=');return dd,kk
 for nm in v:
  for nn,zz,ee in ((2**256,1,1),(2**512,100000,100),(2**2048,1,2**22)):
   px(nm,nn,zz,ee)
 for nm in ('I','F','G','J','L0','Lr','E'):
  hh=len(v[nm]);px(nm,2**48,1,1,2**80);px(nm,2**256,256*prod(ar(nm,hh))*2**20,2**20-hh-522)
 UU=2**160;L=ceil(Q*UU);rf=UU;seed_s=2*ceil(((F(50,7)-F(1,2**100))*UU-lg(UU+4))/2);tag=3*(seed_s+34*rf)
 ee0=3*tag+rf+22+lg(C)+13;om=ee0+2*L+42
 nn=floor(do*om);zz=DX*om;adv=tag+2+lg(C);ee=ee0+4
 wn=max(2**80,ee+44+2*lg(4*nn));wf=max(2**80,ee+2*lg(4*nn)+32,ceil(F(zz,128)))
 m1,kf1=ceil(i['d']*wf),ceil(i['k']*wf);m0,k2=px('Ext2',nn,m1,ee)
 sig=ee+4;rr1=ceil((m0+PB*wn+DR+F(23,4)*ee+56)/(theta-(6-eta/2)*x['c']))
 ss=floor((eta*rr1-11*sig-8)/4);rawv=6*rr1-floor(eta*rr1/2);d01,k01=px('Ext0',nn,rawv,ee)
 el=ceil(F(d01+ee-1,2))+ee;d11=2*(lg(2*(rr1+adv+1)*m0+1)+m0+2*el);d2=d01+d11+m1+k2+ee
 batch('NBexec.',[('seed-lower',m0,SD*wn,'>='),('I1',ss,m0+2*el,'>='),('advice',ss,adv+40,'>='),('density',F(adv+sig,rr1+adv),eta/970,'<='),('AG1',ss,adv+40,'>='),('AG2',eta*rr1,970*(adv+sig+2),'>='),('AG1-margin',ss,500*(adv+40),'>'),('AG2-margin',10*eta*rr1,F(10001,1000)*970*(adv+sig+2),'>'),('source',rr1-sig,(1-eta/970)*(rr1+adv+1)+1,'>='),('joint-source',2*rr1-sig,(2-2*eta/970)*(rr1+adv+1)+1,'>='),('joint-short',2*rr1-floor(eta*rr1/2)-sig,(2-eta/2-eta/485)*(rr1+adv+1)+1,'>='),('capacity',F(5,4)*eta*(rr1+adv),5*ss+11*sig,'>='),('density-error',eta*(rr1+adv),4*ss+11*sig+8,'>='),('prime',rr1+adv,2**90,'>'),('dC',d01+d11+d2,floor(sc*om),'<='),('kC',max(k01,kf1+2*rawv+zz+2*ee),floor(kc*om),'<=')])
 UO=4*C*do+ko;e1=F(2599,7)-9*F(1,2**100);LA=e1/Q+2;A=UO*LA+2*e1/Q;B=F(0);CC=102*UO+102
 batch('refine.',[('source',A,do*LA/gg+Bbase/Q,'>'),('leading',A,H-1,'<='),('constant-guard',CC+1,2**160,'<'),('seed-margin',(3*ER-F(1,2**100))*2**160,165,'>='),('width-margin',(F(50,7)-F(1,2**100)-7)*2**160,162,'>='),('outer-log-margin',9*(UO+2)*160,CC+2,'>='),('walk-base',F(7,16)*2**60,3320,'>'),('walk-slope',20,3320,'<='),('alpha-constant',55*alpha-47*(1-alpha),0,'>'),('alpha-floor',alpha*(370+2*Q),1,'>')])
 ck('refine.constant-improvement',4*(CC+1),134976247231948,'<');ck('result.H',H,226903813738);ck('result.K',K,907615254952)
 ck('result.raw',raw,226903813738);pt(P,f'H_raw=4C*delta_O+kappa_O+T={raw}');ck('result.base',bh,226903813649)
 ck('result.improvement',4*A,F('83450049915423073911796719301359574178203876738232717/18731372364015809820563343750000000000000'),'<');ck('result.loss',[4*A],[F('275641530050784772536619301873011988211525204730975621989593898086872075586743/353603652556121189444782306923675269257838460928000000000000000000')])
 if TEXT:
  for key,val in [('c1',4*A)]:pt(P,key+'='+str(val))
  for nm,tag in [('I','MI'),('F','MF'),('G','MG'),('J','MJ'),('Ext0','NB0'),('Ext2','NB2'),('L0','L0'),('Lr','Lr'),('E','E')]:pt(V,' v'+tag+'=('+VE[nm]+')')
  for f in (P,V,S):pt(f,'Q_S='+str(Q))
  for tok in ('TR=TR0+H(t+1)(2dy+2dx)-(t+1)dx','T=TR+TL+(t+1)dx'):pt(W,tok)
  for f,tok in [(P,f'dy={DY}Om, dx={DX}Om, r={RATE}Om'),(P,f's_C=floor(({sc})Om)'),(P,f'k_Cslot=floor(({kc})Om)'),(P,f'd_O=floor(({do})Om), k_O=floor(({ko})Om)'),(P,f'={tr}Om'),(P,f'={tl}Om'),(P,f'K=4H={K}'),(V,f'B_B=ceil((xi_B*gamma_B+xi_B/2^40)/G_B)={Bbase}')]:pt(f,tok)
  sp=TEXT[S].split('# SP. Finite rings')[1]
  heads=re.findall(r'^#+ (SP\.[0-9.]+)',TEXT[S],re.M);tc('refs.headings',len(heads),len(set(heads)))
  labs=re.findall(r'^(?!#) *.*\((SP-[\w]+)\)\s*$',sp,re.M)
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
  terms=bj['guarantee']['remainder']['terms'];declared=[F(t['coefficient']['numerator'])/F(t['coefficient']['denominator']) for t in terms];new=[4*A]
  tc('meta.declared',declared==new)
  tc('meta.experiment',[e.get('source') for e in sj['experiments']],['verify.py'])
  tc('meta.evidence',all((ROOT/f).is_file() for claim in sj['claims'] for f in claim['evidence']))
  for claim in sj['claims']:tc('frozen.declaration',json.loads(claim['statement'].split('Declaration: ')[1]),bj)
 return dict(refined_loss=[str(4*A)],all_passed=not fail,all_checks_passed=not fail and not meta,checks=ct,arithmetic_checks=ac,text_checks_run=bool(TEXT),skipped_reason=None if TEXT else 'no package; arithmetic ran',H=str(H),K=str(K),H_raw=str(raw),H_base=str(bh),failures=fail,metadata_consistent=(not meta) if TEXT else None,metadata_failures=meta,groups=gs,scope='exact arithmetic+integrity')
HASHES={'sampler-and-disperser.md': '00006a1b4e00b53171e87ae96506559ac647b662fff9b2f4749f33078ec8e4f1', 'submission.json': '17f1509aeb52f56b1e5a845407734bdd2c18a51d42dcf34201a8b9e7dbc274f0', 'affine-and-standard-wrappers.md': '1067909c567c212bebecfca0736367f0d85527e13d80fe2245ac5b34137cf6da', 'restricted-nonmalleability-and-breaking.md': '474fc431bd36a6b9accf21e4578c507bb4e8dfff24883fcc4321ff825300cdad', 'proof.md': '336b799897c87565c4b60510e7d7064b62adc1ecb20cde204829a8460a138803', 'linear-and-bilinear-primitives.md': '7a92c16893a6102ab28e902d7e9034f7b1891f4515e1e53c1abb93b03b072780', 'bound.json': '96c76e470de347427e4deeba99d789edea9c3cc03e1d1ae954bb3ab6f38d9201'}
if __name__=='__main__':
 if json.load(sys.stdin)!={}:raise ValueError('Expected {}')
 roots=(Path(__file__).resolve().parent,Path('submission'),Path('/snapshot/files/submission'))
 ROOT=next((p for p in roots if all((p/f).is_file() for f in (P,W,V,S,R,'bound.json','submission.json'))),None)
 if ROOT:TEXT={f:(ROOT/f).read_text() for f in (P,W,V,S,R)}
 res=run();out=json.dumps(res,sort_keys=True);assert len(out.encode())<65536
 print(out);raise SystemExit(0 if res['all_passed'] else 1)
