from fractions import Fraction as F
from math import ceil,floor,isqrt,comb,prod
from pathlib import Path
import json,sys,re,hashlib
EPS=F(1,2**20);NU=F(1,2**80); NE=F(1,10**18); EB=F(1,10**18);ER=F(1,2**80);N0=2**60;AB=3;DG=F(18114171871,20000000);AP=F(1728,7)+F(1,10**6);AG=3*DG/2;ATG=AG+F(1,10**6)
SL={'a': '3475242658459537263/200000000000000000', 's': '14180938250242493107/1000000000000000000', 's0': '3475222875032359663/200000000000000000', 'hc': '110929330716693065469/1000000000000000000', 'az': '14161026002913383337/1000000000000000000', 'sz': '1407840376069793067/100000000000000000', 'd': '39792781701092586879/50000000000000000', 'k': '620177094375890235546081/250000000000000000'}
VN={'I': '2089387870654562,2089382159435300,2088850614824045', 'F': '2089422056574409,2089413827626494,2088883677848892', 'G': '3388867152809931', 'J': '3218689018843046', 'Jz': '2121824122564418', 'Gz': '1997890759241239', 'Fz': '2013609896154610,3013630897631181', 'Ext0': '28455564889762572,33136593627256344,33106860889940124,32939111506095664,32027739815380100,31214323648437332', 'Ext2': '16791013175945694,17540267024864918', 'L0': '3101630234180164,4223884039482512,4223621378833313,4220333218486985,4175100282751948', 'Lr': '3968992880314711,3980220352932104,5249934268502773,5236618836371856,5114402132460489', 'E': '2013211533595001,3013224931925304'}
VE={k:','.join(str(F(int(x),10**15)) for x in s.split(',')) for k,s in VN.items()}
QS=F('42787674712039/20000000000000');QC=17;DX=F(31463022350519,1099511627776) 
RSHIFT=F(1090064599267,1099511627776)
ZCAP=F(34386103315449,1000000000000)
INIT_GAMMA0=F(400177460987,200000000000)
INIT_XI0=INIT_GAMMA0/(INIT_GAMMA0-1);INIT_G0=1-(INIT_XI0-1)**2;INIT_A0=1/INIT_G0;INIT_B0=INIT_XI0*INIT_GAMMA0/INIT_G0+1
P='proof.md';W='affine-and-standard-wrappers.md';V='linear-and-bilinear-primitives.md'
S='sampler-and-disperser.md';R='restricted-nonmalleability-and-breaking.md'
TAG=dict(I='MI',F='MF',G='MG',J='MJ',Ext0='NB0',Ext2='NB2',L0='L0',Lr='Lr',E='E',Fz='MFZ',Gz='MGZ',Jz='MJZ')
TEXT={};ALIASES={};ROOT=None;fail=[];meta=[];ct=0;ac=0;gs={}
def ck(nm,a,b=True,op='==',text=False):
 global ct,ac
 ok=(a==b if op=='==' else a<b if op=='<' else a>b if op=='>' else a<=b if op=='<=' else a>=b)
 ct+=1;ac+=not text;key=nm.split('.')[0];gs[key]=gs.get(key,0)+1
 if not ok:(meta if nm.startswith('meta.') else fail).append(dict(check=nm,left=str(a),right=str(b),op=op))
def batch(prefix,rows):
 for nm,*a in rows:ck(prefix+nm,*a)
def tc(nm,*a):ck(nm,*a,text=True)
def pt(f,tok):
 for v,k in ALIASES.items():tok=tok.replace(v,k)
 if TEXT:tc('text.'+tok[:70],' '.join(tok.split()) in ' '.join(TEXT[f].split()))
def lg(x):
 x=F(x);k=max(0,x.numerator.bit_length()-x.denominator.bit_length());return k+(2**k<x)
JS=(4,5,5,5,5,6);JY=(5,6);JI=(3,3,3);JG=(3,);JL0=(3,4,4,4,4);JLR=(3,3,4,4,4);JE=(3,4);JFZ=(3,4);OM=2**40;BM=prod(JS);BY=prod(JY)
def envm(vec,js=JS,acc=F(1)):
 assert len(vec)==len(js)
 p=F(1);g=F(0)
 for ix,(j,r) in enumerate(zip(js,vec)):
  lam=j*r/(r-j+1);g+=p*(acc*lam+(1 if ix==0 else acc)*F(j-1,j)*(r+1)+max(0,(j-7)*lam+F(j-1,j))*EPS);p*=lam
 B=prod(js)
 return dict(lam=p/B,g=g+2*acc*p,c=F(1,B),pi=4*acc+vec[0]+1+acc*sum(r+1 for r in vec[1:]))
def ar(nm,h):return JS if nm=='Ext0' else JY if nm=='Ext2' else JI if nm in ('I','F') else JG if nm=='G' else JL0 if nm=='L0' else JLR if nm=='Lr' else JE if nm=='E' else JFZ if nm=='Fz' else (2,)*h
def et(e):return ','.join(k+'='+str(e[k]) for k in ('c','pi','lam','g'))
def omunit(x):return x
def walkextra(s):return 0
def slot(x):return F(ceil(x*10**18)+2,10**18)
def cleaf(N,m,a):
 D=1<<(2*a+1);B=D+1;d=m+4*a+lg(4*(N+1));T=1<<d
 xs=list(range(1,1<<N));ww=[1]*len(xs);tb=[]
 def pot(rs):
  tot=0
  for x,w in zip(xs,ww):
   fr=0;bd=False
   for row in rs:
    fx=0;ud=False
    for j,v in enumerate(row):
     if x>>j&1:
      if v is None:ud=True
      else:fx^=v
    if ud:fr+=1
    elif fx:bd=True;break
   tot+=w*((D<<m)+(0 if bd else 1<<(m-fr)))
  return tot
 for s in range(T):
  rs=[[None]*N for _ in range(m)]
  for i in range(m):
   for j in range(N):
    rs[i][j]=0;p0=pot(rs);rs[i][j]=1;p1=pot(rs);rs[i][j]=int(p1<p0)
  ms=[sum(v<<j for j,v in enumerate(row)) for row in rs];tb.append(ms)
  for t,x in enumerate(xs):ww[t]*=B if all((r&x).bit_count()%2==0 for r in ms) else D
 co=[sum(all((r&x).bit_count()%2==0 for r in ms) for ms in tb) for x in xs]
 assert max(co)<T*F(1,1<<m)*(1+F(1,1<<(2*a)))
 return tb,co

def run():
 table,counts=cleaf(5,1,1)
 ck('CLtable.size',len(table),1024)
 for nzero in counts:ck('CLtable.collision',nzero,640,'<')
 for source in (list(range(8)),[0,1,2,4,8,16,3,31],[1,3,5,7,9,11,13,15]):
  tvnum=[abs(2*sum((row[0]&xx).bit_count()%2 for xx in source)-len(source)) for row in table];excess=sum(tvnum)
  ck('CLtable.strong',F(excess,2*len(source)*len(table)),F(1,2),'<');ck('CLtable.weakseed',F(sum(sorted(tvnum,reverse=True)[:len(table)//2]),len(source)*len(table)),F(1,2),'<=')
 v={k:list(map(F,s.split(','))) for k,s in VE.items()};i={k:F(s) for k,s in SL.items()}
 acc=(F(1728,7)-6*F(1,2**100))/(F(1728,7)-6*F(1,2**100)+QS)+NE
 ev={k:envm(a,ar(k,len(a)),acc if k in ('L0','Lr') else F(1)) for k,a in v.items()}
 ev['E']['g']*=acc;ev['E']['pi']*=acc
 I,FF,G,J=[ev[k] for k in ('I','F','G','J')];FZ,GZ,JZ=[ev[k] for k in ('Fz','Gz','Jz')];x,y=ev['Ext0'],ev['Ext2']
 SD=(1-F(1,2**40))*(sum(rr+1 for rr in v['Ext2'])+i['d']/BY+4);DR=F(0);PB=x['pi']-DR*EPS
 h_R=F(3)+F(2,DG)+2*NU;DD=F(17)+F(8,DG)+5*NU;eta=(2-NU)/(14*DD+h_R+F(5,2));theta=eta/4;gamma=DD*eta;delta_R=gamma+NU;q_R=14*gamma-eta+NU;U=2**40
 eb=F(1728,7)-6*F(1,2**100);alpha=eb/(eb+QS)+NE
 df=i['d']*alpha+NE;kf=i['k']*alpha+NE
 M0=y['c']*df+y['pi']*alpha+NE;K2=y['lam']*df+y['g']*alpha+NE
 rr=(M0+PB*alpha+F(23,4)*alpha+(56+DR)*NE)/(theta-(6-eta/2)*x['c'])+NE
 amin=F(12342857,100000)/(F(1728,7)+QS+F(1,10**6));rr-=(1+(6-eta/2)*x['c']/(theta-(6-eta/2)*x['c']))*amin
 d0=(6-eta/2)*x['c']*rr+PB*alpha+(DR+2)*NE;d1=M0+2*d0+6*alpha+NE
 k00=(6-eta/2)*x['lam']*rr+x['g']*alpha+NE;kother=kf+10*M0+DX+2*alpha+NE
 dc=2*d0+2*d1+df+K2+alpha;sc=slot(dc);kc=slot(max(k00,kother))
 q=QC;C=4*(q+1);t=4*C-1;j=lg(t);Q=QS;beta=F(2*(q+8),7);p=144
 EE,RR,ZZ=[ev[k] for k in ('E','Lr','L0')];RR=dict(c=F(0),pi=1+8*alpha,lam=F(9,8),g=7*alpha);ZZ=dict(c=F(0),pi=1+2*INIT_GAMMA0*alpha,lam=INIT_A0,g=INIT_B0*alpha);RR=dict(c=F(0),pi=1+8*alpha,lam=F(9,8),g=7*alpha);ZZ=dict(c=F(0),pi=1+2*INIT_GAMMA0*alpha,lam=INIT_A0,g=INIT_B0*alpha);RR=dict(c=F(0),pi=1+8*alpha,lam=F(9,8),g=7*alpha);ZZ=dict(c=F(0),pi=1+2*INIT_GAMMA0*alpha,lam=INIT_A0,g=INIT_B0*alpha)
 DY=F(14641446576009,1099511627776);RATE=(t+2)*DY+RSHIFT
 s0=slot(ZZ['c']*sc+ZZ['pi']);tr=(t+1)*(s0+sc+2*j*(DY+DX));tl=(t+1)*(sc+2*j*DY)
 rows=dict(dy_E=EE['c']*DX+EE['pi'],dx_Lr=RR['c']*RATE+RR['pi'],dx_Lm=2*DY+2*alpha+NE,s0_L0=ZZ['c']*sc+ZZ['pi'],r_3a=(t+2)*DY+alpha,r_dy=DY,r_Km=DY+2,
 d_3c=(t+1)*(s0+sc)+kc+alpha,d_3d=tr+EE['lam']*DX+EE['g']+t*DX+alpha,d_s0=s0,d_KE=EE['lam']*DX+EE['g'],d_kC=kc,
 k_3b=tl+RR['lam']*RATE+RR['g']+t*RATE,k_K0=ZZ['lam']*sc+ZZ['g'])
 do=slot(max(val for key,val in rows.items() if key.startswith('d_')));ko=slot(max(val for key,val in rows.items() if key.startswith('k_')))
 ALIASES.update({str(DX):'DX',str(DY):'DY',str(RATE):'RATE'});ALIASES.update({str(do):'Dslot',str(ko):'Kslot',str(s0):'S0'})
 co=(3*p+1)/Q+1+F(1,10**15);delta=ceil(do*co);kappa=ceil(ko*co);T=1
 pt(P,f'T={T}')
 a0=F(N0+1,N0);gb=(beta-a0*(Q+1+EB)-3*ER)/(2*(1+EB));xi=gb/(gb-1);gg=1-(xi-1)**2;aa=1/gg;bb=ceil((xi*gb+xi/U)/gg)
 Bbase=bb;raw=4*C*delta+kappa+T;bh=ceil(aa*delta+bb*(1+EB)/Q+ER);H=max(raw,bh);K=4*H
 kb=aa*delta*Q*(1+F(1,U))+bb*(1+EB)+ER;b0c=a0*(Q+1+EB)+ER
 IR=lambda I,FF,G,J,FZ,GZ,JZ:[('aJ','a',J['c']*(i['s']+NE)+J['pi']),('aF','a',FF['c']*(i['hc']+NE)+FF['pi']),('sG','s',G['c']*(i['a']+NE)+G['pi']),('s0I','s0',I['c']*(i['hc']+NE)+I['pi']),('hJ','hc',(5+J['lam'])*(i['s']+NE)+J['g']+1),('hJz','hc',(5+JZ['lam'])*(i['sz']+NE)+JZ['g']+1),('hcap','hc',F(100)),('h>=s','hc',i['s']+NE),('azF','az',FZ['c']*ZCAP+FZ['pi']),('azJz','az',JZ['c']*(i['sz']+NE)+JZ['pi']),('szGz','sz',GZ['c']*(i['az']+NE)+GZ['pi']),('sz<=s','s',i['sz']+NE),('az<=a','a',i['az']+NE),('dEarly','d',9*(i['s0']+NE)+(27+G['lam'])*(i['a']+NE)+G['g']+1),('dLast','d',9*(i['s0']+NE)+27*(i['a']+NE)+(9+GZ['lam'])*(i['az']+NE)+GZ['g']+1),('d>=s0','d',i['s0']+NE),('kFearly','k',(22+FF['lam'])*(i['hc']+NE)+13*(i['s']+NE)+FF['g']+1),('kFlast','k',18*(i['hc']+NE)+18*(i['s']+NE)+(8+FZ['lam'])*ZCAP+FZ['g']+1),('kI','k',I['lam']*(i['hc']+NE)+I['g']),('k>=d','k',i['d']+NE)];ir=IR(I,FF,G,J,FZ,GZ,JZ)
 for nm,ky,val in ir:
  ck('inner.'+nm,i[ky],val,'>=')
 for ky,val in i.items():pt(W,f'{ky}=ceil(({val})w).')
 for ky,val in rows.items():
  lim=dict(dy=DY,dx=DX,r=RATE,s0=s0,d=do,k=ko)[ky.split('_')[0]]
  ck('outer.'+ky,lim,(F(ceil(val*10**18),10**18)+NE if ky.split('_')[0] in ('s0','d','k') else val),'>=');pt(P,f'{ky}: '+('U18(rhs)+eps_N' if ky.split('_')[0] in ('s0','d','k') else 'rhs')+f' <= {lim}.')
 den=theta-(6-eta/2)*x['c'];runit=(i['d']/BY+y['pi']+PB+6)/den
 batch('Zcap.',[('Q',QS,F(21393851,10**7),'<'),('dyadic',10**6*9*162,2**160,'<'),('E',F(1728,7)-6*F(1,2**100)-F(1,10**6),F(1728,7)-F(2,10**6),'>'),('Om',(OM+116)*10**6,2**160,'<'),('z',DX*(F(1728,7)+QS+F(1,10**6)),ZCAP*(F(1728,7)-F(2,10**6)),'<')])
 pt(P,'Q_S<3 give');pt(P,'Om<=(1728/7+Q_S)S+M+116<=(1728/7+Q_S+1/10^6)S,');pt(P,f'({DX})*(1728/7+Q_S+1/10^6)<({ZCAP})*(1728/7-2/10^6).')
 batch('NBshift.',[('beta',3*(F(50,7)-F(1,2**100))+102,F(12342857,100000),'>='),('capacity',eta/4+eta/DG,1,'<'),('amin',amin,0,'>')])
 batch('CLguard.',[('Om',co*Q*(1+F(1,2**80)),440,'<='),('wF',DX*440,128*246,'<='),('Fguard',(100*246-DX*440)*2**160,1200,'>'),('Ext2cap',i['d']+NU,256*BY,'<='),('Ext0cap',(6-eta/2)*runit+NU,256*BM,'<='),('rround',56+den+F(1,BY),F(2**80,4),'<=')])
 for nm,vec in v.items():
  en=ev[nm];h=len(vec)
  for ridx,rate in enumerate(vec):
   ck('PT.rate',rate,ar(nm,len(vec))[ridx]-1,'>')
   sh=ceil(7*(rate+1));nn=ceil(F(3*(sh-14),2))-2
   bc=ceil(F(2*(nn+2),sh-14));mm=ceil(F(bc*sh,2));kk=mm+nn+1
   ck('PT.boundary',kk,bc*(sh-7),'<')
  f=W if nm in ('I','F','G','J','Fz','Gz','Jz') else P
  if nm in ('I','F','G','J','Fz','Gz','Jz'):pt(W,'ordered rates,pi=4+sum_i(R_i+1)')
  kap=(vec[-1]+1)/(vec[-1]-ar(nm,len(vec))[-1]+1)
  batch('CLguard.',[('kap',kap,505,'<='),('last',vec[-1],64,'<='),('height',h,13,'<=')])
  if nm not in ('Ext0','Ext2'):
   out=dict(I=i['hc']+EPS,F=i['hc']+EPS,G=i['a']+EPS,J=i['s']+EPS,Fz=ZCAP,Gz=i['az']+EPS,Jz=i['sz']+EPS,L0=sc,Lr=RATE,E=DX)[nm]
   ck('CLguard.output',out,256*prod(ar(nm,h)),'<=')
 pt(S,f'D_E={DD-5*NU}+5nu');pt(S,f'relative energy exponent{DD-5*NU}+4nu');pt(S,f'h_R={h_R-2*NU}+2nu');pt(R,f'a_set=c_set={1+F(2,DG)}+nu,b_set={F(3,2)+F(2,DG)}+nu');pt(R,f'2-({F(1,2)+F(2,DG)})eta')
 pt(V,'## 7h. Mixed-arity exact primitive PXTM');pt(W,'Ext0 is PXTM with6 mixed-arity levels');pt(W,'Ext2=PXTM(JY,vNB2)');pt(P,'bases 1067<')
 pt(W,'p0='+str(PB));pt(W,'D0='+str(DR));pt(W,'1000000000000z<=34386103315449w')
 pt(W,'c0='+str(x['c']));pt(W,'pi0='+str(x['pi']))
 pt(W,'bc0='+str((6-eta/2)*x['c'])+'<theta='+str(theta));pt(P,f'delta_O=ceil({do}c_O)={delta}');pt(P,f'kappa_O=ceil({ko}c_O)={kappa}');pt(P,'D_F=alpha*d_inner+eps_N');pt(P,'K_F=alpha*k_inner+eps_N');pt(P,str(co))
 cl=[ceil(z) for z in (M0,K2,rr,d0,d1,dc)];pt(P,'ceil(M0)=%d, ceil(K2bar)=%d,\n ceil(Rnew)=%d, ceil(d0bar)=%d,\n ceil(d1bar)=%d, ceil(DNB)=%d.'%tuple(cl))
 pt(P,'KC='+str(kc));pt(P,'SC=Slot(DNB)='+str(sc))
 eu={k:{c:F(ceil(e[c]*10**40),10**40) for c in e} for k,e in ev.items()}
 for c in ('c','pi','lam','g'):
  for k in ev:ck('AE.enclosure',0<=eu[k][c]-ev[k][c]<F(1,10**40),True)
 xu,yu=eu['Ext0'],eu['Ext2'];M0u=yu['c']*df+yu['pi']*alpha+NE;K2u=yu['lam']*df+yu['g']*alpha+NE
 ck('AE.denominator',theta,(6-eta/2)*xu['c'],'>');ru=(M0u+(xu['pi']-DR*EPS)*alpha+F(23,4)*alpha+(56+DR)*NE)/(theta-(6-eta/2)*xu['c'])+NE-(1+(6-eta/2)*xu['c']/(theta-(6-eta/2)*xu['c']))*amin
 d0u=(6-eta/2)*xu['c']*ru+(xu['pi']-DR*EPS)*alpha+(DR+2)*NE;d1u=M0u+2*d0u+6*alpha+NE;dcu=2*d0u+2*d1u+df+K2u+alpha
 aec=[ceil(z) for z in (M0u,K2u,ru,d0u,d1u,dcu,(6-eta/2)*xu['lam']*ru+xu['g']*alpha+NE,kf+10*M0u+DX+2*alpha+NE)]
 ck('AE.seed',dcu,sc-NE,'<=');ck('AE.source',max((6-eta/2)*xu['lam']*ru+xu['g']*alpha+NE,kf+10*M0u+DX+2*alpha+NE),kc-NE,'<=');
 Iu,Fu,Gu,Ju=[eu[k] for k in ('I','F','G','J')]
 for nm,ky,val in IR(Iu,Fu,Gu,Ju,eu['Fz'],eu['Gz'],eu['Jz']):ck('AE.inner.'+nm,i[ky],val,'>=')
 Eu=eu['E'];Ru={k:F(ceil(v*10**40),10**40) for k,v in RR.items()};Zu={k:F(ceil(v*10**40),10**40) for k,v in ZZ.items()}
 for nm,lim,val in [('dyE',DY,Eu['c']*DX+Eu['pi']),('dxLr',DX,Ru['c']*RATE+Ru['pi']),('s0L0',s0,Zu['c']*sc+Zu['pi']),('d3d',do,tr+Eu['lam']*DX+Eu['g']+t*DX+alpha),('dKE',do,Eu['lam']*DX+Eu['g']),('k3b',ko,tl+Ru['lam']*RATE+Ru['g']+t*RATE),('kK0',ko,Zu['lam']*sc+Zu['g'])]:ck('AE.outer.'+nm,lim,val+(NE if nm in ('s0L0','d3d','dKE','k3b','kK0') else 0),'>=')
 for val in (k00,kother):ck('NB.source',val,kc,'<=')
 ck('NB.seed',dc,sc,'<=');ck('NB.source-dimension',kc,do,'<=');ck('NB.seed-dimension',sc,do,'<=')
 batch('NB.',[('denominator',theta,(6-eta/2)*x['c'],'>'),('p0',PB,16,'>='),('height',max(lg(BM),lg(BY)),38,'<='),('cap',DX,128,'<='),('advice',F(23,4)+F(56,U),28,'<='),('unit-p0',4*PB,64,'>='),('AG2-slope',AG,ATG,'<'),('AP-base',418,F(2**80,10**6),'<'),('AP-growth',2,418,'<='),('Ee-base',290,F(2**80,2**40),'<'),('Ee-growth',2,290,'<='),('AG2-const',(ATG-AG)*2**80,(6*DG),'>'),('AG2-r-lower',eta*(PB+SD+F(23,4)*(1-F(1,2**40))-(6-eta/2)*x['c']/2)/(theta-(6-eta/2)*x['c']),ATG,'>='),('AG2-seed-round',33,F(2**80,2**40),'<='),('seed-aa-base',441,F(2**160,2**40),'<'),('seed-aa-slope',2,441,'<='),('constant',ceil(2*(rr+2)*M0)+1,2**240,'<'),('weaklog',4*(6+1),(6+1)**2,'<='),('Lm',RATE*DY+1,2**33,'<'),('heights',max(lg(prod(ar(k,len(v[k])))) for k in ('I','F','G','J'))+13,21,'<='),('outer-heights',max(lg(prod(ar(k,len(v[k])))) for k in ('L0','Lr','E'))+12,25,'<='),('NB-heights',max(lg(BM),lg(BY))+6,27,'<='),('prime',4*16*2**80/eta,2**90,'>'),('capacity',F(5,4)*4,5,'=='),('raw-blocks',1+1+(2-eta/2)+2,6-eta/2,'=='),('raw-round',max(4,x['lam']),NE*2**80,'<='),('c0',x['c'],1,'<'),('raw-rows',4+1,5,'=='),('matching-charge',2*4+2,10,'=='),('sig-capacity',F(5,4)*11,11,'>'),('floor',56-x['c']-13-1,F(41),'>='),('advice40',F(41),40,'>'),('error',14,16,'<'),('final-error',15,16,'<')])
 for hh in (0,1,2,j):
  rt=s0+sc+DX;lf=sc;mx=0
  for _ in range(hh):
   rt+=DY;lf+=DY;mx=max(mx,rt);rt+=DX+DY;lf+=DY;mx=max(mx,rt);rt+=DX
  if hh:ck('causal.E',mx,s0+sc+2*hh*(DY+DX))
  ck('causal.full',lf+rt,s0+2*sc+DX+hh*(4*DY+2*DX))
 batch('AF.',[('ceil-round',1,ER*2**80,'<='),('gamma',gb,2,'>'),('rank',gg*aa,1,'=='),('rounding',gg*bb-xi*gb,xi/U,'>='),('B-minimal',gg*(bb-1)-xi*gb,xi/U,'<'),('seed',a0*(Q+1+EB)+2*gb*(1+EB)+3*ER,beta,'=='),('Kbits',kb,2**48,'<'),('pad0',(N0+1)*2**48+b0c,2**110,'<'),('pad1',xi*2**48+gb*(1+EB)+EB,2**110,'<'),('error',6,2**AB,'<'),('q',F(7,2)*beta-8,q,'=='),('p',p,3*Q,'>'),('H',H,2*delta,'>='),('runtime-power',ceil(2*p/Q),135,'=='),('incidence',48*U,134217728,'>')])
 for nm,c,b,budget in [('rho',267,10,EB),('rho0',Q+97,9,EB),('outer',10*5768+160,2*5768,EPS*Q),('runtime',48*q*5+6,48*q,p),('NB-log',248,8,NE),('NB-quarter',242,2,F(1,4)),('Lm',84,4,NE)]:
  kk=80 if nm in ('rho','rho0','NB-log','Lm') else 40;ck('induction.'+nm,c+kk*b,budget*2**kk,'<');ck('induction.slope',b,c+kk*b,'<=')
 ck('induction.Om-base',5768*165+OM+107,F(1,10**15)*Q*2**160,'<');ck('induction.Om-slope',5768,5768*165+OM+107,'<=')
 batch('Leaf505.',[('binary-endpoint',(F(253,252)+1)/(F(253,252)-1),505),('mixed-endpoint',(F(477,28)+1)/(F(477,28)-17),505),('large-gap',7*20-13-5,0,'>'),('small-gap',20+20+13+1,2*2**20,'<')]);ck('PXG.legal-small',F(1,8192)+15*EPS,2,'<')
 batch('PXG.',[('log-base',64,64,'<='),('log-slope',2,64,'<='),('mh',505*(256+3)+F(65,2),130828,'<'),('dimension',130828*2**20+9,131072*2**20,'<='),('product',131072*256,2**25,'==')])
 batch('local.',[('alpha',alpha,1,'<'),('output',F(DX,128)+NE,alpha,'<='),('source',do,2**30,'<'),('Omega',co*Q*(1+F(1,2**80)),1024,'<'),('unit-base',490628,NE*(eb+Q)*2**80,'<'),('unit-slope',5770,490628,'<='),('hash-base',804,NE*2**80,'<'),('hash-slope',4,804,'<='),('round',1,NE*2**80,'<='),('advice-charge',F(23,4),F(23,4)),('AG1-coef',9*F(1,3),1,'>='),('AG2-coef',DG*F(1,2)+DG,F((3*DG),2),'==')])
 pt(P,f'R0={RSHIFT}');pt(P,'eps=2^-18');pt(P,f'z=ceil(({DX})Om)');pt(P,'5770*80+29028=490628');pt(P,'2a=2s_tag+18=Ee-16 EXACTLY');pt(W,'(AG2) eta*(r+a)>=(18114171871/20000000)(a+sig+2)');pt(P,'Rnew=(M0+p0*alpha+(23/4)*alpha+');pt(V,'253/252<=Rh<=64')
 for ss in sorted({s for k in range(160,201) for s in (2**k,2**(k+1)-1)}|{2**160+j for j in range(8)}):
  k=ss.bit_length()-1;s=2*ceil(((F(50,7)-F(1,2**100))*ss)/2);ellc=lg(s+4);rf=ss+walkextra(s);ab=rf+AB;u1=2**(110+lg(ss));rho=ab+2*lg(u1)+2
  walk=s+34*rf;tag=3*walk
  batch('AFexec.',[('rho',rho,rf+225+2*lg(ss),'=='),('rho-dyadic',rho-ss,10*k+267,'<='),('rho-bound',rho,(1+EB)*ss,'<='),('log',ellc,k+5,'<='),('walk',F(7,16)*F(ss,2**100),100+20*(k+1),'>'),('parity',walk%2,0,'==')])
  for L in (ceil(Q*ss),ceil(Q*(ss+1))-1):
   ee=2*tag+12;om=omunit(ee+18+L+37);d=floor(do*om);kk=ceil(aa*d+bb*rho)
   eeN=ee+18+4;np=d;adv=tag+2+lg(C)
   wn=max(2**80,eeN+44+2*lg(4*np));rlo=ceil(((PB+SD)*wn+DR+F(23,4)*eeN+56-(6-eta/2)*x['c']*adv)/(theta-(6-eta/2)*x['c']));wf=max(2**80,eeN+2*lg(4*np)+32,ceil(F(DX*om,128)))
   batch('localexec.',[('Omega',om,1024*ss,'<='),('nlog',lg(4*np),42+lg(ss),'<='),('wN',wn,alpha*om,'<='),('wF',wf,alpha*om,'<='),('advice',2*adv,eeN,'<='),('advice-exact',2*adv,eeN-16,'=='),('AP',wn,AP*ss,'<='),('Ee-lower',eeN,(1-F(1,2**40))*wn,'>='),('AG2-upper',DG*(adv+eeN+6),AG*wn,'<='),('AG1-term',PB*wn+3*eeN+F(41),adv+40,'>='),('AG2-term',eta*(PB+SD+F(23,4)*(1-F(1,2**40))-(6-eta/2)*x['c']/2)*wn/(theta-(6-eta/2)*x['c']),DG*(adv+eeN+4+2),'>='),('AG2-wN',eta*(PB+SD+F(23,4)*(1-F(1,2**40))-(6-eta/2)*x['c']/2)*wn/(theta-(6-eta/2)*x['c']),AG*wn+(6*DG),'>'),('amin-exec',adv,amin*om,'>='),('raw-positive',rlo-adv,1,'>='),('raw-capacity',rlo-adv,floor((eta*rlo-11*(eeN+4)-8)/4)+eeN+4,'>='),('AG1-exec',floor((eta*rlo-11*(eeN+4)-8)/4),adv+40,'>='),('AG2-exec',eta*rlo,DG*(adv+eeN+4+2),'>='),('hash',480+4*lg(om),NE*om,'<=')])
   for shift in (2,3):
    rho0=ab+L+shift+lg(kk)+2;b0=ceil(a0*rho0);u0=b0-rho0;b=ceil(gb*rho);uu=b-rho
    v0=b0*ceil(F(kk,u0));v1=b*ceil(F(kk,uu));cc=v1-kk;v2=b*ceil(F(cc,uu));D=b0+2*b
    ee=2*tag+12;om=omunit(ee+18+L+37)
    batch('AFexec.',[('floor',floor(L/Q),ss,'=='),('rho0',rho0-(Q+1)*ss,9*k+Q+97,'<='),('rho0-bound',rho0,(Q+1+EB)*ss,'<='),('N0',rho0,N0,'>='),('u0',u0,1,'>='),('ratio0',F(b0,u0),N0+1,'<='),('u',uu,1,'>='),('ratio',F(b,uu),xi,'<='),('Kbound',kk,kb*ss,'<='),('Klog',lg(kk),48+lg(ss),'<='),('v0',v0,u1,'<='),('v1',v1,u1,'<='),('c',cc,1,'>='),('out',v1-v2,d,'>='),('rank-bound',v1-v2,gg*kk-xi*b,'>='),('seed',D,s,'<='),('source',kk,bh*L,'<='),('tag',tag,2*L,'>='),('runtime',tag,2*p*L/Q,'<='),('unit',om,co*L,'<='),('tail',8,T*L,'<='),('combined',4*C*d+ko*om+8,H*L,'<=')])
 batch('',[('SP.gamma',q_R-14*gamma+eta,NU,'=='),('SP.guard',223,NU*2**90,'<'),('SP.log-slope',2,223,'<='),('SP.gamma-fourteen',gamma,F(1,7),'<'),('SP.ring-log',29,F(2*2**90,7),'<'),
 ('SP.subset',(delta_R-gamma)*2**90,3,'>'),('SP.delta-cap',delta_R,F(1,7),'<'),('SP.projection-degree',q_R,1,'>'),('SP.hereditary-cap',h_R,7,'<'),('SP.dimension',12*gamma+F(29,2**90),F(2),'<'),
 ('SP.power',2+12*gamma+F(29,2**90),4,'<'),('SP.fiber',14*gamma+(h_R+F(5,2))*eta,2-NU,'=='),('SP.eta',((F(17)+F(8,DG))+5*NU)*eta,gamma,'<='),('SP.cap',eta,F(1,120),'<='),('SP.theta',4*theta,eta,'=='),('CZ.quotient',4*F(5,4)+2*((F(1)+F(2,DG))+NU)+2*((F(3,2)+F(2,DG))+NU)+7,(F(17)+F(8,DG))+4*NU,'=='),('CZ.exponent',4*F(5,4)+4*((F(1)+F(2,DG))+NU)+8,(F(17)+F(8,DG))+4*NU,'=='),('CZ.hered',(h_R-((F(1)+F(2,DG))+NU)-2)*64/NU,6,'>'),('CZ.zero',1-9*eta,0,'>'),('CZ.T',2+F(9,4)*eta,4,'<='),('CZ.smoothing-exponent',F(5,4),1,'>='),
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
 batch('rank.',[('bad-density',2,5,'<'),('budget',F(2**9,2**12)*(F(2,3)+F(3,256)),F(3*255,35*256),'<')])
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
 def px(nm,nn,zz,ee,unit=0,cd=None,lu0=None):
  w=v[nm];h=len(w);q=ar(nm,h);BB=prod(q);e=ee+lg(BB)+(2 if cd is not None else 4)
  l=cd is not None;r=([e+lu0+22]+[e+2*lg(4*cd)+2]*(h-1)) if l else [e+2*lg(4*nn)+2]*h
  u=unit if l else max(unit,2**20,r[0]);o=ceil(F(zz,BB));xi=alpha if l else F(1)
  g=[F(1)]+[xi]*(h-1)
  s=[ceil((rr+1)*rh) for rr,rh in zip(w,r)];kk=o+2*e
  if l:
   ck('LC.first',r[0]-e-lu0,22,'>=');ck('LC.gap',4*(xi*u-e),20+lg(u)+h+1,'>=')
   for ix,rh in enumerate(r):
    ck('LC.rho',rh,g[ix]*u,'<=');ck('LC.low',u,2*rh,'<=')
    if ix:ck('LC.later',rh-e-2*lg(4*cd),2,'>=')
  ck('PT.leaf-threshold',kk,F(zz,BB)+2*u-11,'<=')
  if not l:ck('PT.expansion-guard',4*(4*nn)**2,2**(r[0]-e),'<=')
  ck('PT.a-unit',e,xi*u-6,'<=')
  for ix in range(h-1,-1,-1):
   p=r[ix];a=kk;j=q[ix];b=ceil(F(j*(a+e+j-1),s[ix]-j*p));m=ceil(F(b*s[ix],j));kk=(j-1)*m+a+e
   ck('PX.capacity',kk,b*(s[ix]-p),'<');ck('PX.count',b,kk,'<=')
   ck('PX.split',kk,(j-1)*m+a+e,'==');ck('PX.child',m,a,'>')
   if l:
    ck('LC.bcap',b,2**20,'<=')
    if ix:ck('LC.padding',m,cd,'<=')
   if ix==h-1:
    f=m;k=(w[-1]+1)/(w[-1]-j+1)
    ck('PT.terminal',f,k*(ceil(F(zz,BB))+3*e+j-1)+F(s[-1],j)+1,'<')
  j=q[-1];t=o+4*e+lg(4*(f+1))
  d=sum(s)+t;n=envm(w,ar(nm,len(w)),xi)
  k=(w[-1]+1)/(w[-1]-j+1)
  ck('PT.seed-slices',sum(s),sum(g[ix]*(rr+1) for ix,rr in enumerate(w))*u+h,'<=')
  if zz<=256*BB*u:
   batch('CL.',[('ell',o,256*u,'<='),('kap',k,505,'<='),('rate',w[-1],64,'<='),('height',h,13,'<='),('arity',j,18,'<='),('mh',f,131072*u,'<='),('log',lg(4*(f+1)),20+lg(u),'<=')])
  ck('PX.threshold',kk,n['lam']*zz+n['g']*u,'<=')
  if zz<=256*BB*u:ck('PX.seed',d,n['c']*zz+n['pi']*u,'<=')
  if nm=='Ext0' and zz<=256*BB*u:ck('PX.finite-seed',d,n['c']*zz+PB*u+DR,'<=')
  if l:
   ck('LC.Kd',kk,cd,'<=');ck('LC.Ku',lg(kk),lu0-3,'<=')
  else:ck('PX.padding',kk,nn,'<=')
  return d,kk
 for nm in v:
  for nn,zz,ee in ((2**256,1,1),(2**512,100000,100),(2**2048,1,2**22)):
   px(nm,nn,zz,ee)
 for nm in ('I','F','G','J','Fz','Gz','Jz','L0','Lr','E'):
  hh=len(v[nm]);px(nm,2**48,1,1,2**80);px(nm,2**256,256*prod(ar(nm,hh))*2**20,2**20-hh-522)
 UU=2**160;L=ceil(Q*UU);rf=UU;seed_s=2*ceil(((F(50,7)-F(1,2**100))*UU)/2);tag=3*(seed_s+34*rf)
 ee0=2*tag+12+18;om=omunit(ee0+L+37)
 nn=floor(do*om);zz=ceil(DX*om);adv=tag+2+lg(C);ee=ee0+4
 wn=max(2**80,ee+44+2*lg(4*nn));wf=max(2**80,ee+2*lg(4*nn)+32,ceil(F(zz,128)))
 m1,kf1=ceil(i['d']*wf),ceil(i['k']*wf);m0,k2=px('Ext2',nn,m1,ee)
 sig=ee+4;rr1=ceil((m0+PB*wn+DR+F(23,4)*ee+56-(6-eta/2)*x['c']*adv)/(theta-(6-eta/2)*x['c']))
 rr0=rr1;rr1=rr0-adv
 ss=floor((eta*rr0-11*sig-8)/4);rawv=6*rr1-floor(eta*rr1/2);d01,k01=px('Ext0',nn,rawv,ee)
 el=ceil(F(d01+ee-1,2))+ee;d11=m0+4*el+lg(4*(2*(rr1+adv+1)+1));d2=d01+d11+m1+k2+ee
 batch('NBexec.',[('raw-capacity',rr1,ss+sig,'>='),('raw-upper',rr1,rr*om,'<='),('seed-lower',m0,SD*wn,'>='),('I1',ss,m0+2*el,'>='),('advice',ss,adv+40,'>='),('density',F(adv+sig,rr1+adv),eta/DG,'<='),('AG1',ss,adv+40,'>='),('AG2',eta*rr0,DG*(adv+sig+2),'>='),('AG1-margin',ss,500*(adv+40),'>'),('AG2-margin',10*eta*rr0,F(1000000001,100000000)*DG*(adv+sig+2),'>'),('source',rr1-sig,(1-eta/DG)*(rr1+adv+1)+1,'>='),('joint-source',2*rr1-sig,(2-2*eta/DG)*(rr1+adv+1)+1,'>='),('joint-short',2*rr1-floor(eta*rr1/2)-sig,(2-eta/2-2*eta/DG)*(rr1+adv+1)+1,'>='),('capacity',F(5,4)*eta*(rr1+adv),5*ss+11*sig,'>='),('density-error',eta*(rr1+adv),4*ss+11*sig+8,'>='),('prime',rr1+adv,2**90,'>'),('dC',d01+d11+d2,floor(sc*om),'<='),('kC',max(k01,kf1+10*m0+3+zz+2*ee),floor(kc*om),'<=')])
 for ds0 in (2**160,2**161-1,2**200):
  dl0=ceil(Q*ds0);ds=2*ceil(((F(50,7)-F(1,2**100))*ds0)/2);tag1=3*(ds+34*ds0)
  e0=2*tag1+12+18;unit=omunit(e0+dl0+37);dim=floor(do*unit)
  for lu in (dl0+2,dl0+3):
   for nm,out in (('L0',floor(sc*unit)),('Lr',ceil(RATE*unit))):
    px(nm,None,out,e0,unit,cd=dim,lu0=lu)
 chi=max(F(j,(j-1))/(rr+1) for nm in ('L0','Lr') for rr,j in zip(v[nm],ar(nm,len(v[nm]))))
 srpad=slot(RR['pi'])
 for ssr in (2**160,2**161-1,2**200):
  lr0=ceil(QS*ssr);sr=2*ceil((F(50,7)-F(1,2**100))*ssr/2);er=6*sr+204*ssr+30;ur=er+lr0+37;rr=er+225+2*lg(ur)
  modes=[(floor(sc*ur),floor(s0*ur),floor(ko*ur),INIT_GAMMA0,INIT_XI0,INIT_G0,INIT_A0,INIT_B0),(ceil(RATE*ur),floor(srpad*ur),ceil(F(9,8)*ceil(RATE*ur)+7*alpha*ur),F(4),F(4,3),F(8,9),F(9,8),F(7))]
  for zr,padlim,srclim,gam,xx,ggg,aaa,bbb in modes:
   kr=ceil(aaa*zr+bbb*rr)
   for shift in (2,3):
    rr0=er+5+lr0+shift+lg(kr);b0r=ceil(F(N0+1,N0)*rr0);u0r=b0r-rr0;br=ceil(gam*rr);u1r=br-rr;vr0=b0r*ceil(F(kr,u0r));vr1=br*ceil(F(kr,u1r));cr=vr1-kr;vr2=br*ceil(F(cr,u1r));amb=2**(110+lg(ur))
    batch('AInit.',[('rho',rr,alpha*ur,'<='),('source',kr,srclim,'<='),('N0',rr0,N0,'>='),('v0',vr0,amb,'<='),('v1',vr1,amb,'<='),('positive',cr,1,'>='),('out',vr1-vr2,zr,'>='),('seed',b0r+2*br,padlim,'<='),('Kcap',kr,2**24*ur,'<'),('ceiling-gap',ggg*rr,xx,'>=')])
 batch('AInitScalar.',[('gamma',2<INIT_GAMMA0<4,True),('identity',INIT_G0*INIT_A0,1),('B-gap',INIT_G0*INIT_B0-INIT_XI0*INIT_GAMMA0,INIT_G0),('sourcecoef',max(INIT_A0*sc+INIT_B0*alpha,RR['lam']*(RATE+NE)+RR['g'])+NE,2**24,'<'),('rho-base',243+2*160,NE*2**160-14,'<'),('seed-base',2171+17*160,NE*246*2**160/4,'<'),('seed-growth',17,2171+17*160,'<='),('seed-ratio',F(2,N0)+NE/4,2*NE,'<'),('Lr-padding',DX,srpad,'>'),('round-gap',min(INIT_G0,F(8,9))*246*2**160,max(INIT_XI0,F(4,3)),'>')])
 batch('LC.s.',[('Kdim',ko,do,'<'),('Kcap',ko,2**24,'<'),('ratio',440,2*246,'<'),('bcap',chi*ko*(F(1728,7)+QS+F(1,10**6))/246,2**20,'<')])
 for nm in ('L0','Lr'):
  for rr,jj in zip(v[nm],ar(nm,len(v[nm]))):ck('LC.rate',rr,jj,'>=')
 # Exact guards for individually rounded outer ports and retained slot margins.
 batch('RoundPorts.',[('lambda',max(ZZ['lam'],RR['lam'],ev['E']['lam']),2**20,'<'),('margin',NE*246*2**160,2**31,'>'),('DY-gap',(DY-EE['c']*DX-EE['pi'])*OM,16,'>'),('DX-gap',(DX-RR['c']*RATE-RR['pi'])*OM,16,'>'),('R-gap',(RSHIFT-alpha)*OM,16,'>'),('gap-base',F(16,OM)*246*2**160,2**20,'>'),('cap-ceil',(ZCAP*(F(1728,7)-F(2,10**6))-DX*(F(1728,7)+QS+F(1,10**6)))*2**160,1,'>')])
 for ssr in (2**160,2**161-1,2**200):
  sr=2*ceil((F(50,7)-F(1,2**100))*ssr/2);tagr=3*(sr+34*ssr);er=2*tagr+12+18
  for lr0 in (ceil(QS*ssr),ceil(QS*(ssr+1))-1):
   ur=er+lr0+37;dyr=ceil(DY*ur);dxr=ceil(DX*ur);rrr=ceil(RATE*ur);scr=floor(sc*ur);kcr=floor(kc*ur);s0r=floor(s0*ur);dor=floor(do*ur);kor=floor(ko*ur)
   trr=(t+1)*(s0r+scr)+2*j*(t+1)*(dyr+dxr);tlr=(t+1)*scr+2*j*(t+1)*dyr
   ker=ceil(EE['lam']*dxr+EE['g']*ur);krr=ceil(RR['lam']*rrr+RR['g']*ur);k0r=ceil(ZZ['lam']*scr+ZZ['g']*ur)
   batch('RoundedCaller.',[('E-seed',dyr,EE['c']*dxr+EE['pi']*ur,'>='),('Lr-seed',dxr,floor(slot(RR['pi'])*ur),'>='),('Lm-seed',dxr,2*(lg(rrr*dyr+1)+dyr+er+1),'>='),('L0-seed',s0r,ZZ['c']*scr+ZZ['pi']*ur,'>='),('r3a',rrr,(t+2)*dyr+er+1,'>='),('d3c',dor,(t+1)*(s0r+scr)+kcr+er,'>='),('d3d',dor,trr+ker+t*dxr+er,'>='),('k3b',kor,tlr+krr+t*rrr,'>='),('K0',kor,k0r,'>='),('dimensions',0<dyr<=rrr<=kor<=dor,True)])
 batch('AffineRank.',[('init',5,5),('recurrence',3*5+8,23),('closed',9*3**9-4,177143),('error',177143,2**18,'<'),('AGHP-affine',F(1,2)+F(1,2),1,'<=')])
 UO=4*C*do+ko;e1=F(1728,7)-6*F(1,2**100);LA=e1/Q+1;AT=UO*LA;LOSS=4*AT;A=AT;CL=4*(UO*79+9);B=F(0);CC=(OM+113)*(UO+1)
 batch('refine.',[('source',A,do*LA/gg+Bbase/Q,'>'),('leading',A,H-1,'<='),('constant-guard',CC+1,2**160,'<'),('seed-margin',(3*ER-F(1,2**100))*2**160,165,'>='),('width-margin',(F(50,7)-F(1,2**100)-7)*2**160,162,'>='),('log-slack-positive',NE,0,'>'),('walk-base',F(7,16)*2**60,3320,'>'),('walk-slope',20,3320,'<='),('alpha-constant',55*alpha-47*(1-alpha),0,'>'),('alpha-floor',alpha*(246+Q),1,'>')])
 batch('UnitM.',[('offset',RSHIFT,alpha,'>='),('DY',OM*DY,14641446576009),('RATE',OM*RATE,4232468125065868),('DX',OM*DX,31463022350519),('Elog',lg(prod(JE))+6,48,'<='),('Ecap',DX*440,256*prod(JE)*246,'<=')])
 batch('OAcc.',[('positive',alpha,0,'>'),('cap',alpha,1,'<'),('constant',55*alpha-47*(1-alpha),18,'>'),('slope',eb+QS,1,'>'),('base',NE*2**160,262,'>'),('dyadic-slope',1,262,'<='),('logOm',440,512,'<'),('depth',max(len(JL0),len(JLR),len(JE)),5,'<='),('accuracy',max(lg(prod(jj)) for jj in (JL0,JLR,JE))+4,14,'<=')])
 ck('refine.round-slack',NE*Q*2**160,CC+2,'>');ck('result.H',H,216489736906);ck('result.K',K,865958947624)
 ck('result.raw',raw,216489736906);pt(P,f'H_raw=4C*delta_O+kappa_O+T={raw}');ck('result.base',bh,216489736494)
 ck('result.improvement',4*A,F('83450049915423073911796719301359574178203876738232717/18731372364015809820563343750000000000000'),'<');ck('result.tight_loss',[4*AT],[F('1435385536760996550382576784873316272232291153737375474466012346518307/2896718984356729864124679454581151760384000000000000000000')]);ck('result.loss',[4*A],[F('1435385536760996550382576784873316272232291153737375474466012346518307/2896718984356729864124679454581151760384000000000000000000')]);ck('refine.roundup',4*A,4*AT,'>=');ck('refine.ceil',4*A,4*AT+F(1,10**6),'<')
 if TEXT:
  for nm,vv in [('DX',DX),('DY',DY),('RATE',RATE),('gamma_0',INIT_GAMMA0)]:tc('ports.definition',nm+'='+str(vv) in TEXT[P],True)
  for ff in (P,W):
   for a,b in re.findall(r'(\d+)z[<>]=?(\d+)w',TEXT[ff]):tc('cap.guard',F(int(b),int(a)),ZCAP)
  pt(V,'Om=6s+204S+L+67');pt(P,'r>= (t+2)dy+e+1')
  for key,val in [('c1',4*A),('c0',CL)]:pt(P,key+'='+str(val))
  for nm,tag in TAG.items():pt(V,' v'+tag+'=('+VE[nm]+')')
  for f in (P,V,S):pt(f,'Q_S='+str(Q))
  for tok in ('TR=TR0+H(t+1)(2dy+2dx)-(t+1)dx','T=TR+TL+(t+1)dx'):pt(W,tok)
  for f,tok in [(P,f'dy=ceil(({DY})Om), dx=ceil(({DX})Om), r=ceil(({RATE})Om)'),(P,f's_C=floor(({sc})Om)'),(P,f'k_Cslot=floor(({kc})Om)'),(P,f'd_O=floor(({do})Om), k_O=floor(({ko})Om)'),(P,f'={tr}Om'),(P,f'={tl}Om'),(P,f'K=4H={K}'),(V,f'B_B=ceil((xi_B*gamma_B+xi_B/2^40)/G_B)={Bbase}')]:pt(f,tok)
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
  terms=bj['guarantee']['remainder']['terms'];declared=[F(t['coefficient']['numerator'])/F(t['coefficient']['denominator']) for t in terms];new=[4*A,CL]
  tc('meta.declared',declared==new)
  tc('meta.experiment',[e.get('source') for e in sj['experiments']],['verify.py'])
  tc('meta.evidence',all((ROOT/f).is_file() for claim in sj['claims'] for f in claim['evidence']))
  for claim in sj['claims']:tc('frozen.declaration',json.loads(claim['statement'].split('Declaration: ')[1]),bj)
 return dict(refined_loss=[str(4*A),str(CL)],tight_arithmetic_loss=[str(4*AT)],all_passed=not fail,all_checks_passed=not fail and not meta,checks=ct,arithmetic_checks=ac,text_checks_run=bool(TEXT),skipped_reason=None if TEXT else 'no package; arithmetic ran',H=str(H),K=str(K),H_raw=str(raw),H_base=str(bh),failures=fail,metadata_consistent=(not meta) if TEXT else None,metadata_failures=meta,groups=gs,scope='exact arithmetic+integrity')
HASHES={'sampler-and-disperser.md': '7b5f45defc10b48115efb6bee650ab10fa116563764221c904f467bfe416efb6', 'submission.json': '67da64d350acf789a1b9ff277b6907939529bda76c7c70d246d0c985c3e9af99', 'affine-and-standard-wrappers.md': '6a7cf80218da26927e8b81dc57b8e5608dc3f552d90c868df305dbe2c5b5d631', 'restricted-nonmalleability-and-breaking.md': 'c87fe8cb3035647b16ccde29add5e6e40d5d864b5a47cff70e8dfa83fc252185', 'proof.md': '1e50a1af79486e7b910d17d4006252318294485fde9688e9dcc7096e3fa90685', 'linear-and-bilinear-primitives.md': 'a70cfb657d1a43ed589d87095758ac1cafb6614cfe78102884478ab6d656e38f', 'bound.json': '22f3e3f2d2bab156719e68078889904dfda5d085ac0db192f999c9c328ea1e8e'}
if __name__=='__main__':
 if json.load(sys.stdin)!={}:raise ValueError('Expected {}')
 roots=(Path(__file__).resolve().parent,Path('submission'),Path('/snapshot/files/submission'))
 ROOT=next((p for p in roots if all((p/f).is_file() for f in (P,W,V,S,R,'bound.json','submission.json'))),None)
 if ROOT:TEXT={f:(ROOT/f).read_text() for f in (P,W,V,S,R)}
 res=run();out=json.dumps(res,sort_keys=True);assert len(out.encode())<65536
 print(out);raise SystemExit(0 if res['all_checks_passed'] else 1)
