from fractions import Fraction as F
from math import ceil,floor,isqrt,comb,prod
from pathlib import Path
import json,sys,re,hashlib
EPS=F(1,2**20);NU=F(1,2**80); NE=F(1,10**18);AD=F(1,10**40); EB=F(1,10**18);ER=F(1,2**80);N0=2**60;AB=3;DG=F(4697,100);AP=F(1728,7)+F(1,10**6);AG=3*DG/2;ATG=AG+F(1,10**12)
SL={'a': '1025159208321000000500001/250000000000000000000000','s': '1653827960948000000500001/250000000000000000000000','s0': '2024891758030000000500001/500000000000000000000000','s0s': '2001433313178000000500001/500000000000000000000000','hc': '87035787326359839978663443/500000000000000000000000','h2': '14761032727755167994532687/100000000000000000000000','az': '2050318416642000000500001/500000000000000000000000','sz': '3307655921896000000500001/500000000000000000000000','d': '26552686387019180912931861/100000000000000000000000','k': '1855515966979300916395324671/500000000000000000000000'}
WR=[('11000', 4, 0, 14, 1), ('11001', 0, 4, 32, 2), ('11002', 0, 8, 37, 2), ('11003', 0, 8, 38, 3), ('10100', 4, 0, 14, 2), ('10101', 4, 4, 19, 1), ('10102', 4, 6, 33, 2), ('10103', 4, 6, 34, 3), ('10010', 4, 0, 14, 2), ('10011', 4, 4, 19, 2), ('10012', 4, 8, 24, 1), ('10013', 4, 6, 34, 3), ('01100', 6, 0, 3, 0), ('01101', 6, 2, 19, 1), ('01102', 6, 4, 33, 2), ('01103', 6, 4, 34, 3), ('01010', 6, 0, 3, 0), ('01011', 6, 2, 19, 2), ('01012', 6, 6, 24, 1), ('01013', 6, 4, 34, 3)]
VN={'E': '1000230392928287'}
AFF={'X': '5975545601/2500000000','I': '202489175803/100000000000','Is': '1000716656589/500000000000','F': '1025159208321/500000000000','Fz': '41/20','J': '1025159208321/500000000000','Jz': '1025159208321/500000000000','G': '27/10','Gz': '413456990237/125000000000','Ext2': '3544080630229/1000000000000'}
VE={k:','.join(str(F(int(x),10**15)) for x in s.split(',')) for k,s in VN.items()}
QS=F(213941337231821,100000000000000);QC=17;DX=F(5106975461447733700219,1180591620717411303424) 
RSHIFT=F(40779241983047950732003,1180591620717411303424)
ZCAP=F(1771192802333,62500000000)
GAM0=F(200035107332033,100000000000000)
XI0=GAM0/(GAM0-1);G0=1-(XI0-1)**2;A0=1/G0;B0=XI0*GAM0/G0
GAMR=F(218163313656157,100000000000000)
XIR=GAMR/(GAMR-1);GR=1-(XIR-1)**2;AR=1/GR;BR=XIR*GAMR/GR
M_GAM=F(218163313656157,100000000000000)
M_XI=M_GAM/(M_GAM-1);M_G=1-(M_XI-1)**2;M_A=1/M_G;M_B=M_XI*M_GAM/M_G
P='proof.md';W='affine-wrappers-and-sampler.md';V='linear-and-bilinear-primitives.md'
S='affine-wrappers-and-sampler.md';R='restricted-nonmalleability-and-breaking.md'
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
JS=(4,5,5,5,5,6);JY=(3,4);JI=(3,3,3);JG=(3,);JL0=(3,4,4,4,4);JLR=(3,3,4,4,4);JE=(2,);JFZ=(3,4);OM=2**70;BM=prod(JS);BY=prod(JY)
def ak(gamma):
 gamma=F(gamma);xi=gamma/(gamma-1);gg=1-(xi-1)**2
 return dict(c=F(0),pi=2*gamma+NE,lam=1/gg,g=xi*gamma/gg+NE)
def envm(vec,js=JE,acc=F(1)):
 assert js==(2,) and len(vec)==1
 r=vec[0];lam=2*r/(r-1)
 return dict(lam=lam/2,g=3*acc*lam+(r+1)/2,c=F(1,2),pi=4*acc+r+1)
def ar(nm,h):
 assert nm=='E'
 return JE
def omunit(x):return x
def slot(x):return F(ceil(x*10**18)+2,10**18)
def lm_cost(N,z,e,O,g):
 xi=g/(g-1);G=1-(xi-1)**2;rho=e+225+2*lg(O);K=ceil(z/G+xi*g/G*rho+xi/G)
 return dict(seed=2*ceil(g*rho),K=K)
def run():
 v={k:list(map(F,s.split(','))) for k,s in VE.items()};i={k:F(s) for k,s in SL.items()}
 acc=(F(1728,7)-6*F(1,2**79))/(F(1728,7)-6*F(1,2**79)+QS)+AD
 ev={k:envm(a,ar(k,len(a)),acc if k in ('L0','Lr') else F(1)) for k,a in v.items()}
 ev['E']['g']*=acc;ev['E']['pi']*=acc
 ev.update({nm:ak(AFF[tag]) for nm,tag in [('Ext0','X'),('I','I'),('Is','Is'),('F','F'),('Fz','Fz'),('J','J'),('Jz','Jz'),('G','G'),('Gz','Gz'),('Ext2','Ext2')]})
 I,FF,G,J=[ev[k] for k in ('I','F','G','J')];FZ,GZ,JZ=[ev[k] for k in ('Fz','Gz','Jz')];x,y=ev['Ext0'],ev['Ext2'];IS=ev['Is']
 SD=(1-F(1,2**40))*ev['Ext2']['pi'];DR=F(0);PB=x['pi']-DR*EPS
 h_R=F(3)+F(2,DG)+2*NU;DD=F(17)+F(8,DG)+5*NU;eta=(2-NU)/(14*DD+h_R+F(5,2));theta=eta/4;gamma=DD*eta;delta_R=gamma+NU;q_R=14*gamma-eta+NU;U=2**40
 eb=F(1728,7)-6*F(1,2**79);alpha=eb/(eb+QS)+AD
 df=i['d']*alpha+NE;kf=i['k']*alpha+NE
 M0=y['c']*df+y['pi']*alpha+NE;K2=y['lam']*df+y['g']*alpha+NE
 rr=(M0+PB*alpha+F(23,4)*alpha+(56+DR)*NE)/(theta-(6-eta/2)*x['c'])+NE
 amin=(F(864,7)-3*F(1,2**79))/(F(1728,7)-6*F(1,2**79)+QS);rr-=(1+(6-eta/2)*x['c']/(theta-(6-eta/2)*x['c']))*amin
 d0=(6-eta/2)*x['c']*rr+PB*alpha+(DR+2)*NE;d1=M0+2*d0+6*alpha+NE
 k00=(6-eta/2)*x['lam']*rr+x['g']*alpha+NE;kother=kf+2*(6-eta/2)*rr+(5*i['s']+i['sz'])*alpha+DX+2*alpha+100*NE
 dc=2*d0+2*d1+df+K2+alpha;sc=slot(dc);kc=slot(max(k00,kother))
 q=QC;C=4*(q+1);t=4*C-1;j=lg(t);Q=QS;beta=F(2*(q+8),7);p=144
 EE=ev['E'];RR=dict(c=F(0),pi=2*GAMR*alpha,lam=AR,g=BR*alpha);ZZ=dict(c=F(0),pi=1+2*GAM0*alpha,lam=A0,g=B0*alpha)
 DY=F(9576444264106363929335,1180591620717411303424);SM=slot(2*M_GAM*alpha);KM=slot(M_A*(DY+NE)+M_B*alpha);RATE=(t+2)*DY+RSHIFT
 s0=slot(ZZ['c']*sc+ZZ['pi']);tr=(t+1)*(s0+sc+(2*j-1)*(DY+DX));tl=(t+1)*(sc+2*(j-1)*DY)
 rows=dict(dy_E=EE['c']*DX+EE['pi'],dx_Lr=RR['c']*RATE+RR['pi'],dx_Lm=SM,s0_L0=ZZ['c']*sc+ZZ['pi'],r_3a=(t+1)*DY+KM,r_dy=DY,r_Km=KM,
 d_3c=(t+1)*(s0+sc)+kc+alpha,d_3d=tr+EE['lam']*DX+EE['g']+t*DX+alpha,d_s0=s0,d_KE=EE['lam']*DX+EE['g'],d_kC=kc,
 k_3b=tl+RR['lam']*RATE+RR['g']+t*RATE,k_K0=ZZ['lam']*sc+ZZ['g'])
 do=slot(max(val for key,val in rows.items() if key.startswith('d_')));ko=slot(max(val for key,val in rows.items() if key.startswith('k_')))
 ALIASES.update({str(DX):'DX',str(DY):'DY',str(RATE):'RATE'});ALIASES.update({str(do):'Dslot',str(ko):'Kslot',str(s0):'S0'})
 co=(3*p+1)/Q+1+F(1,10**15);delta=ceil(do*co);kappa=ceil(ko*co);T=1
 a0=F(N0+1,N0);gb=(beta-a0*(Q+1+EB)-3*ER)/(2*(1+EB));xi=gb/(gb-1);gg=1-(xi-1)**2;aa=1/gg;bb=ceil((xi*gb+xi/U)/gg)
 Bbase=bb;raw=4*C*delta+kappa+T;bh=ceil(aa*delta+bb*(1+EB)/Q+ER);H=max(raw,bh);K=4*H
 kb=aa*delta*Q*(1+F(1,U))+bb*(1+EB)+ER;b0c=a0*(Q+1+EB)+ER
 IR=lambda I,IS,FF,G,J,FZ,GZ,JZ:[('aJ','a',J['c']*(i['s']+NE)+J['pi']),('aF','a',FF['c']*(i['hc']+NE)+FF['pi']),('sG','s',G['c']*(i['a']+NE)+G['pi']),('s0I','s0',I['c']*(i['hc']+NE)+I['pi']),('s0Is','s0s',IS['pi']),('s0order','s0',i['s0s']+NE),('hJ','hc',(9+J['lam'])*(i['s']+NE)+J['g']+1),('hJ2','h2',(5+J['lam'])*(i['s']+NE)+J['g']+1),('horder','hc',i['h2']+NE),('hJz','h2',(5+JZ['lam'])*(i['sz']+NE)+JZ['g']+1),('hcap','h2',F(100)),('h>=s','h2',i['s']+NE),('azF','az',FZ['c']*ZCAP+FZ['pi']),('azJz','az',JZ['c']*(i['sz']+NE)+JZ['pi']),('szGz','sz',GZ['c']*(i['az']+NE)+GZ['pi']),('sz<=s','s',i['sz']+NE),('az<=a','a',i['az']+NE),('dEarly','d',6*(i['s0']+NE)+3*(i['s0s']+NE)+(40+G['lam'])*(i['a']+NE)+G['g']+1),('dLast','d',6*(i['s0']+NE)+3*(i['s0s']+NE)+42*(i['a']+NE)+(11+GZ['lam'])*(i['az']+NE)+GZ['g']+1),('d>=s0','d',i['s0']+NE),('kI','k',I['lam']*(i['hc']+NE)+I['g']),('kIs','k',IS['lam']*(i['s']+NE)+IS['g']),('k>=d','k',i['d']+NE)]+[(nm,'k',aa*(i['hc']+NE)+bb*(i['h2']+NE)+cc*(i['s']+NE)+(FF['lam']*([i['hc'],i['h2'],i['s']][tt]+NE)+FF['g']+1 if tt<3 else (8+FZ['lam'])*ZCAP+FZ['g']+1)) for nm,aa,bb,cc,tt in WR];ir=IR(I,IS,FF,G,J,FZ,GZ,JZ)
 for nm,ky,val in ir:
  ck('inner.'+nm,i[ky],val,'>=')
 for ky,val in rows.items():
  lim=dict(dy=DY,dx=DX,r=RATE,s0=s0,d=do,k=ko)[ky.split('_')[0]]
  ck('outer.'+ky,lim,(F(ceil(val*10**18),10**18)+NE if ky.split('_')[0] in ('s0','d','k') else val),'>=')
 den=theta-(6-eta/2)*x['c'];runit=(i['d']/BY+y['pi']+PB+6)/den
 batch('Zcap.',[('Q',QS,F(3),'<'),('dyadic',10**12*9*162,2**160,'<'),('E',F(1728,7)-6*F(1,2**79)-F(1,10**12),F(1728,7)-F(2,10**12),'>'),('Om',(OM+116)*10**12,2**160,'<'),('z',DX*(F(1728,7)+QS+F(1,10**12)),ZCAP*(F(1728,7)-F(2,10**12)),'<')])
 batch('NBshift.',[('beta',3*(F(50,7)-F(1,2**79))+102,F(12342857,100000),'>='),('capacity',eta/4+eta/DG,1,'<'),('amin',amin,0,'>')])
 batch('CLguard.',[('Om',co*Q*(1+F(1,2**80)),440,'<='),('wF',DX*440,128*246,'<='),('Fguard',(100*246-DX*440)*2**160,1200,'>'),('Ext2cap',i['d']+NU,256*BY,'<='),('AKambient',do+NE,2**24,'<'),('rround',56+den+F(1,BY),F(2**80,4),'<=')])
 for nm,vec in v.items():
  en=ev[nm];h=len(vec)
  for ridx,rate in enumerate(vec):
   ck('PT.rate',rate,ar(nm,len(vec))[ridx]-1,'>')
   sh=ceil(7*(rate+1));nn=ceil(F(3*(sh-14),2))-2
   bc=ceil(F(2*(nn+2),sh-14));mm=ceil(F(bc*sh,2));kk=mm+nn+1
   ck('PT.boundary',kk,bc*(sh-7),'<')
  kap=(vec[-1]+1)/(vec[-1]-ar(nm,len(vec))[-1]+1)
  batch('CLguard.',[('kap',kap,2**20,'<='),('last',vec[-1],64,'<='),('height',h,13,'<=')])
  if nm not in ('Ext0','Ext2'):
   out=dict(I=i['hc']+EPS,F=i['hc']+EPS,G=i['a']+EPS,J=i['s']+EPS,Fz=ZCAP,Gz=i['az']+EPS,Jz=i['sz']+EPS,L0=sc,Lr=RATE,E=DX)[nm]
   ck('CLguard.output',out,256*prod(ar(nm,h)),'<=')
 cl=[ceil(z) for z in (M0,K2,rr,d0,d1,dc)]
 eu={k:{c:F(ceil(e[c]*10**40),10**40) for c in e} for k,e in ev.items()}
 for c in ('c','pi','lam','g'):
  for k in ev:ck('AE.enclosure',0<=eu[k][c]-ev[k][c]<F(1,10**40),True)
 xu,yu=eu['Ext0'],eu['Ext2'];M0u=yu['c']*df+yu['pi']*alpha+NE;K2u=yu['lam']*df+yu['g']*alpha+NE
 ck('AE.denominator',theta,(6-eta/2)*xu['c'],'>');ru=(M0u+(xu['pi']-DR*EPS)*alpha+F(23,4)*alpha+(56+DR)*NE)/(theta-(6-eta/2)*xu['c'])+NE-(1+(6-eta/2)*xu['c']/(theta-(6-eta/2)*xu['c']))*amin
 d0u=(6-eta/2)*xu['c']*ru+(xu['pi']-DR*EPS)*alpha+(DR+2)*NE;d1u=M0u+2*d0u+6*alpha+NE;dcu=2*d0u+2*d1u+df+K2u+alpha
 aec=[ceil(z) for z in (M0u,K2u,ru,d0u,d1u,dcu,(6-eta/2)*xu['lam']*ru+xu['g']*alpha+NE,kf+2*(6-eta/2)*ru+(5*i['s']+i['sz'])*alpha+DX+2*alpha+100*NE)]
 ck('AE.seed',dcu,sc-NE,'<=');ck('AE.source',max((6-eta/2)*xu['lam']*ru+xu['g']*alpha+NE,kf+2*(6-eta/2)*ru+(5*i['s']+i['sz'])*alpha+DX+2*alpha+100*NE),kc-NE,'<=');
 Iu,Fu,Gu,Ju=[eu[k] for k in ('I','F','G','J')]
 for nm,ky,val in IR(Iu,eu['Is'],Fu,Gu,Ju,eu['Fz'],eu['Gz'],eu['Jz']):ck('AE.inner.'+nm,i[ky],val,'>=')
 Eu=eu['E'];Ru={k:F(ceil(v*10**40),10**40) for k,v in RR.items()};Zu={k:F(ceil(v*10**40),10**40) for k,v in ZZ.items()}
 for nm,lim,val in [('dyE',DY,Eu['c']*DX+Eu['pi']),('dxLr',DX,Ru['c']*RATE+Ru['pi']),('s0L0',s0,Zu['c']*sc+Zu['pi']),('d3d',do,tr+Eu['lam']*DX+Eu['g']+t*DX+alpha),('dKE',do,Eu['lam']*DX+Eu['g']),('k3b',ko,tl+Ru['lam']*RATE+Ru['g']+t*RATE),('kK0',ko,Zu['lam']*sc+Zu['g'])]:ck('AE.outer.'+nm,lim,val+(NE if nm in ('s0L0','d3d','dKE','k3b','kK0') else 0),'>=')
 for val in (k00,kother):ck('NB.source',val,kc,'<=')
 ck('NB.seed',dc,sc,'<=');ck('NB.source-dimension',kc,do,'<=');ck('NB.seed-dimension',sc,do,'<=')
 batch('NB.',[('denominator',theta,(6-eta/2)*x['c'],'>'),('seed-mass',PB+SD+F(23,4)*(1-F(1,2**40)),16,'>='),('cap',DX,128,'<='),('advice',F(23,4)+F(56,U),28,'<='),('unit-mass',4*(PB+SD+F(23,4)*(1-F(1,2**40))),64,'>='),('AG2-slope',AG,ATG,'<'),('AP-base',418,F(2**80,10**6),'<'),('AP-growth',2,418,'<='),('Ee-base',290,F(2**80,2**40),'<'),('Ee-growth',2,290,'<='),('AG2-const',(ATG-AG)*2**80,(6*DG),'>'),('AG2-r-lower',eta*(PB+SD+F(23,4)*(1-F(1,2**40))-(6-eta/2)*x['c']/2)/(theta-(6-eta/2)*x['c']),ATG,'>='),('AG2-seed-round',33,F(2**80,2**40),'<='),('seed-aa-base',441,F(2**160,2**40),'<'),('seed-aa-slope',2,441,'<='),('constant',ceil(2*(rr+2)*M0)+1,2**240,'<'),('weaklog',4*(6+1),(6+1)**2,'<='),('Lm',RATE*DY+1,2**33,'<'),('heights',lg(prod(JE))+13,21,'<='),('outer-heights',max(lg(prod(ar(k,len(v[k])))) for k in ('E',))+12,25,'<='),('NB-heights',max(lg(BM),lg(BY))+6,27,'<='),('prime',4*16*2**80/eta,2**90,'>'),('capacity',F(5,4)*4,5,'=='),('raw-blocks',1+1+(2-eta/2)+2,6-eta/2,'=='),('raw-round',max(4,x['lam']),NE*2**80,'<='),('c0',x['c'],1,'<'),('raw-rows',4+1,5,'=='),('matching-charge',2*4+2,10,'=='),('sig-capacity',F(5,4)*11,11,'>'),('floor',56-1-13-1,F(41),'>='),('advice40',F(41),40,'>'),('error',14,16,'<'),('final-error',15,16,'<')])
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
 batch('local.',[('alpha',alpha,1,'<'),('output',F(DX,128)+NE,alpha,'<='),('source',do,2**30,'<'),('Omega',co*Q*(1+F(1,2**80)),1024,'<'),('unit-base',520,AD*(eb+Q)*2**160,'<'),('unit-slope',2,520,'<='),('hash-base',804,NE*2**80,'<'),('hash-slope',4,804,'<='),('round',1,NE*2**80,'<='),('advice-charge',F(23,4),F(23,4)),('AG1-coef',9*F(1,3),1,'>='),('AG2-coef',DG*F(1,2)+DG,F((3*DG),2),'==')])
 for ss in sorted({s for k in range(160,201) for s in (2**k,2**(k+1)-1)}|{2**160+j for j in range(8)}):
  k=ss.bit_length()-1;s=2*ceil(((F(50,7)-F(1,2**79))*ss)/2);ellc=lg(s+4);rf=ss+0;ab=rf+AB;u1=2**(110+lg(ss));rho=ab+2*lg(u1)+2
  walk=s+34*rf;tag=3*walk
  batch('AFexec.',[('rho',rho,rf+225+2*lg(ss),'=='),('rho-dyadic',rho-ss,10*k+267,'<='),('rho-bound',rho,(1+EB)*ss,'<='),('log',ellc,k+5,'<='),('walk',F(7,16)*F(ss,2**79),100+20*(k+1),'>'),('parity',walk%2,0,'==')])
  for L in (ceil(Q*ss),ceil(Q*(ss+1))-1):
   ee=2*tag+12;om=omunit(ee+18+L-43);d=floor(do*om);kk=ceil(aa*d+bb*rho)
   eeN=ee+18+4;np=d;adv=tag+2+lg(C)
   wn=max(2**80,eeN+44+2*lg(4*np));rlo=ceil(((PB+SD)*wn+DR+F(23,4)*eeN+56-(6-eta/2)*x['c']*adv)/(theta-(6-eta/2)*x['c']));wf=max(2**80,eeN+2*lg(4*np)+32,ceil(F(DX*om,128)))
   batch('localexec.',[('Omega',om,1024*ss,'<='),('nlog',lg(4*np),42+lg(ss),'<='),('wN',wn,alpha*om,'<='),('wF',wf,alpha*om,'<='),('advice',2*adv,eeN,'<='),('advice-exact',2*adv,eeN-16,'=='),('AP',wn,AP*ss,'<='),('Ee-lower',eeN,(1-F(1,2**40))*wn,'>='),('AG2-upper',DG*(adv+eeN+6),AG*wn,'<='),('AG1-term',PB*wn+3*eeN+F(41),adv+40,'>='),('AG2-term',eta*(PB+SD+F(23,4)*(1-F(1,2**40))-(6-eta/2)*x['c']/2)*wn/(theta-(6-eta/2)*x['c']),DG*(adv+eeN+4+2),'>='),('AG2-wN',eta*(PB+SD+F(23,4)*(1-F(1,2**40))-(6-eta/2)*x['c']/2)*wn/(theta-(6-eta/2)*x['c']),AG*wn+(6*DG),'>'),('amin-exec',adv,amin*om,'>='),('raw-positive',rlo-adv,1,'>='),('raw-capacity',rlo-adv,floor((eta*rlo-11*(eeN+4)-8)/4)+eeN+4,'>='),('AG1-exec',floor((eta*rlo-11*(eeN+4)-8)/4),adv+40,'>='),('AG2-exec',eta*rlo,DG*(adv+eeN+4+2),'>='),('hash',480+4*lg(om),NE*om,'<=')])
   for shift in (2,3):
    rho0=ab+L+shift+lg(kk)+2;b0=ceil(a0*rho0);u0=b0-rho0;b=ceil(gb*rho);uu=b-rho
    v0=b0*ceil(F(kk,u0));v1=b*ceil(F(kk,uu));cc=v1-kk;v2=b*ceil(F(cc,uu));D=b0+2*b
    ee=2*tag+12;om=omunit(ee+18+L-43)
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
 ck('F.distinct',len(set(codes)),5)
 for aa in codes:
  first=aa.index('1');ck('F.prefix5',1+2*sum(bb[first]=='0' for bb in codes if bb!=aa),5,'<=')
  ck('F.weight',aa.count('1'),2)
  for bb in codes:
   if aa!=bb:ck('F.antichain',any(x=='1' and y=='0' for x,y in zip(aa,bb)))
 for a in codes:
  others=[x for x in codes if x!=a]
  for mask in range(81):
   cc=[a]+[x for j,x in enumerate(others) for _ in range(mask//3**j%3)];qv=[0,0,0]
   for qt in range(4):
    ss=[x for x in cc[1:] if any(a[j]=='1' and x[j]=='0' for j in range(qt))];nu=len(cc)-len(ss) if a[qt]=='1' else 0
    prev=len(cc)*qt+sum(x[:qt].count('1') for x in cc);cur=len(ss)+sum(x[qt]=='1' for x in ss)+nu
    ck('BiRight.ledger',prev+cur if qt<3 else prev+max(0,cur-11),40 if qt<3 else 42,'<=')
    ck('Widths.nu',nu,9 if '1' not in a[:qt] else 5,'<=')
    for qcopy in cc:
     if a[qt]=='1':qv[2]+=2 if qcopy[qt]=='1' else 1
     else:qv[2 if qcopy[qt]=='0' else int('1' in qcopy[:qt])]+=1
    ns=[x for x in cc[1:] if any(a[j]=='1' and x[j]=='0' for j in range(qt+1))];qq=qv[:]
    if qt<3:
     for qcopy in ns:qq[2 if qcopy[qt+1]=='0' else int('1' in qcopy[:qt+1])]+=1
    row=next(z for z in WR if z[0]==a+str(qt));ck('Widths.left',all(x<=y for x,y in zip(qq,row[1:4])) and (qt<3 or len(ns)<=8))
 for code in codes:
  b=code[:3].count('1');last=int(code[3]);ck('Widths.flat',last==0 and b==2 or last==1 and b==1)
 for d in range(1,20):
  for nn in range(1,40):
   vv=2*nn//(d+2);ck('SP.dimension-count',(vv+1)*comb(d+2,2),nn*(d+1),'>')
 for d0 in range(1,12):
  for d1 in range(d0):
   TT=d0+d1+7
   ck('SP.codimension',comb(TT+2,2)-comb(TT-d0+2,2)-comb(TT-d1+2,2)+comb(TT-d0-d1+2,2),d0*d1)
 def px(nm,nn,zz,ee,unit=0):
  assert nm=='E';rr=v[nm][0];a=ee+5;rh=a+2*lg(4*nn)+2;u=max(unit,2**20,rh);ell=ceil(F(zz,2));S=ceil((rr+1)*rh);K1=ell+2*a;b=ceil(F(2*(K1+a+1),S-2*rh));m=ceil(F(b*S,2));K0=m+K1+a;d=S+ell+4*a+lg(4*(m+1));en=envm(v[nm],JE)
  batch('PX.',[('capacity',K0,b*(S-rh),'<'),('count',b,K0,'<='),('child',m,K1,'>'),('a-unit',a,u-6,'<='),('expansion',4*(4*nn)**2,2**(rh-a),'<='),('source',K0,en['lam']*zz+en['g']*u,'<='),('padding',K0,nn,'<=')])
  if zz<=512*u:batch('PX.',[('seed',d,en['c']*zz+en['pi']*u,'<='),('fieldcap',m,2**29*u,'<='),('logcap',lg(4*(m+1)),32+lg(u),'<=')])
  return d,K0
 for nm in v:
  for nn,zz,ee in ((2**256,1,1),(2**512,100000,100),(2**2048,1,2**22)):
   px(nm,nn,zz,ee)
 for nm in ('E',):
  hh=len(v[nm]);px(nm,2**48,1,1,2**80);px(nm,2**256,256*prod(ar(nm,hh))*2**20,2**20-hh-522)
 def akc(role,N,z,e,O,w):
  gam=F(AFF[role]);xi=gam/(gam-1);gg=1-(xi-1)**2;aa=1/gg;bb=xi*gam/gg;rh=e+225+2*lg(O);kk=ceil(aa*z+bb*rh+xi/gg);b=ceil(gam*rh);u=b-rh;v1=b*ceil(F(kk,u));c=v1-kk;v2=b*ceil(F(c,u));amb=2**(110+lg(O));pad=ceil((2*gam+NE)*w)
  batch('AK.',[('gamma',2<gam<=16,True),('Bcap',bb,2**12,'<='),('xiGcap',xi/gg,2**12,'<='),('unit',rh,w+200,'<='),('reserve',NE*w,2**20,'>='),('source',kk,N,'<='),('ambient',N,amb,'<='),('field',rh,2**40,'>='),('v1',v1,amb,'<='),('codim',c,1,'>='),('output',v1-v2,z,'>='),('seed',2*b,pad,'<='),('source-envelope',kk,aa*z+(bb+NE)*w,'<='),('rho-error',rh,e+2*lg(amb)+5,'=='),('kernel-margin',gg*kk-xi*b,z,'>=')])
  return pad,kk
 UU=2**160;L=ceil(Q*UU);rf=UU;seed_s=2*ceil(((F(50,7)-F(1,2**79))*UU)/2);tag=3*(seed_s+34*rf)
 ee0=2*tag+12+18;om=omunit(ee0+L-43)
 nn=floor(do*om);zz=ceil(DX*om);adv=tag+2+lg(C);ee=ee0+4
 wn=max(2**80,ee+44+2*lg(4*nn));wf=max(2**80,ee+2*lg(4*nn)+32,ceil(F(zz,128)))
 m1,kf1=ceil(i['d']*wf),ceil(i['k']*wf);m0,k2=akc('Ext2',nn,m1,ee,om,wn)
 sig=ee+4;rr1=ceil((m0+PB*wn+DR+F(23,4)*ee+56-(6-eta/2)*x['c']*adv)/(theta-(6-eta/2)*x['c']))
 rr0=rr1;rr1=rr0-adv
 ss=floor((eta*rr0-11*sig-8)/4);rawv=6*rr1-floor(eta*rr1/2);d01,k01=akc('X',nn,rawv,ee,om,wn)
 el=ceil(F(d01+ee-1,2))+ee;d11=m0+4*el+lg(4*(2*(rr1+adv+1)+1));d2=d01+d11+m1+k2+ee
 batch('NBexec.',[('raw-capacity',rr1,ss+sig,'>='),('raw-upper',rr1,rr*om,'<='),('seed-lower',m0,SD*wn,'>='),('I1',ss,m0+2*el,'>='),('advice',ss,adv+40,'>='),('density',F(adv+sig,rr1+adv),eta/DG,'<='),('AG1',ss,adv+40,'>='),('AG2',eta*rr0,DG*(adv+sig+2),'>='),('AG2-margin',10*eta*rr0,F(1000000000001,100000000000)*DG*(adv+sig+2),'>'),('source',rr1-sig,(1-eta/DG)*(rr1+adv+1)+1,'>='),('joint-source',2*rr1-sig,(2-2*eta/DG)*(rr1+adv+1)+1,'>='),('joint-short',2*rr1-floor(eta*rr1/2)-sig,(2-eta/2-2*eta/DG)*(rr1+adv+1)+1,'>='),('capacity',F(5,4)*eta*(rr1+adv),5*ss+11*sig,'>='),('density-error',eta*(rr1+adv),4*ss+11*sig+8,'>='),('prime',rr1+adv,2**90,'>'),('dC',d01+d11+d2,floor(sc*om),'<='),('kC',max(k01,kf1+2*rawv+5*ceil(i['s']*wf)+ceil(i['sz']*wf)+zz+2*ee+100),floor(kc*om),'<=')])
 h2f=ceil(i['h2']*wf);hf,sf,af,s0f,azf,szf=[ceil(i[k]*wf) for k in ('hc','s','a','s0','az','sz')];qf=ee+7
 isp,ik=akc('I',nn,hf,qf,om,wf);fsp,fk=akc('F',nn,hf,qf,om,wf);f2sp,f2k=akc('F',nn,h2f,qf,om,wf);fssp,fsk=akc('F',nn,sf,qf,om,wf);issp,isk=akc('Is',nn,sf,qf,om,wf);zsp,zk=akc('Fz',nn,zz,qf,om,wf)
 gsp,gk=akc('G',m1,af,qf,om,wf);jsp,jk=akc('J',hf,sf,qf,om,wf);j2sp,j2k=akc('J',h2f,sf,qf,om,wf);gzsp,gzk=akc('Gz',m1,azf,qf,om,wf);jzsp,jzk=akc('Jz',h2f,szf,qf,om,wf)
 batch('AffineInner.',[('Iseed',isp,s0f,'<='),('Fseed',fsp,af,'<='),('F2seed',f2sp,af,'<='),('Fsseed',fssp,af,'<='),('Isseed',issp,ceil(i['s0s']*wf),'<='),('Fzseed',zsp,azf,'<='),('Gseed',gsp,sf,'<='),('Jseed',jsp,af,'<='),('Gzseed',gzsp,szf,'<='),('Jzseed',jzsp,azf,'<='),('source',kf1,max(ik,isk,*[aa*hf+bb*h2f+cc*sf+([fk,f2k,fsk][tt] if tt<3 else 8*zz+zk)+qf for _,aa,bb,cc,tt in WR]),'>='),('seed',m1,max(6*s0f+3*ceil(i['s0s']*wf)+40*af+gk+qf,6*s0f+3*ceil(i['s0s']*wf)+42*af+11*azf+gzk+qf),'>='),('h',hf,9*sf+jk+qf,'>='),('h2',h2f,5*sf+j2k+qf,'>='),('hz',h2f,5*szf+jzk+qf,'>='),('horder',hf,h2f,'>='),('h2s',h2f,sf,'>='),('h2z',h2f,zz,'>='),('flat-summary',5*sf+szf,(5*i['s']+i['sz'])*alpha*om+8,'<='),('mass',m0+PB*wn+F(23,4)*ee,16*wn,'>=')])
 srpad=slot(RR['pi'])
 for ssr in (2**160,2**161-1,2**200):
  lr0=ceil(QS*ssr);sr=2*ceil((F(50,7)-F(1,2**79))*ssr/2);er=6*sr+204*ssr+30;ur=er+lr0-43;rr=er+225+2*lg(ur)
  modes=[(floor(sc*ur),floor(s0*ur),floor(ko*ur),GAM0,XI0,G0,A0,B0),(ceil(RATE*ur),floor(srpad*ur),ceil(AR*ceil(RATE*ur)+BR*alpha*ur),GAMR,XIR,GR,AR,BR),(ceil(DY*ur),floor(SM*ur),floor(KM*ur),M_GAM,M_XI,M_G,M_A,M_B)]
  for mi,(zr,padlim,srclim,gam,xx,ggg,aaa,bbb) in enumerate(modes):
   kr=ceil(aaa*zr+bbb*rr+xx/ggg)
   for shift in ((2,) if mi==2 else (2,3)):
    al=lg(4*ceil(RATE*ur)) if mi==2 else lr0+shift
    rr0=er+5+al+lg(floor(ko*ur) if mi<2 else kr);b0r=ceil(F(N0+1,N0)*rr0);u0r=b0r-rr0;br=ceil(gam*rr);u1r=br-rr;vr0=b0r*ceil(F(kr,u0r));vr1=br*ceil(F(kr,u1r));cr=vr1-kr;vr2=br*ceil(F(cr,u1r));amb=2**(110+lg(ur))
    if mi==0:vr0=b0r*ceil(F(floor(ko*ur),u0r))
    if mi==1:b0r=0;u0r=1;vr0=2**(110+lg(ur))
    if mi==2:b0r=0;u0r=1;vr0=ceil(RATE*ur)
    batch('AInit.',[('rho',rr,alpha*ur,'<='),('source',kr,srclim,'<='),('N0',rr0,N0,'>='),('v0',vr0,amb,'<='),('v1',vr1,amb,'<='),('positive',cr,1,'>='),('out',vr1-vr2,zr,'>='),('seed',b0r+2*br,padlim,'<='),('Kcap',kr,2**24*ur,'<'),('rho0cap',rr0,2*ur,'<='),('ceiling-gap',ggg*rr,xx,'>=')])
 ck('AInitScalar.preK',ko+NE,2**24,'<');ck('AInitScalar.Mambient',RATE+NE,2**12,'<')
 batch('AInitScalar.',[('gamma0',2<GAM0<4,True),('gammaR',2<GAMR<16,True),('gammaM',2<M_GAM<16,True),('Midentity',M_G*M_A,1),('MB-gap',M_G*M_B-M_XI*M_GAM,0),('Ridentity',GR*AR,1),('RB-gap',GR*BR-XIR*GAMR,0),('identity',G0*A0,1),('B-gap',G0*B0-XI0*GAM0,0),('sourcecoef',max(A0*sc+B0*alpha,RR['lam']*(RATE+NE)+RR['g'],M_A*(DY+NE)+M_B*alpha)+NE,2**24,'<'),('rho-base',313+2*160,AD*2**160,'<'),('seed-base',8687+65*160,NE*246*2**160/4,'<'),('seed-growth',65,8687+65*160,'<='),('seed-ratio',F(2,N0)+NE/4,2*NE,'<'),('Lr-padding',DX,srpad,'>'),('Lm-padding',DX,SM,'>'),('Lm-input',RATE,KM,'>'),('Lm-Km',KM,M_A*(DY+NE)+M_B*alpha,'>'),('fixed-cap',max(XI0/G0,XIR/GR,M_XI/M_G),2**12,'<'),('fixed-floor',NE*246*2**160/4,2**12+2,'>')])
 batch('RoundPorts.',[('lambda',max(ZZ['lam'],RR['lam'],M_A,ev['E']['lam']),2**20,'<'),('margin',NE*246*2**160,2**31,'>'),('DY-gap',(DY-EE['c']*DX-EE['pi'])*OM,16,'>'),('DX-gap',(DX-RR['c']*RATE-RR['pi'])*OM,16,'>'),('R-gap',(RATE-(t+1)*DY-KM)*OM,16,'>'),('M-seed',(DX-SM)*OM,16,'>'),('gap-base',F(16,OM)*246*2**160,2**20,'>'),('cap-ceil',(ZCAP*(F(1728,7)-F(2,10**12))-DX*(F(1728,7)+QS+F(1,10**12)))*2**160,1,'>')])
 for ssr in (2**160,2**161-1,2**200):
  sr=2*ceil((F(50,7)-F(1,2**79))*ssr/2);tagr=3*(sr+34*ssr);er=2*tagr+12+18
  for lr0 in (ceil(QS*ssr),ceil(QS*(ssr+1))-1):
   ur=er+lr0-43;dyr=ceil(DY*ur);dxr=ceil(DX*ur);rrr=ceil(RATE*ur);scr=floor(sc*ur);kcr=floor(kc*ur);s0r=floor(s0*ur);dor=floor(do*ur);kor=floor(ko*ur)
   trr=(t+1)*(s0r+scr)+(2*j-1)*(t+1)*(dyr+dxr);tlr=(t+1)*scr+2*(j-1)*(t+1)*dyr
   ker=ceil(EE['lam']*dxr+EE['g']*ur);krr=ceil(RR['lam']*rrr+RR['g']*ur+XIR/GR);k0r=ceil(ZZ['lam']*scr+ZZ['g']*ur+XI0/G0)
   batch('RoundedCaller.',[('E-seed',dyr,EE['c']*dxr+EE['pi']*ur,'>='),('Lr-seed',dxr,floor(slot(RR['pi'])*ur),'>='),('Lm-seed',dxr,lm_cost(rrr,dyr,er,ur,M_GAM)['seed'],'>='),('Lm-source',floor(KM*ur),lm_cost(rrr,dyr,er,ur,M_GAM)['K'],'>='),('L0-seed',s0r,ZZ['c']*scr+ZZ['pi']*ur,'>='),('r3a',rrr,(t+1)*dyr+floor(KM*ur),'>='),('d3c',dor,(t+1)*(s0r+scr)+kcr+er+1,'>='),('d3d',dor,trr+ker+t*dxr+er+1,'>='),('k3b',kor,tlr+krr+t*rrr,'>='),('K0',kor,k0r,'>='),('dimensions',0<dyr<=rrr<=kor and s0r<=dor and kcr<=dor and scr<=dor,True)])
 batch('AffineRank.',[('init',5,5),('recurrence',3*5+8,23),('closed',3*(9*3**8-4)+4,177139),('error',F(708557,4),2**18,'<'),('AGHP-affine',F(1,2)+F(1,2),1,'<=')])
 UO=4*C*do+ko;e1=F(1728,7)-6*F(1,2**79);LA=e1/Q+1;AT=UO*LA;LOSS=4*AT;A=AT;CL=F(0);B=F(0);CC=(OM+113)*(UO+1)
 batch('refine.',[('source',A,do*LA/gg+Bbase/Q,'>'),('leading',A,H-1,'<='),('constant-guard',CC+1,2**160,'<'),('seed-margin',(3*ER-F(1,2**79))*2**160,165,'>='),('width-margin',(F(50,7)-F(1,2**79)-7)*2**160,162,'>='),('log-slack-positive',NE,0,'>'),('walk-base',F(7,16)*2**60,3320,'>'),('raw-intercept',UO,9,'>'),('walk-slope',20,3320,'<='),('alpha-constant',-25*alpha-47*(1-alpha),-47,'>'),('alpha-floor',alpha*246,1,'>')])
 batch('UnitM.',[('offset',RSHIFT,alpha,'>='),('DY',OM*DY,9576444264106363929335),('RATE',OM*RATE,2808371634309787126309818),('DX',OM*DX,5106975461447733700219),('Elog',lg(prod(JE))+6,48,'<='),('Ecap',DX*440,256*prod(JE)*246,'<=')])
 batch('OAcc.',[('positive',alpha,0,'>'),('cap',alpha,1,'<'),('constant',-25*alpha-47*(1-alpha),-47,'>'),('slope',eb+QS,1,'>'),('base',AD*2**160,633,'>'),('dyadic-slope',2,633,'<='),('logOm',440,512,'<'),('depth',max(len(JL0),len(JLR),len(JE)),5,'<='),('accuracy',max(lg(prod(jj)) for jj in (JL0,JLR,JE))+4,14,'<=')])
 ck('refine.round-slack',NE*Q*2**160,CC+2,'>');ck('result.H',H,17247453189);ck('result.K',K,68989812756)
 ck('result.raw',raw,17247453189);ck('result.base',bh,17247452783)
 ck('result.improvement',4*A,F('136418604801305797502302825932641527130897694686838312269953483/3453205918805690895515985951850496000000000000000000'),'<');ck('refine.roundup',4*A,4*AT,'>=');ck('refine.ceil',4*A,4*AT+F(1,10**6),'<')
 if TEXT:
  pt(P,f'Q_S={QS},DG={DG}')
  tc('text.ports',f'Ports: DX={DX},DY={DY},RATE={RATE}.' in TEXT[P])
  tc('text.outer-slots',f'Dslot={do},Kslot={ko}' in TEXT[P])
  pt(P,f'SC=Slot(DNB)={sc},KC=Slot(max(K00bar,Kotherbar))={kc}')
  pt(P,'c1='+str(4*A))
  pt(P,'delta_O=ceil(Dslot*c_O)='+str(delta));pt(V,'delta_O='+str(delta)+',')
  pt(P,'4(p0+SD+(23/4)f)-3DG/2-10^-12='+str(4*(PB+SD+F(23,4)*(1-F(1,2**40)))-3*DG/2-F(1,10**12))+'>0.')
  pt(V,f'B_B=ceil((xi_B*gamma_B+xi_B/2^40)/G_B)={Bbase}')
  pt(W,'eta*(r+a)>=4*(m0+p0*w_N+(23/4)Ee)>=64*w_N>=64*2^80')
  for role,tag in [('Ext0','X'),('I','I'),('Is','Is'),('F','F'),('Fz','Fz'),('J','J'),('Jz','Jz'),('G','G'),('Gz','Gz'),('Ext2','Ext2')]:
   en=ev[role];pt(P,f'AK{tag}: gamma={AFF[tag]}.')
  for nm,rs in VE.items():pt(P,f'{nm}: {rs}.')
  for key,val in i.items():pt(P,f'{key}=ceil(({val})w).')
  for f,tok in [(V,'## 7f-AK.'),(V,"rho=e'+225+2lg Om"),(V,'pi=2gamma+eps_N,Lam=A,g=B+eps_N'),(R,'Affine real-history invariant.'),(R,'the new WHOLE state tuple is uniform affine'),(R,'rank_aux(W0|H2)>=h-nu*s>=KJ+q'),(R,'## 8.5 Affine fibers and generated auxiliary records'),(W,"Ext2's affine AK guards are paid."),(W,"Ext0's requested AK output is v"),(W,'J_l=AK_J(inputl,outputs,seeda)'),(W,'Jz=AK_Jz(inputh2,outputsz,seedaz)'),(W,'rank loss at most2v'),(W,'Call this resulting history H5.'),(W,'5s_F+sz_F'),(W,'m0+p0*w_N+(23/4)Ee>=16w_N'),(W,'G=AK_G(inputd,outputa,seeds)'),(W,'Gz=AK_Gz(inputd,outputaz,seedsz)'),(R,'## 8.6 Selected two-sided-affine variant'),(W,'I_h1=AK_I(inputn,outputh1,seeds0),I_s=AK_Is(inputn,outputs,seeds0s)'),(W,'F_l=AK_F(inputn,outputl,seeda)'),(P,'The original tuple(B,Y) is affine in an independent uniform affine U')]:pt(f,tok)
 if TEXT:
  for nm,vv in [('DX',DX),('DY',DY),('RATE',RATE),('gamma_0',GAM0),('gamma_R',GAMR),('gamma_M',M_GAM)]:tc('ports.definition',nm+'='+str(vv) in TEXT[P],True)
  for ff in (P,W):
   for a,b in re.findall(r'(\d+)z[<>]=?(\d+)w',TEXT[ff]):tc('cap.guard',F(int(b),int(a)),ZCAP)
  tc('selection.Lm',re.search(r'Selected\s+Lm\s+uses\s*AGHPaff',TEXT[V]) is None)
  sp=TEXT[S].split('# SP. Finite rings')[1]+TEXT[V].split('# Relocated sampler appendices')[1]
  heads=re.findall(r'^#+ (SP\.[0-9.]+)',sp,re.M);tc('refs.headings',len(heads),len(set(heads)))
  labs=re.findall(r'^(?!#) *.*\((SP-[\w]+)\)\s*$',sp,re.M)
  tc('refs.labels',set(labs),set('SP-1 SP-2 SP-3 SP-4 SP-5 SP-5s SP-Mi0 SP-Mi1 SP-Mi2 SP-Mi3 SP-Mi4 SP-17 SP-18 SP-19m SP-21m SP-S2 SP-M SP-I'.split()))
  tc('refs.unique',len(labs),len(set(labs)))
  for lab in re.findall(r'\((SP-[\w]+)\)',sp):tc('refs.'+lab,lab in labs)
  for f in TEXT:
   tc('refs.removed','sum-product-and-ordinary-condensation.md' not in TEXT[f])
   for sec in re.findall(r'Section (SP\.\d+(?:\.\d+)?)',TEXT[f]):tc('refs.section',bool(re.search(r'^#+ '+re.escape(sec)+r'(?:[. ]|$)',sp,re.M)))
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
 return dict(refined_loss=[str(4*A)],tight_arithmetic_loss=[str(4*AT)],all_passed=not fail,all_checks_passed=not fail and not meta,checks=ct,arithmetic_checks=ac,text_checks_run=bool(TEXT),skipped_reason=None if TEXT else 'no package; arithmetic ran',H=str(H),K=str(K),H_raw=str(raw),H_base=str(bh),failures=fail,metadata_consistent=(not meta) if TEXT else None,metadata_failures=meta,groups=gs,scope='exact arithmetic+integrity')
HASHES={'linear-and-bilinear-primitives.md':'f35a49d3696ed7dfa9136994b4cee8aed803b1046501c155ec7eb2d9eecaa7eb','bound.json':'3939b1ab3bbe3928677034c168559f56357d3094bc442ba699dd44967c8cdee4','proof.md':'ebc1b3327af70d5ea979d24b20e48467fe825cdd794dae15e1ac4c9872bc176c','restricted-nonmalleability-and-breaking.md':'b864fd7143cfbc03a6b3aa94545e74ebdfe316283ad1c50d8b8b8458d6a6680c','affine-wrappers-and-sampler.md':'a900b0639f59d4ce31909b914adfc54f2274f7927661145d5f6ac36cd9092e0d','submission.json':'7742dfa73e151a20a345b1f8e526139e4eac8464ad5cb3ccce76f35f496f1ebc'}
if __name__=='__main__':
 if json.load(sys.stdin)!={}:raise ValueError('Expected {}')
 roots=(Path(__file__).resolve().parent,Path('submission'),Path('/snapshot/files/submission'))
 ROOT=next((p for p in roots if all((p/f).is_file() for f in (P,W,V,S,R,'bound.json','submission.json'))),None)
 if ROOT:TEXT={f:(ROOT/f).read_text() for f in (P,W,V,S,R)}
 res=run();out=json.dumps(res,sort_keys=True);assert len(out.encode())<65536
 print(out);raise SystemExit(0 if res['all_checks_passed'] else 1)
