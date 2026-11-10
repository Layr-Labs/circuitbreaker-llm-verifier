from hashlib import sha256
import importlib.util,sys
from pathlib import Path
path=Path('/snapshot/files/family.py')
if not path.is_file():path=Path(__file__).with_name('family.py')
assert sha256(path.read_bytes()).hexdigest()=='197762268ba0c8091e19561e27af0f92f6e84cfba9ea8a4a95c75dad0956f8b8'
spec=importlib.util.spec_from_file_location('candidate_family',path)
f=importlib.util.module_from_spec(spec);spec.loader.exec_module(f)
assert __import__('json').loads(sys.stdin.readline())=={}
def parameter_certificate():
 n=1<<16384;K=1<<1024;take=1<<160;mode=(2,5);L=f.ceil_log2(n+1)
 cp=f.constant_parameters(n,K,mode=mode);rp=f.response_parameters(n,K,take,mode=mode)
 assert cp is not None and rp is not None and rp[0]==take
 retuned=max(2,f.ceil_root((take*(L+1)**12)**5,2));assert retuned<K
 assert(retuned-1)**2<(take*(L+1)**12)**5<=retuned**2
 def entropy(c):
  m,s,d,h,_=c;assert(1<<s)>(1<<30)*h*d*s*(32*m)
  return 400*m*h*d*s*s+f.ceil_log2(32*m)+3
 td,tr=entropy(cp),entropy(rp);assert td<K and tr<=retuned and tr<(1<<700)
 sn,sd=f.regime(mode)['seed_guard']
 for c in(cp,rp):assert c[-1]*sd<sn*f.ceil_log2(n)
 return dict(n_log2=16384,K_log2=1024,take_log2=160,default_entropy_bits=td.bit_length(),response_entropy_bits=tr.bit_length(),retuned_bits=retuned.bit_length(),default_seed=cp[-1],response_seed=rp[-1],passed=True)
def small_certificate():
 n=1<<1470;r=n;mode=(2,5);L=f.ceil_log2(n+1)
 pp=f.small_parameters(n,r,mode=mode);assert pp is not None
 m,z,cs,ell,p,Q,A,d=pp;limit=2*(r.bit_length()-1)//10-24*f.ceil_log2(L+1)
 assert (limit,p,Q,z,cs,ell,A,d)==(30,29,536870912,247,11047,22094,47,288230376151711744)
 assert all(not f.is_prime(j) for j in range(p+1,limit+1))
 assert m==r//(L+1)**12 and Q==1<<p and z==Q//(L+1)**2 and d==Q*Q
 assert ell<=Q and A<=ell and 2*(d+m*(1<<(A-1))+cs+z+f.ceil_log2(4*m))<=r
 near=f.near_parameters(n,r,mode=mode);M,J,seed=near
 assert M==r//4 and J==(M+m-1)//m and seed==J*d and J*(1<<100)<(1<<z)
 return dict(n_log2=1470,limit=limit,p=p,z=z,code_prime=cs,design_seed=str(d),near_seed_bits=seed.bit_length(),near_error_below_2minus100=True,passed=True)
from fractions import Fraction as F
from math import lcm,isqrt
from pathlib import Path
import json,time
start=time.monotonic();a=F(6249,10000);w0=F(25419,10000);D=F(5593059,2000000);wm=D/2+F(1,10**6);X=D-w0;pc=a/2-F(1,10**10);t=F(4213,10000);c=F(4358,10000);cg=(1+c*t*t)/t
fb=(2*t**3/3+2*c*t**5/5)/cg;offset=fb-t**4/2;F1=F(1,2)+offset
# Alternating exact integral Taylor enclosure; actual expansion variable c*t²<.08.
def H(u,terms):return 2*cg*sum(((-c)**j*u**(2*j+3)/F(2*j+3)for j in range(terms)),F(0))
Hb_lo=H(t,16);Hb_hi=H(t,15);H1_hi=1-t*t+Hb_hi;assert Hb_lo<Hb_hi and Hb_hi-Hb_lo<F(1,10**18)
M=1048576;scale=10**15;uscale=10**12;step=(w0-wm)/M;V=lcm(wm.denominator,(step/2).denominator,D.denominator);U0=int((wm+step/2)*V);inc=int(step*V);dn=int(D*V);ap,aq=a.numerator,a.denominator;mid=0;right=[]
for j in range(M):
 U=U0+j*inc;z=U*aq-ap*V;num=scale*z*z;denf=aq*aq*U*V;q=num//denf;assert q*denf<=num<(q+1)*denf;mid+=q
 U+=inc//2;z=U*aq-ap*V;nr=scale*z*z;dr=aq*aq*U*V;right.append(nr//dr+1)
 if j in(0,1,M//2,M-1):
  for pos,value in [(U-inc//2,q),(U,right[-1]-1)]:
   w=F(pos,V);beta=a/w;alpha=1/(1-beta);exact=scale*w*(1-beta)/(1+alpha*beta);assert value==exact.numerator//exact.denominator
assert mid==969265479979700582795 and sum(right)==969265988972203121754
I=step*F(mid,scale);need=F1/(4*I);h=F(-(-need.numerator*10**14//need.denominator)+10000,10**14);suffix=sum(right);Jsum=0;branches=[0,0,0];yscale=uscale
# q=F1-4h*rightsuffix; fixed common rational denominator.
den=lcm(F1.denominator,(4*h*step/scale).denominator,fb.denominator,offset.denominator);Fnum=int(F1*den);coef=int(4*h*step/scale*den);fbnum=int(fb*den);offnum=int(offset*den)
# F(u²)= (2u³/3+2cu5/5)/cg. Compare using integers only.
cp,cq=c.numerator,c.denominator;gp,gq=cg.numerator,cg.denominator
# For floor u on 1/uscale, multiply left polynomial and target by common factors.
polyA=10*cq*uscale**2;polyB=6*cp;polyden=15*cq*uscale**5*gp;polymul=gq*den
hb_floor=(Hb_lo*yscale).numerator//(Hb_lo*yscale).denominator;boundary=int(t*uscale)
for j in range(M):
 rn=max(0,Fnum-coef*suffix)
 if not rn:y=0;branches[0]+=1
 elif rn>=fbnum:
  sq=isqrt(2*(rn-offnum)*yscale**2//den);assert sq*sq*den<=2*(rn-offnum)*yscale**2<(sq+1)**2*den
  y=sq-int(t*t*yscale)+hb_floor;branches[2]+=1
 else:
  lo,hi=0,boundary+1
  while lo+1<hi:
   u=(lo+hi)//2
   if u**3*(polyA+polyB*u*u)*polymul<=rn*polyden:lo=u
   else:hi=u
  assert lo**3*(polyA+polyB*lo*lo)*polymul<=rn*polyden<hi**3*(polyA+polyB*hi*hi)*polymul
  # H(u)=2Cg sum_j(-c)^j u^(2j+3)/(2j+3);16terms giveslower.
  y=0;upow=lo**3;cpow=1;qpow=1;spow=uscale**3
  for k in range(16):
   nr=2*gp*cpow*upow*yscale;dr=gq*qpow*spow*(2*k+3)
   q=nr//dr;assert q*dr<=nr<(q+1)*dr
   y+=q if k%2==0 else (-nr)//dr
   upow*=lo*lo;cpow*=cp;qpow*=cq;spow*=uscale*uscale
  y=max(0,y);branches[1]+=1
 U=U0+j*inc-inc//2;gap=dn-U;Jsum+=scale*y*V*V//(yscale*gap*gap);suffix-=right[j]
assert suffix==0 and branches==[2,84277,964297] and Jsum==2200333850470845559336
energy=4*h*I/F1;confidence=2*h*(H1_hi/X-step*F(Jsum,scale));row=dict(D=str(D),a=str(a),w0=str(w0),wm=str(wm),X=str(X),t=str(t),c=str(c),cg=str(cg),F1=str(F1),fb=str(fb),offset=str(offset),Hb_lo=str(Hb_lo),Hb_hi=str(Hb_hi),h=str(h),p=str(pc),M=M,mid=str(mid),right=str(sum(right)),Jsum=str(Jsum),branches=branches,energy=str(energy),confidence=str(confidence),margin=str(pc-confidence),passed=energy>1 and confidence<pc,elapsed=time.monotonic()-start)

assert row['passed'] and f.analytic_prefilter((h,a,wm,w0,X,pc)) and D<F(279653,100000)
W=w0-wm;mu=energy-1-8*h*W/10**12;assert mu>0;target=176*h*W*W/mu;leftmesh=1+-(-target.numerator//target.denominator)
assert energy-8*h*W/10**12-88*h*W*W/leftmesh>1+mu/2>1
numeric=dict(terms=M,MID=str(mid),RIGHT=str(sum(right)),Jsum=str(Jsum),branches=branches,left_mesh=str(leftmesh),margin=str(pc-confidence),passed=True)
print(json.dumps(dict(stage='numeric',result=numeric)),file=sys.stderr,flush=True)
su=parameter_certificate();print(json.dumps(dict(stage='su_response',result=su)),file=sys.stderr,flush=True)
small=small_certificate();print(json.dumps(dict(stage='small_near',result=small)),file=sys.stderr,flush=True)
L=10**12;B=f.ceil_log2(L+1);Lpad=10**6;ell_log=f.ceiling(h*f.ceil_root(L,2)*f.ceil_root(B,2));previous=ell_log
for v in (0,L//2,4*L//5,9*L//10,99*L//100,L-1):
 value=f.profile_log(h,L,v,Lpad,B,ell_log);u=max(L-v,Lpad)
 z=h*cg*L/(L+c*u);square=z*z*u*B;exact=f.ceil_root(f.ceiling(square),2)
 assert (exact-1)**2<square<=exact**2
 assert value==min(ell_log,exact) and 0<value<=previous;previous=value
print(json.dumps(dict(all_passed=True,numeric=numeric,su_response=su,small_near=small,profile_probe=True)))
