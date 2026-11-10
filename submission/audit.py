from fractions import Fraction as F
from math import lcm
from pathlib import Path
from hashlib import sha256
import json,importlib.util,sys
path=Path('/snapshot/files/family.py')
if not path.is_file():path=Path(__file__).with_name('family.py')
assert sha256(path.read_bytes()).hexdigest()=='5d8a9fed61872e8eb583352ff178ec9aa3a087927bb5bde82a8e4825e9fe5326'
spec=importlib.util.spec_from_file_location('candidate_family',path)
f=importlib.util.module_from_spec(spec);spec.loader.exec_module(f)
def midpoint_certificate():
 from math import isqrt
 a=F(6249,10000);D=F(559309,200000);wm=D/2+F(1,10**6);w0=F(25419,10000);X=D-w0;h=F(11933947745977,10**14)
 pc=a/2-F(1,10**10);assert f.analytic_prefilter((h,a,wm,w0,X,pc))and D==w0+X<F(55931,20000)
 M=1048576;scale=10**15;yscale=18*10**12;step=(w0-wm)/M
 V=lcm(wm.denominator,(step/2).denominator,D.denominator);U0=int((wm+step/2)*V);inc=int(step*V);dn=int(D*V);ap,aq=a.numerator,a.denominator
 assert inc%2==0 and F(U0,V)==wm+step/2 and F(inc,V)==step
 total=0;right=[];cross=0
 for j in range(M):
  U=U0+j*inc;z=U*aq-ap*V;num=scale*z*z;den=aq*aq*U*V;q=num//den;assert q*den<=num<(q+1)*den;total+=q
  U+=inc//2;z=U*aq-ap*V;nr=scale*z*z;dr=aq*aq*U*V;upper=nr//dr+1;assert(upper-1)*dr<=nr<upper*dr;right.append(upper)
  if j in(0,1,M//2,M-1):
   for pos,value in[(U-inc//2,q),(U,upper-1)]:
    w=F(pos,V);beta=a/w;alpha=1/(1-beta);exact=scale*w*(1-beta)/(1+alpha*beta);assert value==exact.numerator//exact.denominator
   cross+=1
 assert total==969269008899573753702 and sum(right)==969269517888975238659
 suffix=sum(right);A=216*h.denominator*step.denominator*scale;B=864*h.numerator*step.numerator;Jsum=0;branches=[0,0,0]
 for j in range(M):
  rn=max(0,109*h.denominator*step.denominator*scale-B*suffix)
  if rn==0:y=0;branches[0]+=1
  elif 54*rn<=A:
   nr=6*rn*yscale;y=nr//A;assert y*A<=nr<(y+1)*A;branches[1]+=1
  else:
   nr=(216*rn-A)*yscale**2;dr=108*A;sq=isqrt(nr//dr);assert sq*sq*dr<=nr<(sq+1)**2*dr;y=sq-yscale//18;branches[2]+=1
  assert y>=0
  U=U0+j*inc-inc//2;gap=dn-U;num=scale*y*V*V;den=yscale*gap*gap;assert den>0;q=num//den;assert q*den<=num<(q+1)*den;Jsum+=q;suffix-=right[j]
 assert branches==[2,77041,971533]
 assert suffix==0 and Jsum==2200334831055428006970
 ratio=F(864,109)*h*step*F(total,scale);confidence=2*h*(F(17,18)/X-step*F(Jsum,scale))
 assert ratio==F(357171200299305012948303173302838822157837,357171200000000000000000000000000000000000)>1
 assert confidence==F(1501713933296693195895002093931429573818924323,4806263439360000000000000000000000000000000000)<pc
 W=w0-wm;mu=ratio-1-8*h*W/10**12;assert mu>0;target=176*h*W*W/mu;J=1+-(-target.numerator//target.denominator)
 assert ratio-8*h*W/10**12-88*h*W*W/J>1+mu/2>1
 return dict(terms=M,total=str(total),right=str(sum(right)),Jsum=str(Jsum),crosschecks=cross,left_mesh=str(J),branches=branches)

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
request=json.loads(sys.stdin.readline())
assert request=={}
results={}
for name,fn in [('numeric',midpoint_certificate),('su_response',parameter_certificate),('small_near',small_certificate)]:
 results[name]=fn();print(json.dumps(dict(stage=name,result=results[name])),file=sys.stderr,flush=True)
print(json.dumps(dict(all_passed=True,results=results)))
