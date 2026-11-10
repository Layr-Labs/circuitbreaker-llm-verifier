from functools import lru_cache
import json,itertools,time
def affine_flat_indices(legal,n,d):
 legal_set=set(legal);index={v:i for i,v in enumerate(legal)};planes=[]
 def extend(offset,points,pivots,low):
  if len(pivots)==d:
   planes.append(tuple(index[v]for v in points));return
  k=len(pivots)
  for pivot in range(low,n-(d-k)+1):
   if offset>>pivot&1:continue
   free=[j for j in range(pivot)if j not in pivots]
   for bits in range(1<<len(free)):
    v=(1<<pivot)|sum(((bits>>j)&1)<<b for j,b in enumerate(free));more=[x^v for x in points]
    if all(x in legal_set for x in more):
     extend(offset,points+more,pivots+[pivot],pivot+1)
 for a in legal:extend(a,[a],[],0)
 return planes
def transform_coloring(value,permutation,flip):
 return sum((((value>>i)&1)^flip)<<j for i,j in enumerate(permutation))
def canonical_colorings(group):
 remaining=set(range(4096));representatives=[]
 while remaining:
  value=min(remaining);orbit={transform_coloring(value,p,f)for p,f in group};assert value==min(orbit)
  representatives.append(value);remaining.difference_update(orbit)
 return representatives
def ramsey84_normalized_roots():
 legal=[a for a in range(256)if not(a&1 and a&4)
 and not(a&2 and a&8)and not(a&16 and a&32)]
 planes=affine_flat_indices(legal,8,3);assert len(legal)==108 and len(planes)==14391;ids={};values=[]
 for a in legal:
  key=(int(not(a&1)and not(a&2)),a>>2)
  if key not in ids:ids[key]=len(ids)+1
  values.append(ids[key])
 assert len(ids)==84;edges=sorted({tuple(sorted({values[j]for j in plane}))for plane in planes});assert len(edges)==11349;base=[]
 for edge in edges:base.extend([list(edge),[-v for v in edge]])
 small=[a for a in range(16)if not(a&1 and a&2)]
 small_planes=[sum(1<<small.index(x)for x in points)
 for points in itertools.combinations(small,4)
 if points[0]^points[1]^points[2]^points[3]==0]
 bad={v for v in range(4096)if not any(v&m in(0,m)for m in small_planes)};assert len(small_planes)==39 and len(bad)==48;core=238;orbit=set();stabilizer=[]
 for pair in itertools.permutations((0,1,2)):
  for r1,r2 in itertools.permutations((1,2,3),2):
   for shifts in itertools.product(range(4),repeat=3):
    permutation=[]
    for a in small:
     b=a>>2;free=((b&r1).bit_count()&1)|(((b&r2).bit_count()&1)<<1);free^=shifts[a&3];permutation.append(small.index(pair[a&3]|(free<<2)))
    for flip in(0,1):
     result=transform_coloring(core,permutation,flip);orbit.add(result)
     if result==core:stabilizer.append((permutation,flip))
 assert orbit==bad and len(stabilizer)==96;zeros=[i for i in range(12)if not(core>>i&1)];ones=[i for i in range(12)if core>>i&1];allowed=set()
 for iz,io in itertools.product([-1]+zeros,[-1]+ones):
  code=4095^core
  for j in(iz,io):
   if j>=0:code^=1<<j
  allowed.add(code)
 assert len(allowed)==49;cases=[];remaining=set(allowed)
 while remaining:
  code=min(remaining);current={transform_coloring(code,p,f)for p,f in stabilizer};assert current<=allowed
  cases.append(code);remaining.difference_update(current)
 assert cases==[1809,1811,1813,3857];core0=[ids[(0,t<<2)]for t in small];core1=[ids[(1,t<<2)]for t in small];A=[ids[(0,1|(t<<2))]for t in small]
 B=[ids[(0,2|(t<<2))]for t in small];translations=[[small.index(a^(shift<<2))for a in small]for shift in range(4)]
 canonical_B=canonical_colorings([(p,0)for p in translations]);assert len(canonical_B)==1072;outputs=[]
 for code in cases:
  current=[(p,f)for p,f in stabilizer if transform_coloring(code,p,f)==code];group=[]
  for p,f in current:
   for shift in range(4):
    group.append(([small.index(small[j]^(shift<<2))for j in p],f))
  canonical_A=canonical_colorings(group);clauses=[list(q)for q in base]
  for j,v in enumerate(core0):clauses.append([v if core>>j&1 else-v])
  for j,v in enumerate(core1):clauses.append([v if code>>j&1 else-v])
  def prefix_exclusions(codes,variables):
   def visit(xs,bit,prefix):
    if not xs:
     clauses.append([-variables[j]if b else variables[j]for j,b in prefix]);return
    if len(xs)==1<<(bit+1):return
    assert bit>=0
    for b in(0,1):
     visit([x for x in xs if(x>>bit&1)==b],bit-1,prefix+[(bit,b)])
   visit(codes,11,[])
  prefix_exclusions(canonical_A,A);prefix_exclusions(canonical_B,B);previous=None;nv=84
  for j in range(11,-1,-1):
   a,b=A[j],B[j];prefix=[]if previous is None else[-previous];clauses.append(prefix+[-a,b])
   if j:
    nv+=1;clauses.extend([prefix+[a,b,nv],prefix+[-a,-b,nv]]);previous=nv
  assert nv==95;units={};work=[tuple(q)for q in clauses]
  while True:
   new=[];found={}
   for clause in work:
    if any(abs(x)in units and units[abs(x)]==(x>0)for x in clause):continue
    clause=tuple(x for x in clause if abs(x)not in units);assert clause
    if len(clause)==1:
     variable=abs(clause[0]);value=clause[0]>0;assert variable not in found or found[variable]==value;found[variable]=value
    else:new.append(clause)
   work=new
   if not found:break
   for variable,value in found.items():
    assert variable not in units or units[variable]==value
   units.update(found)
  roots=sorted(set(tuple(sorted(q,key=lambda x:(abs(x),x)))for q in work),
  key=lambda q:(len(q),q))
  outputs.append((code,roots,len(canonical_A),len(units)))
 assert[len(x[1])for x in outputs]==[5130,5049,4773,4408];assert[x[2]for x in outputs]==[270,364,262,49];return outputs

def resolution_refutation(roots,nv,cap=100000):
 clauses=[];watch={v:[]for v in range(-nv,nv+1)if v};values=[0]*(nv+1);levels=[0]*(nv+1);reason=[-1]*(nv+1);activity=[0.]*(nv+1);phase=[-1]*(nv+1);increment=1.
 trail=[];limits=[];head=0;conflicts=0;resolutions=0;decisions=0
 def enqueue(lit,why):
  v=abs(lit);sign=1 if lit>0 else-1
  if values[v]:return values[v]==sign
  values[v]=sign;phase[v]=sign;levels[v]=len(limits);reason[v]=why;trail.append(lit);return True
 def add(clause):
  assert all(1<=abs(l)<=nv for l in clause);assert len({abs(l)for l in clause})==len(clause);ix=len(clauses);clauses.append(list(clause))
  if len(clause)>1:
   watch[clause[0]].append(ix);watch[clause[1]].append(ix)
  return ix
 for row in roots:
  if not row:return{'sat':False,'conflicts':0,'resolutions':0,'decisions':0}
  ix=add(row)
  for l in row:activity[abs(l)]+=2.**-len(row)
  if len(row)==1:
   if not enqueue(row[0],ix):
   # Two contradictory root units resolve to the empty clause.
    assert reason[abs(row[0])]>=0 and clauses[reason[abs(row[0])]]==[-row[0]];return{'sat':False,'conflicts':1,'resolutions':1,'decisions':0}
 def propagate():
  nonlocal head
  while head<len(trail):
   false=-trail[head];head+=1;ws=watch[false];i=0
   while i<len(ws):
    ix=ws[i];row=clauses[ix]
    if row[0]==false:row[0],row[1]=row[1],row[0]
    assert row[1]==false;other=row[0];value=values[abs(other)]*(1 if other>0 else-1)
    if value==1:i+=1;continue
    for j in range(2,len(row)):
     lit=row[j]
     if values[abs(lit)]*(1 if lit>0 else-1)!=-1:break
    else:
     if value==-1:return ix
     assert enqueue(other,ix);i+=1;continue
    row[1],row[j]=row[j],row[1];ws[i]=ws[-1];ws.pop();watch[row[1]].append(ix)
  return-1
 def backtrack(level):
  nonlocal head
  if len(limits)<=level:return
  end=limits[level]
  for lit in trail[end:]:values[abs(lit)]=0;reason[abs(lit)]=-1
  del trail[end:];del limits[level:];head=min(head,end)
 while True:
  bad=propagate()
  if bad>=0:
   conflicts+=1
   if conflicts>cap:raise RuntimeError('finite conflict cap exceeded')
   learnt=set(clauses[bad]);dl=len(limits);assert all(values[abs(l)]*(1 if l>0 else-1)==-1 for l in learnt)
   def resolve(lit):
    nonlocal resolutions,learnt;ri=reason[abs(lit)];assert ri>=0 and lit in clauses[ri]and-lit in learnt
    # This exact set equality is the resolution inference rule.
    learnt=(learnt-{-lit})|(set(clauses[ri])-{lit});assert not any(-x in learnt for x in learnt)
    assert all(values[abs(x)]*(1 if x>0 else-1)==-1 for x in learnt);resolutions+=1
   if dl==0:
    for lit in reversed(trail):
     if-lit in learnt:resolve(lit)
    assert not learnt;return{'sat':False,'conflicts':conflicts,'resolutions':resolutions,'decisions':decisions,'clauses':len(clauses)}
   cursor=len(trail)-1
   while sum(levels[abs(l)]==dl for l in learnt)>1:
    while-trail[cursor]not in learnt:cursor-=1
    lit=trail[cursor];cursor-=1;resolve(lit)
   current=[l for l in learnt if levels[abs(l)]==dl];assert len(current)==1
   first=current[0];others=sorted(learnt-{first},key=lambda l:(-levels[abs(l)],abs(l),l));jump=max((levels[abs(l)]for l in others),default=0);assert jump<dl
   for l in learnt:activity[abs(l)]+=increment
   increment/=.95
   if increment>1.e100:
    activity[:]=[a*1.e-100 for a in activity];increment*=1.e-100
   backtrack(jump);ix=add([first]+others);assert enqueue(first,ix)
   # A fixed restart schedule affects search only, never a proof inference.
   if conflicts%128==0:backtrack(0)
  else:
   free=[v for v in range(1,nv+1)if not values[v]]
   if not free:
    assert all(any(values[abs(l)]*(1 if l>0 else-1)==1 for l in row)for row in roots)
    return{'sat':True,'model':values[1:],'conflicts':conflicts,'resolutions':resolutions,'decisions':decisions}
   variable=max(free,key=lambda v:(activity[v],-v));limits.append(len(trail));decisions+=1;assert enqueue(variable*phase[variable],-1)

def ramsey_direct_audit():
 started=time.monotonic();cases=[]
 for code,roots,canonical_A,units in ramsey84_normalized_roots():
  result=resolution_refutation(roots,95);assert result['sat']is False
  cases.append(dict(code=code,roots=len(roots),canonical_A=canonical_A,units=units,**result))
 assert len(cases)==4;return dict(passed=True,cases=cases,elapsed_seconds=time.monotonic()-started)


CACHE={}
def generated_cover(values,flats,common,hints,full):
 assert len(set(common))==len(common)and all(type(i)is int and 0<=i<len(flats)for i in common);left=full;used=set();steps=0
 for choice in common:
  zero=one=full
  for j in flats[choice]:one&=values[j];zero&=full^values[j]
  left&=full^(zero|one);steps+=1;used.add(choice)
  if not left:return steps
 key=id(flats)
 if key not in CACHE:CACHE[key]=(flats,[sum(1<<j for j in flat)for flat in flats])
 masks=CACHE[key][1]
 while left:
  bit=(left&-left).bit_length()-1;signature=sum(((v>>bit)&1)<<j for j,v in enumerate(values))
  for choice,mask in enumerate(masks):
   hit=signature&mask
   if hit==0 or hit==mask:break
  else:raise AssertionError('no monochromatic witness')
  assert choice not in used;zero=one=full
  for j in flats[choice]:one&=values[j];zero&=full^values[j]
  hit=(zero|one)&left;assert hit>>bit&1;left&=full^hit;used.add(choice);steps+=1;assert steps<=len(flats)+len(common)
 return steps
check_common_cover=generated_cover;check_union_cover=generated_cover

def readonce_seven_audit(envelope):
 started=time.monotonic();legal=[a for a in range(128)if not(a&1 and a&2)and not(a&4 and a&8)];planes=affine_flat_indices(legal,7,3)
 assert len(legal)==72 and len(planes)==5049
 @lru_cache(None)
 def trees(labels):
  if len(labels)==1:return[labels[0]]
  out=[];first=labels[0];rest=labels[1:]
  for k in range(len(rest)):
   for more in itertools.combinations(rest,k):
    left=(first,)+more;right=tuple(x for x in labels if x not in left)
    for a in trees(left):
     for b in trees(right):out.append((a,b))
  return out
 @lru_cache(16384)
 def serial(tree,permutation):
  if isinstance(tree,int):return str(permutation[tree])
  a,b=serial(tree[0],permutation),serial(tree[1],permutation);return '('+','.join(sorted((a,b)))+')'
 permutations=[]
 for swap,a,b in itertools.product((0,1),repeat=3):
  for free in itertools.permutations((4,5,6)):
   p=[0,1,2,3]+list(free)
   if a:p[0],p[1]=p[1],p[0]
   if b:p[2],p[3]=p[3],p[2]
   if swap:p[:4]=p[2:4]+p[:2]
   permutations.append(tuple(p))
 all_trees=trees(tuple(range(7)));representatives=[];seen=set();identity=tuple(range(7))
 for tree in all_trees:
  if serial(tree,identity)in seen:continue
  representatives.append(tree);seen.update(serial(tree,p)for p in permutations)
 assert len(permutations)==48 and len(all_trees)==10395 and len(representatives)==350;serial.cache_clear();trees.cache_clear()
 common=[];assignments=16*3**6;full=(1<<assignments)-1;repeat=full//65535
 signs=[sum(((v>>i)&1)<<v for v in range(16))*repeat for i in range(4)];op_masks=[]
 for i in range(6):
  stride=16*3**i;a=full//((1<<(3*stride))-1)*((1<<stride)-1);op_masks.append((a,a<<stride,a<<(2*stride)))
 for tree in representatives:
  operations=[]
  def compile_tree(t):
   if isinstance(t,int):return t
   a,b=compile_tree(t[0]),compile_tree(t[1]);operations.append((a,b));return 6+len(operations)
  root=compile_tree(tree);assert root==12 and len(operations)==6;values={}
  for j in range(len(legal)):
   a=legal[j];v=[signs[i]if not(a>>i&1)else full^signs[i]for i in range(4)];v.extend(full if a>>i&1 else 0 for i in(4,5,6))
   for(left,right),(aa,oo,xx)in zip(operations,op_masks):
    x,y=v[left],v[right];v.append(((x&y)&aa)|((x|y)&oo)|((x^y)&xx))
   values[j]=v[root]
  generated_cover([values[j]for j in range(len(legal))],planes,common,[],full)
 return{'passed':True,'shape_orbits':350,'assignments_each':assignments,
 'assignments_verified':350*assignments,
 'shared_planes':len(common),
 'elapsed_seconds':time.monotonic()-started}


def modern_scalar_audit():
 import importlib.util,sys
 from pathlib import Path
 sys.dont_write_bytecode=True
 path=Path('/snapshot/files/family.py')
 if not path.exists():path=Path(__file__).with_name('family.py')
 spec=importlib.util.spec_from_file_location('modern_scalar_family',path)
 eng=importlib.util.module_from_spec(spec);sys.modules[spec.name]=eng
 spec.loader.exec_module(eng)
 L=1<<16384;z=eng.ceil_log(L);h=eng.next_pow2(L*L)
 e=eng.next_pow2(L*L*z);ell=eng.next_pow2(8*z*e)
 a=eng.next_pow2(z**3//eng.ceil_log(z));step=eng.ceil_log(a)
 cap=eng.ceil_div(337*L,100*step)
 k,_,ep,q,d,_,Y=eng.split_threshold_resource_candidate(L,h,e,ell,step,cap)
 CW,CY=2749,8491
 assert 31**31>7**10*10**10*21**21 and 15**15>19**2*4*13**13
 kr,_,_,qr,dr,_,Yr=eng.split_threshold_resource_candidate(L,h,e,ell,step,cap,31500,86000)
 assert k<kr and q<=qr and d<=dr and Y<Yr
 assert CW*6250>2*8587831 and CY*625>2*2653383
 plans=[eng.DyadicPrefixPlan(n,C*e*z*z*eng.next_pow2(m),m,e,ell,
  dimension=eng.ceil_div(10*z,7),primitive=True)
  for n,C,m in ((d,CW,ell),(Y,CY,q),(Y,CY,d))]
 Kw=CW*e*z*z*ell;logK=eng.ceil_log(Kw)
 pv=eng.next_pow2(L+2*logK+2*e+33)
 pn=eng.ceil_div(Kw,logK+e+20)*pv
 wsu=eng.seeded_params(pn,Kw//2,eng.ceil_div(10*z,7),ell,e+logK+10)
 degree=max([wsu.s]+[p.su.s for p in plans])
 wsu=eng.SeedParams(wsu.n,wsu.K,wsu.D,wsu.h,wsu.m,wsu.errlog,degree)
 assert wsu.valid() and pv+wsu.d<8*z*e and pv+wsu.d<=ell
 for p in plans:
  p.use_shared_field(degree)
  assert p.valid and p.calls==1 and p.M==eng.next_pow2(p.prefix_output)
  assert p.su.D==eng.ceil_div(10*z,7) and p.su.h<=19*p.su.D
  assert 50*degree<277*e and p.actual_seed<8*z*e and p.actual_seed<=ell
 assert plans[0].su.h<=7*plans[0].su.D and wsu.h<=7*wsu.D
 assert Y-4*cap*h*d==CY*e*z*z*eng.next_pow2(d)
 assert q-2*h*(a-1)*ell==CW*e*z*z*ell
 assert k-h*Y-100*cap*h*a*ell==10*d+CW*e*z*z*ell
 for hh in (4,8,16,128):
  for aa in (2,hh):
   dd=(2*hh*aa+10)*16
   assert 3*hh*dd+8*hh*aa*16+8*hh*16<=4*hh*dd
 return {'passed':True,'L_log2':16384,'shared_ports':4,
  'degree_bit_length':degree.bit_length(),'k_bit_length':k.bit_length(),
  'scope':'One exact scalar shared-port case and bounded ledger examples; no ambient field/input allocation or universal theorem execution.'}


def reciprocal_codec_audit():
 import sys
 eng=sys.modules['modern_scalar_family'];contexts=0;total=0
 def gcd(a,b):
  while b:a,b=b,eng.polynomial_remainder(a,b)
  return a
 for cutoff in (1,2):
  M,parts=eng.small_factor_code(cutoff);d=M.bit_length()-1
  for q in range(4,8):
   t=d+q;J=q//2;A=(1<<(q-1))-(1<<(q-J-1))
   for H in range(1<<d):
    if gcd(H,M)!=1:continue
    expected=set()
    for low in range(1<<q):
     R=eng.binary_polynomial_product((1<<q)|low,M)^H
     for j in range(J):
      lo=(R>>(j+1))&1;hi=(R>>(t-1-j))&1
      if lo!=hi:
       if hi==0:expected.add(R)
       break
    actual=[]
    for index in range(1<<(q-1)):
     Q=eng.canonical_quotient(t,M,H,index)
     if Q is not None:actual.append(eng.binary_polynomial_product(Q,M)^H)
    assert set(actual)==expected and len(actual)==len(set(actual))==A
    contexts+=1;total+=1<<q
 widths=[];edges=0
 for t in (128,160,256,512):
  M,parts,q,A,phi,width=eng.reciprocal_unit_layout(t)
  old=t-eng.primitive_package_saving(t)
  assert width==eng.reciprocal_unit_bits(t)<=old-1
  widths.append({'degree':t,'reference_bits':old,'actual_bits':width})
  for I in (0,1,phi-1,phi):
   for ix in (0,1,A-1,A):
    code=I*A+ix;R=eng.decode_reciprocal_unit_polynomial(t,code)
    if code>=phi*A:assert R is None
    else:
     assert R is not None and R.bit_length()==t+1 and R&1
     assert R.bit_count()%2==1 and gcd(R,M)==1
     for j in range(q//2):
      lo=(R>>(j+1))&1;hi=(R>>(t-1-j))&1
      if lo!=hi:assert lo==1 and hi==0;break
     else:raise AssertionError('missing canonical unequal pair')
    edges+=1
 return {'passed':True,'toy_contexts':contexts,'raw_quotients_enumerated':total,
  'edge_codes':edges,'widths':widths,
  'scope':'Exact toy quotient coverage and bounded code properties; primitive existence for all large degrees is analytic,not a finite search.'}

def run():
 result={'ramsey84_resolution':ramsey_direct_audit(),'readonce7':readonce_seven_audit(None),'modern_shared_scalar':modern_scalar_audit(),'reciprocal_codec':reciprocal_codec_audit()}
 result['all_passed']=True
 result['scope']='Four exact R84 derivations, exhaustive Q7 coverage, and one scoped modern shared-port scalar case and bounded codec checks; complete universal arguments are supplied in Markdown.'
 return result
if __name__=='__main__':
 request=json.loads(input())
 if request!={}:raise ValueError('empty object required')
 print(json.dumps(run(),sort_keys=True))
