"""Exact scoped component checks. Universal claims rest on supplied proofs."""
import json,random,itertools,time
from fractions import Fraction
import importlib.util
from pathlib import Path
import sys
sys.dont_write_bytecode=True
family_path=Path('/snapshot/files/family.py')
if not family_path.exists():family_path=Path(__file__).with_name('family.py')
spec=importlib.util.spec_from_file_location('candidate_family',family_path)
f=importlib.util.module_from_spec(spec)
spec.loader.exec_module(f)
TEST_MODES=tuple(f.rational(i,100)for i in range(15,54))

def rank(rows):
 rows=list(rows);out=0
 while rows:
  a=rows.pop()
  if not a:continue
  pivot=a.bit_length()-1;out+=1
  rows=[b^a if b>>pivot&1 else b for b in rows]
 return out

def xor_basis(basis,mask):
 out=0
 for i,b in enumerate(basis):
  if mask>>i&1:out^=b
 return out

def restrict_form(rows,basis):
 return tuple(sum(((sum(((row&v).bit_count()%2)<<i for i,row in enumerate(rows))&u).bit_count()%2)<<j
      for j,v in enumerate(basis)) for u in basis)

def active_witness(basis,N,S,products=(0,2,4,6)):
 basis=tuple(basis);leaf=(N+(1<<products[-1])-1)//(1<<products[-1]);counts=[1<<(b-c)for c,b in zip(products,products[1:])]
 assert S>leaf and rank(basis)==len(basis)and len(basis)>=S+(S-1)*sum(T-1 for T in counts)
 def entropy(rows,left,right):return rank([(v>>left)&((1<<(right-left))-1)for v in rows])
 index=0;total_prefix_rank=0
 for u,T in enumerate(counts):
  parent_left=index*N//(1<<products[u]);parent_right=(index+1)*N//(1<<products[u]);parent_before=entropy(basis,parent_left,parent_right)
  required=S+(S-1)*sum(T-1 for T in counts[u:]);next_required=S+(S-1)*sum(T-1 for T in counts[u+1:]);assert parent_before>=required
  chosen=None
  for j in range(T):
   child_index=index*T+j;left=child_index*N//(1<<products[u+1]);right=(child_index+1)*N//(1<<products[u+1])
   conditional=entropy(basis,parent_left,right)-entropy(basis,parent_left,left)
   if conditional>=S:chosen=(j,child_index,left,right,conditional);break
  assert chosen is not None
  j,child_index,left,right,conditional=chosen
  forms=[sum(((v>>bit)&1)<<k for k,v in enumerate(basis))for bit in range(left)]
  pivots={}
  for row in forms:
   while row:
    lead=row.bit_length()-1
    if lead not in pivots:pivots[lead]=row;break
    row^=pivots[lead]
  independent=tuple(pivots.values())
  ker=f.kernel(independent+(0,)*(len(basis)-len(independent)),len(basis));spent=len(basis)-len(ker);assert spent<=j*(S-1)
  total_prefix_rank+=spent;basis=tuple(xor_basis(basis,v)for v in ker)
  child_entropy=entropy(basis,left,right);assert child_entropy==conditional
  if child_entropy<next_required:
   outer=entropy(basis,parent_left,parent_right)-child_entropy
   assert child_entropy>=S and outer>=S
   return dict(stopDepth=u,prefixRank=total_prefix_rank,childEntropy=child_entropy,outerConditionalEntropy=outer,S=S)
  index=child_index
 raise AssertionError('continued into impossible leaf')

def run():
 rng=random.Random(20261106);result={};root_cases=0
 for n in range(257):
  for degree in (1,2,3,5,8,10,100,1000):
   a=f.root(n,degree);c=f.ceil_root(n,degree)
   assert a**degree<=n<(a+1)**degree
   assert c==a+(a**degree<n);root_cases+=1
 result['exact_integer_roots']=root_cases
 fields=[]
 for k in (3,5,7,11,13,17,29,37):
  q=f.irreducible_prime(k);assert q.bit_length()==k+1
  value=2
  for _ in range(k):value=f.mul(value,value,q)
  assert value==2 and f.gcd(q,0b110)==1
  fields.append([k,q])
 result['prime_field_constructors']=fields
 linear_cases=0
 for k in (3,5,7):
  q=f.irreducible_prime(k)
  for _ in range(100):
   a=rng.randrange(1<<k);b=rng.randrange(1<<k)
   assert f.trace(a^b,q)==f.trace(a,q)^f.trace(b,q)
   linear_cases+=1
 n=80;m=1;z=2;s,ell=f.code_parameters(n,m,z)
 assert s==37 and (1<<s)>1024*(n+1)*(2*m)**6*(1<<(6*z))
 for _ in range(100):
  x=tuple(rng.randrange(2) for _ in range(n));y=tuple(rng.randrange(2) for _ in range(n))
  index=rng.getrandbits(ell)
  a=f.code_bit(f.coefficients(x,s),s,index);b=f.code_bit(f.coefficients(y,s),s,index)
  xy=tuple(a^b for a,b in zip(x,y))
  assert f.code_bit(f.coefficients(xy,s),s,index)==a^b
  linear_cases+=1
 result['trace_and_actual_local_code_linearity']=linear_cases
 result['local_code_parameters']={'n':n,'m':m,'epsilon_exponent':z,'field_degree':s,'index_bits':ell}
 q=f.irreducible_prime(7);Q=128;sets=[]
 for a in range(4):
  for b in range(4):sets.append({alpha*Q+(f.mul(b,alpha,q)^a) for alpha in range(74)})
 graph_cases=0
 for i,j in itertools.combinations(range(len(sets)),2):
  assert len(sets[i]&sets[j])<=1;graph_cases+=1
 result['finite_linear_graph_intersections']=graph_cases
 spaces=0;rank_checks=0
 for trial in range(64):
  n=12;c=1+trial%3;generators=[]
  for _ in range(4):
   rows=[0]*n
   for i in range(c):
    for j in range(i+1,n):
     if rng.randrange(2):rows[i]^=1<<j;rows[j]^=1<<i
   generators.append(tuple(rows))
  forms=[]
  for mask in range(16):
   rows=[0]*n
   for j,g in enumerate(generators):
    if mask>>j&1:rows=[a^b for a,b in zip(rows,g)]
   forms.append(tuple(rows))
  initial=max(rank(a) for a in forms);basis=tuple(1<<i for i in range(n))
  while True:
   restricted=[restrict_form(a,basis) for a in forms]
   ranks=[rank(a) for a in restricted];maximum=max(ranks)
   if maximum==0:break
   kernel=f.kernel(restricted[ranks.index(maximum)],len(basis))
   nextbasis=tuple(xor_basis(basis,v) for v in kernel)
   for a in forms:
    assert 2*rank(restrict_form(a,nextbasis))<=maximum;rank_checks+=1
   basis=nextbasis
  assert n-len(basis)<=2*initial
  assert all(rank(restrict_form(a,basis))==0 for a in forms);spaces+=1
 result['common_radical_spaces']=spaces;result['rank_halving_checks']=rank_checks
 permutations=list(itertools.permutations(range(3)));stats={};support_cases=0
 for trial in range(1000):
  k=rng.randrange(83,501);d=(2*k+4)//5+10;end=k-d;L=k-2*d+4;U=d-4
  support={0,end}
  for i in range(1,end):
   if trial%2==0 and L<=i<=U:continue
   if rng.randrange(2):support.add(i)
  ordered=sorted(support)
  gaps=[(ordered[(i+1)%len(ordered)]-x-1)%k for i,x in enumerate(ordered)]
  if max(gaps)!=d-1 or sum(g==d-1 for g in gaps)!=1:continue
  middle=sorted(i for i in support if L<=i<=U)
  if middle:columns=[0,middle[0]+1,end+2];case='middle'
  else:
   a=max(i for i in support if i<L);b=min(i for i in support if i>U);A=a+1;C=end-b+1
   if A>C:columns=[0,end+1,a+2];case='first_island'
   elif C>A:columns=[b,a+1,end+2];case='second_island'
   else:columns=[0,b+1,a+2];case='equal_islands'
  assert len(set(columns))==3
  for rho in range(k):
   matrix=[[((b-u-rho)%k) in support for b in columns] for u in range(3)]
   match=sum(all(matrix[u][p[u]] for u in range(3)) for p in permutations)
   assert match==(1 if rho==0 else 0)
  support_cases+=1;stats[case]=stats.get(case,0)+1
 assert len(stats)==4
 result['cubic_support_patterns']=support_cases;result['isolation_cases']=stats
 interval_cases=0
 for n in range(1,257):
  K=f.ceil_root(n,100)
  assert (K-1)**100<n<=K**100
  for t in (2,4,8):
   level=0
   while t**level<=n:
    intervals=[f.interval_bounds(n,t,level,i) for i in range(t**level)]
    assert intervals[0][0]==0 and intervals[-1][1]==n
    assert all(a[1]==b[0] for a,b in zip(intervals,intervals[1:]))
    lengths=[b-a for a,b in intervals]
    assert min(lengths)>=1 and max(lengths)-min(lengths)<=1
    if level:
     for i in range(t**(level-1)):
       parent=f.interval_bounds(n,t,level-1,i)
       children=intervals[i*t:(i+1)*t]
       assert children[0][0]==parent[0] and children[-1][1]==parent[1]
    interval_cases+=len(intervals);level+=1
 result['actual_interval_partition_nodes']=interval_cases
 rounding_cases=0
 for L in range(1,2049):
  B=f.ceil_log2(L+1);U=f.ceil_root(L**3,4);V=f.ceil_root(B,4)
  Z=U*V;E=13*Z;base=f.ceil_root(L,4)*f.ceil_root(B**3,4)
  A=(27*base+24)//25
  assert 25*(A-1)<27*base<=25*A
  ell_log=f.ceil_root(L,2)*f.ceil_root(B,2)
  ratio=(L+B-1)//B;qbase=f.ceil_root(ratio,2);q=(22*qbase+4)//5
  assert (U-1)**4<L**3<=U**4 and (V-1)**4<B<=V**4
  assert (1<<(B-1))<L+1<=(1<<B)
  assert (qbase-1)**2<ratio<=qbase*qbase and 5*(q-1)<22*qbase<=5*q and E==13*Z
  pexp=(6301*Z+9999)//10000
  assert 10000*(pexp-1)<6301*Z<=10000*pexp
  assert A>=1 and ell_log>=1
  rounding_cases+=1
 for denominator in (1,2,3,10,100,1000,10000,100000):
  L=(1<<160)*denominator**8;B=f.ceil_log2(L+1)
  Z=f.ceil_root(L**3,4)*f.ceil_root(B,4)
  assert 3+13*Z<=(L-1)//denominator
  rounding_cases+=1
 result['quarter_entropy_rounding_cases']=rounding_cases
 extension_cases=0;cycles=[]
 for degree in (1,2,3,4):
  ff,modulus,g=f.singer(3,degree);Q=8**degree
  assert f.extension_irreducible(modulus,3,ff)
  seen=set();v=1
  for _ in range(Q-1):
   seen.add(v);v=f.extension_mul(v,g,3,modulus,ff)
  assert v==1 and len(seen)==Q-1
  cycles.append([degree,Q-1])
  for _ in range(200):
   a,b,c=(rng.randrange(Q) for _ in range(3))
   assert f.extension_mul(a,b^c,3,modulus,ff)==f.extension_mul(a,b,3,modulus,ff)^f.extension_mul(a,c,3,modulus,ff)
   assert f.extension_power(a,Q,3,modulus,ff)==a
   extension_cases+=1
 result['singer_nonzero_cycles']=cycles
 result['extension_field_identity_cases']=extension_cases
 monomial_cases=0
 for degree in range(1,7):
  for h in range(1,8):
   exps=list(f.monomial_exponents(degree,h))
   assert len(exps)==f.comb(h+degree-1,degree)==len(set(exps))
   assert all(len(e)==degree and sum(e)<h for e in exps)
   monomial_cases+=1
 result['graded_monomial_enumeration_cases']=monomial_cases
 optimal_linearity=0
 for degree in (1,2,3):
  for _ in range(100):
   n=12;h=14;seed=rng.randrange(1<<(3*(degree+1)))
   x=tuple(rng.randrange(2) for _ in range(n));y=tuple(rng.randrange(2) for _ in range(n))
   xy=tuple(a^b for a,b in zip(x,y))
   a=f.optimal_extract_with_parameters(x,seed,3,degree,h,4)
   b=f.optimal_extract_with_parameters(y,seed,3,degree,h,4)
   c=f.optimal_extract_with_parameters(xy,seed,3,degree,h,4)
   assert c==tuple(v^w for v,w in zip(a,b))
   optimal_linearity+=1
 result['optimal_extractor_finite_linearity_cases']=optimal_linearity
 F=Fraction
 for mode in TEST_MODES:
  r=f.regime(mode);beta=F(*r['beta']);alpha=F(*r['alpha'])
  zeta=F(*r['zeta']);sigma=F(2000001,1000000)*zeta
  near=sigma+F(1,1000000);group=F(*r['group']);slice_exp=F(*r['slice'])
  guard=F(*r['seed_guard']);aa=F(*r['arity'])
  assert group+slice_exp==F(999999,1000000) and group==1-beta
  assert near<slice_exp/2 and slice_exp<beta
  assert beta+1/alpha<1 and 1+4*alpha*beta<guard<6
  assert guard*aa/group==F(9999,10000)
  assert zeta-F(1,100000000)>zeta-F(2,100000000)>F(0)
 result['rational_rate_sample_contexts']=len(TEST_MODES)
 output_roots=0
 for r in (0,1,2,3,4,7,8):
  m=f.root(r**999999,1000000)
  assert m**1000000<=r**999999<(m+1)**1000000
  output_roots+=1
 result['larger_smalloutput_integer_roots']=output_roots
 reservoir_plans=0;reservoir_rounds=0
 for mode in TEST_MODES:
  for bits in (4000000,):
   r=f.regime(mode);bn,bd=r['beta'];gn,gd=r['group']
   K=1<<bits;g=1<<((bits*gn)//gd);w=1<<((bits*r['slice'][0])//r['slice'][1])
   R=1<<((bits*bn)//bd);M=K//4
   assert w<min(R,M)
   for initial in (1,2,g):
    J=f.ceil_log2(initial);count=initial;spent=0;steps=0
    assert (2*initial+J)*w<=K//4
    if initial<=3:
     while count>1:
       spent+=count*w;count=(count+1)//2;steps+=1
    else:
     spent=(2*initial+J-2-(initial-1).bit_count())*w
     steps=J
     for j in (0,1,J//2,J-1):
       power=1<<j
       count=(initial+power-1)//power
       assert (count+1)//2==(initial+2*power-1)//(2*power)
    assert steps==J and spent<=(2*initial+J)*w
    assert K-spent>=K//2+M
    reservoir_plans+=1;reservoir_rounds+=steps
 result['reservoir_row_count_budget_cases']=reservoir_plans
 result['reservoir_row_count_rounds']=reservoir_rounds
 kernel_fibers=0;kernel_sources=0
 for trial in range(32):
  dimension=8;observations=tuple(rng.randrange(1<<dimension) for _ in range(1+trial%5))
  kernel=f.kernel(observations+(0,)*(dimension-len(observations)),dimension)
  kernel_values={xor_basis(kernel,m) for m in range(1<<len(kernel))}
  fibers={}
  for x in range(1<<dimension):
   z=sum(((x&row).bit_count()%2)<<j for j,row in enumerate(observations))
   fibers.setdefault(z,set()).add(x)
  assert len(kernel)==dimension-rank(observations)
  assert len(fibers)*len(kernel_values)==1<<dimension
  for points in fibers.values():
   representative=min(points)
   assert {x^representative for x in points}==kernel_values
   kernel_fibers+=1
  kernel_sources+=1
 result['reservoir_affine_kernel_sources']=kernel_sources
 result['reservoir_affine_kernel_fibers']=kernel_fibers
 confidence_geometry=0
 for ell_log in list(range(1,601))+[1999999,2000000,2000001,3999999,4000000,4000001,6000000]:
  unit=(ell_log+1999999)//2000000;separation=2000001*unit
  gap=separation-ell_log
  assert gap>=unit and 2000000*(unit-1)<ell_log<=2000000*unit
  if ell_log>=2000001:assert gap<=2*ell_log
  for remaining in range(1,21):
   log1=4000002*remaining*unit;log2=log1+separation
   logD=ell_log+log2;ancestor=4000002*(remaining+1)*unit
   assert log2-log1-ell_log==gap
   assert ancestor-logD==gap
   assert 3*ell_log+log2==(4000002*remaining+2000001)*unit+3*ell_log
   confidence_geometry+=1
 result['quantized_confidence_exponent_geometry_cases']=confidence_geometry
 degree_rounding=0
 for mode in TEST_MODES:
  an,ad=f.regime(mode)['alpha']
  for L in range(1,2049):
   for w in sorted(set((1,2,3,L))):
    d=(an*L+ad*w-1)//(ad*w)
    assert ad*w*(d-1)<an*L<=ad*w*d
    degree_rounding+=1
 result['sharpordinary_degree_rounding_cases']=degree_rounding
 fallback=0
 for n in range(257):
  assert f.tree_parameters(n) is None
  for pattern in (0,1,2):
   x=tuple(0 if pattern==0 else 1 if pattern==1 else rng.randrange(2) for _ in range(n))
   assert f.evaluate(x)==0;fallback+=1
 result['small_total_fallback_inputs']=fallback
 import sys
 component_cases=0;accepted_parameter_cases=0
 target=f.optimal_parameters.__wrapped__.__code__
 for mode in TEST_MODES:
  for tree_log in (2,3,4,8,16,32,64,128):
   n=1<<tree_log
   for w in sorted(set((1,max(1,tree_log//2),tree_log))):
    seen={}
    def trace_component(frame,event,arg):
     if frame.f_code is target:
       if event=='return':
        seen.update({k:frame.f_locals[k] for k in ('L','w','d')})
       return trace_component
     return None
    f.optimal_parameters.cache_clear()
    sys.settrace(trace_component)
    try:
     actual=f.optimal_parameters(n,1<<w,mode=mode)
    finally:
     sys.settrace(None)
    local_log=f.ceil_log2(n+1)
    assert local_log==tree_log+1 and seen['L']==local_log and seen['w']==w
    an,ad=f.regime(mode)['alpha']
    expected=(an*local_log+ad*w-1)//(ad*w)
    assert seen['d']==expected
    if actual is not None:
     assert actual[2]==expected
     accepted_parameter_cases+=1
    component_cases+=1
 result['actual_component_local_log_power_two_cases']=component_cases
 result['actual_component_local_log_accepted_cases']=accepted_parameter_cases
 size_guard_cases=0;saved_root=f.root
 def forbidden_root(*args,**kwargs):raise AssertionError('power/root reached before context guard')
 f.root=forbidden_root
 try:
  bad=(1000000,2000000)
  for constructor,args in ((f.optimal_parameters,(16,8)),(f.small_parameters,(16,8)),(f.reservoir_parameters,(16,8,1,1,1)),(f.block_parameters,(16,8))):
   constructor.cache_clear()
   assert constructor(*args,mode=bad)is None;size_guard_cases+=1
 finally:f.root=saved_root
 result['rational_context_size_guard_cases']=size_guard_cases
 public_cases=0
 for Z in (5,13,29):
  low=(527*Z+499)//500;high=(33743*Z+9999)//10000
  for W in (low,(low+high)//2,high):
   mode=f.level_mode(W,Z);n=1<<256;K=1<<W;seen={}
   def trace_public(frame,event,arg):
    if frame.f_code is target:
     if event=='return':seen.update({k:frame.f_locals[k]for k in ('L','w','d')})
     return trace_public
    return None
   f.optimal_parameters.cache_clear();sys.settrace(trace_public)
   try:f.optimal_parameters(n,K,mode=mode)
   finally:sys.settrace(None)
   an,ad=f.regime(mode)['alpha'];expected=(an*257+ad*W-1)//(ad*W)
   assert seen==dict(L=257,w=W,d=expected);public_cases+=1
 result['actual_public_rational_rate_parameter_cases']=public_cases
 accepted_parameters=0;rejected_parameters=0
 # Implicit ambient-size INTEGER computations only: never allocate n input bits,
 for logarithm,entropy,mode,should_accept in ((20000,1800,(3,20),True),(40000,3600,(19,100),True),(20000,4000,(3,10),True),(20000,3000,(2,5),True),(40000,4000,(12,25),True),(8192,2100,(53,100),True),(1024,900,(1,2),False),(2048,1500,(1,2),True),(4096,3000,(1,2),True),(1024,900,(13,25),False),(2048,1500,(13,25),True),(4096,3000,(13,25),True),(1024,1000,(3,20),False),(1024,1000,(13,25),False),(1024,1000,(3,5),False)):
  n=1<<logarithm;K=1<<entropy;params=f.optimal_parameters(n,K,mode=mode)
  assert (params is not None)==should_accept
  if params is None:rejected_parameters+=1;continue
  m,s,d,h,seedlen=params;bn,bd=f.regime(mode)['beta'];an,ad=f.regime(mode)['alpha']
  assert m**bd<=K**bn<(m+1)**bd
  assert d==(an*(logarithm+1)+ad*entropy-1)//(ad*entropy)
  assert f.comb(h+d-1,d)>=n and (1<<s)**d<=(n+1)**f.field_cap(mode) and h**d<=(n+1)**f.field_cap(mode)
  assert (1<<s)>(1<<30)*h*d*s*(32*m)**4 and h<(1<<s)-1
  assert K>=400*m*h*d*s*s+f.ceil_log2(32*m)+3 and seedlen==(d+1)*s
  accepted_parameters+=1
 result['implicit_large_su_parameter_positive_cases']=accepted_parameters
 result['implicit_large_su_parameter_negative_cases']=rejected_parameters
 variable_cases=0;variable_branches=[0,0];partition_cases=0
 for n in range(2,100):
  for sequence in ((1,2,1,3),(2,1,2,1),(3,1,1,2)):
   logs=[0];ar=[];j=0
   while True:
    a=sequence[j%len(sequence)];ar.append(a)
    if logs[-1]+a>n.bit_length()-1:break
    logs.append(logs[-1]+a);j+=1
   for level in range(len(logs)):
    pieces=[f.level_interval_bounds(n,logs,level,i) for i in range(1<<logs[level])]
    assert pieces[0][0]==0 and pieces[-1][1]==n
    assert all(pieces[i][1]==pieces[i+1][0] for i in range(len(pieces)-1))
    assert max(b-a for a,b in pieces)-min(b-a for a,b in pieces)<=1
    if level+1<len(logs):
     for i,(start,end) in enumerate(pieces):
       children=[f.level_interval_bounds(n,logs,level+1,(i<<ar[level])+r) for r in range(1<<ar[level])]
       assert children[0][0]==start and children[-1][1]==end
    partition_cases+=1
   for gapvalue in range(max(ar)+1,max(ar)+4):
    indices=f.active_level_indices(logs,gapvalue)
    assert indices[0]==0 and len(indices)==len(set(indices))
    assert all(logs[indices[i+1]]-logs[indices[i]]<=gapvalue+max(ar) for i in range(len(indices)-1))
    for threshold in range(1,n+1):
     eligible=[j for j in indices if n//(1<<logs[j])>=threshold]
     h=eligible[-1];pos=indices.index(h);maximum=(n+(1<<logs[h])-1)//(1<<logs[h])
     if pos+1<len(indices):
       nxt=indices[pos+1];assert n//(1<<logs[nxt])<threshold;variable_branches[0]+=1
     else:
       assert logs[-1]-logs[h]<gapvalue
       assert n<(1<<(logs[-1]+ar[-1]));variable_branches[1]+=1
     assert maximum<=(1<<(gapvalue+max(ar)))*(threshold+1)+1
     variable_cases+=1
 assert all(variable_branches)
 result['local_variable_partition_levels']=partition_cases
 result['local_variable_active_leaf_cases']=variable_cases
 result['local_variable_leaf_branch_counts']=variable_branches
 interleaved_cases=0;positive_cases=0
 for trial in range(32):
  basis=tuple((1<<j)|(rng.randrange(16)<<12) for j in range(12))
  def conditional_dimension(vectors):
   return rank([v&4095 for v in vectors])-rank([v&255 for v in vectors])
  for step in range(5):
   before=conditional_dimension(basis);positive_cases+=before>0
   form=1<<rng.randrange(8)
   row=sum(((v&form).bit_count()%2)<<j for j,v in enumerate(basis))
   kernel=f.kernel((row,)+(0,)*max(0,len(basis)-1),len(basis));basis=tuple(xor_basis(basis,v) for v in kernel)
   assert conditional_dimension(basis)==before
   form=rng.randrange(1<<16)
   row=sum(((v&form).bit_count()%2)<<j for j,v in enumerate(basis))
   kernel=f.kernel((row,)+(0,)*max(0,len(basis)-1),len(basis));basis=tuple(xor_basis(basis,v) for v in kernel)
   after=conditional_dimension(basis);assert before-1<=after<=before
   interleaved_cases+=1
 assert positive_cases>32
 result['interleaved_affine_restriction_cases']=interleaved_cases
 result['nonzero_conditional_dimension_cases']=positive_cases
 # Symbolic certificate fixtures do not assert primitive availability or a
 toy=dict(n=1 << 100, H=4, ell=2, P=512, active=(0, 2, 4), products=(0, 20, 40, 60, 80), arities=(20,) * 5, blocks=(1,) * 4, thresholds=(1024,) * 4, response_rows=(2,) * 4)
 assert f.geometry(toy,1<<99) is not None
 assert f.geometry(toy,1<<70) is None
 assert f.geometry(dict(toy,thresholds=(1<<90,)*4),1<<99) is None
 assert f.geometry(dict(toy,blocks=(1<<20,)*4),1<<99) is None
 result['finite_geometry_positive_symbolic_fixtures']=1
 result['finite_geometry_negative_symbolic_fixtures']=3
 good=((40,96),(2,12));tiny=F(1,1<<1024)
 assert f.output_events(4,1,2,512,good,tiny) is not None
 assert f.output_events(4,1,2,512,((20,96),(2,12)),tiny) is None
 assert f.output_events(4,1,2,512,((40,60),(2,12)),tiny) is None
 assert f.output_events(4,1,2,64,good,tiny) is None
 assert f.output_events(4,1,2,512,good,F(1,100)) is None
 result['finite_output_positive_symbolic_fixtures']=1
 result['finite_output_negative_symbolic_fixtures']=4
 coverage=0
 for lower,upper in ((F(-2,3),F(4,3)),(F(3,50),F(2,25))):
  for cap in range(1,20):
   expected={F(n,d)for d in range(1,cap+1)for n in range(-2*cap,2*cap+1)if lower<=F(n,d)<=upper}
   actual=f.rationals(lower,upper,cap)
   assert tuple(sorted(expected))==actual and len(set(actual))==len(actual)
   coverage+=1
 assert tuple(f.configurations(5))==(f.BASE,) and f.analytic_prefilter(f.BASE)
 assert not f.analytic_prefilter(f.BASE[:6]+(F(4),))
 result['public_grid_denominator_coverage_cases']=coverage
 ee=f.primitive_error(1024,2,100000000,4,3)
 assert ee['near']==F(1,1<<1023) and ee['cubic']==F(1,1<<199)
 assert ee['block']==F(5,1<<1022)+F(1,1<<197)
 assert f.primitive_error(1024,2,100000000,4,1)['block']==0
 assert f.primitive_error(0,2,100000000,4,3) is None
 result['finite_primitive_error_conditional_arithmetic_fixtures']=4
 assert f.gcd(0b1010,0b110)==0b110
 saved_configs=f.configurations;f.tree_parameters.cache_clear()
 def forbidden_grid(*args):raise AssertionError('early fallback reached grid enumeration')
 try:
  f.configurations=forbidden_grid
  for N in range(129):assert f.evaluate(tuple(i&1 for i in range(N)))==0
 finally:f.configurations=saved_configs
 result['public_grid_fallback_before_enumeration_lengths']=129
 h,a,wm,w0,X,pc,D=F(89753,1250000),F(550791,1000000),F(1042927,1000000),F(168307,50000),F(521463,500000),F(1376943,10000000),F(4409067,1000000)
 assert f.analytic_prefilter((h,a,wm,w0,X,pc)) and D>w0+X
 M=1048576;step=(w0-wm)/M
 V=wm.denominator*step.denominator;U=wm.numerator*step.denominator;inc=step.numerator*wm.denominator
 ap,aq=a.numerator,a.denominator;total=0;cross=0
 for j in range(M):
  z=U*aq-ap*V
  num=10**12*100000*U*z*z
  den=V*(100000*aq*aq*U*U+300000*ap*aq*U*V+4*ap*V*z)
  assert den>0
  integer=num//den;assert integer*den<=num<(integer+1)*den
  if j in (0,1,M//2,M-1):
   w=wm+j*step;beta=a/w;alpha=1/(1-beta)+F(1,100000)
   value=10**12*w*(1-beta)/(1+4*alpha*beta)
   assert integer==value.numerator//value.denominator;cross+=1
  total+=integer;U+=inc
 assert total==785745745700283704
 ratio=8*h*step*F(total,10**12)
 assert ratio==F(20480004810115038364036863107,20480000000000000000000000000)>1
 result['current_exact_cost_grid_explicit_box_certificate_terms']=M
 result['current_exact_cost_grid_explicit_box_certificate_sum']=str(total)
 result['integer_box_independent_rational_crosschecks']=cross
 endpoints=((19,100),(3,10),(2,5),(12,25),(53,100))
 for mode in endpoints:
  sn,sd=f.regime(mode)['seed_guard']
  assert F(sn,sd)<f.field_cap(mode)-F(1,20)
 result['adaptive_field_cap_endpoint_margin_cases']=len(endpoints)
 result['adaptive_field_cap_additional_implicit_positive_cases']=5
 minimum_cases=0
 for width in (64,512,2048):
  toy=dict(n=1 << 100, H=4, ell=2, P=width, active=(0, 2, 4), products=(0, 20, 40, 60, 80), arities=(20,) * 5, blocks=(1,) * 4, thresholds=(1024,) * 4, response_rows=(2,) * 4)
  dim=f.minimum_geometry_dimension(toy)
  assert 1<=dim<=toy['n'] and f.geometry(toy,dim)is not None
  assert f.geometry(toy,dim-1)is None
  assert (dim-1)%(1+2*((1<<40)-1))==0
  minimum_cases+=1
 result['minimum_geometry_positive_symbolic_cases']=minimum_cases
 result['minimum_geometry_boundary_negative_cases']=minimum_cases
 result['improved_active_coefficient_symbolic_cases']=minimum_cases
 rank_cases=[active_witness(tuple(1<<i for i in range(576)),576,S,(0,2,4,6,8))for S in (4,13,40)]
 base=tuple(1<<i for i in range(3))+tuple(1<<i for i in range(36,66))
 rank_cases.append(active_witness(base,144,4))
 rank_cases.append(active_witness(tuple(v^(1<<(100+j%44))for j,v in enumerate(base)),144,4))
 last=tuple(1<<(start+j)for start in (0,36,72)for j in range(3))+tuple(1<<i for i in range(108,130))
 rank_cases.append(active_witness(last,144,4))
 assert {q['stopDepth']for q in rank_cases}=={0,1,2} and rank_cases[-1]['prefixRank']==9
 result['integer_active_budget_affine_rank_cases']=len(rank_cases)
 result['integer_active_budget_nonzero_prefix_cases']=sum(q['prefixRank']>0 for q in rank_cases)
 result['integer_active_budget_distinct_stop_depths']=len({q['stopDepth']for q in rank_cases})
 cost_cases=0
 for minimum in (1,17,100):
  for maximum in (minimum,minimum+1,minimum+73):
   current=maximum;cost=0;floor_seen=current==minimum
   for b in (1,2,3,4,7,16,31,64):
    nxt=f.next_entropy_log(current,minimum,b);q=(2*b-1).bit_length()
    assert nxt==max(minimum,current-q) and minimum<=nxt<=current
    cost+=q
    if not floor_seen and nxt>minimum:assert cost==maximum-nxt
    floor_seen=floor_seen or nxt==minimum
    current=nxt;cost_cases+=1
 result['actual_integer_entropy_cost_recurrence_cases']=cost_cases
 local_product_cases=0
 for width in (64,512,2048):
  toy=dict(n=1 << 100, H=4, ell=2, P=width, active=(0, 2, 4), products=(0, 20, 40, 60, 80), arities=(20,) * 5, blocks=(1,) * 4, thresholds=(1024,) * 4, response_rows=(2,) * 4)
  current=f.minimum_geometry_dimension(toy)
  rank_cost=4*(1<<20)*2*width;res=5*2*width;Q=rank_cost+1025;factor=1
  for j in range(4):
   nxt=factor*2
   Q=max(Q,nxt*(rank_cost+1025),factor*(rank_cost+2*(1<<20)*(1024+res+1)))
   factor=nxt
  Q=max(Q,factor*(rank_cost+(1<<20)+1))
  old=(1+2*((1<<40)-1))*(Q-1)+1
  assert current<old and f.geometry(toy,current)is not None and f.geometry(toy,current-1)is None
  local_product_cases+=1
 result['first_possible_descent_strict_dimension_cases']=local_product_cases
 rank_budget_cases=0
 for arities in ((3,3,3,3,3),(2,3,4,5,20),(10,11,12,13,99)):
  for width in (64,512):
   p0=dict(H=4, ell=2, P=width, arities=arities, response_rows=(2,) * 4)
   expected=sum((1<<a)-1 for a in arities[:3])*2*width
   assert f.response_rank_budget(p0)==expected<4*(1<<max(arities))*2*width
   rank_budget_cases+=1
 result['exact_sibling_response_rank_budget_cases']=rank_budget_cases
 least_degree_cases=0;strict_seed_cases=0
 for L0,W0,mode in ((20000,1800,(3,20)),(40000,3600,(19,100)),(20000,4000,(3,10)),(20000,3000,(2,5)),(40000,4000,(12,25)),(2048,1500,(1,2)),(4096,3000,(1,2)),(2048,1500,(13,25)),(4096,3000,(13,25))):
  n0=1<<L0;K0=1<<W0;pp=f.optimal_parameters(n0,K0,mode=mode);assert pp
  m0,s0,d0,h0,seed0=pp;C0=(1<<30)*h0*d0*(32*m0)**4
  assert f.is_prime(s0) and s0%2 and (1<<s0)>C0*s0
  previous=s0-1
  while previous>=3 and not f.is_prime(previous):previous-=1
  assert previous<3 or (1<<previous)<=C0*previous
  old_s=f.prime_at_least((C0*(f.ceil_log2(n0+1)+1)**2).bit_length())
  assert s0<=old_s and seed0==(d0+1)*s0
  strict_seed_cases+=s0<old_s;least_degree_cases+=1
 assert strict_seed_cases==8
 result['least_valid_odd_prime_degree_cases']=least_degree_cases
 result['least_valid_degree_strict_seed_savings_cases']=strict_seed_cases
 response_target_cases=0;response_seed_strict=0
 for L0,W0,mode in ((20000,1800,(3,20)),(40000,3600,(19,100)),(20000,4000,(3,10)),(20000,3000,(2,5)),(40000,4000,(12,25)),(2048,1500,(1,2)),(4096,3000,(1,2)),(2048,1500,(13,25)),(4096,3000,(13,25))):
  n0=1<<L0;K0=1<<W0;default=f.constant_parameters(n0,K0,mode=mode);assert default
  target=f.response_parameters(n0,K0,1<<64,mode=mode);assert target
  assert target[0]==1<<64 and target[2:4]==default[2:4] and target[-1]<=default[-1]
  assert f.constant_parameters(n0,K0,mode=mode)==default
  m0,s0,d0,h0,seed0=target;assert (1<<s0)>(1<<30)*h0*d0*s0*(32*m0)**4
  assert K0>=400*m0*h0*d0*s0*s0+f.ceil_log2(32*m0)+3
  response_seed_strict+=target[-1]<default[-1];response_target_cases+=1
  assert f.response_parameters(n0,K0,default[0]+1,mode=mode)is None
  assert f.response_parameters(n0,K0,0,mode=mode)is None
 assert response_seed_strict==9
 result['actual_response_output_target_cases']=response_target_cases
 result['actual_response_target_strict_seed_cases']=response_seed_strict
 terminal_rank_cases=0
 for width in (64,128,512):
  pp=dict(n=1 << 40, H=2, ell=2, P=width, active=(0, 1, 2), products=(0, 10, 20), arities=(10, 10, 99), blocks=(1, 1), thresholds=(1024, 1024), response_rows=(2,) * 2)
  dim=f.minimum_geometry_dimension(pp);assert f.geometry(pp,dim)is not None and f.geometry(pp,dim-1)is None
  assert f.response_rank_budget(pp)==1023*2*width
  terminal_rank_cases+=1
 result['last_internal_level_response_rank_exclusion_cases']=terminal_rank_cases
 rr=f.regime((53,100));bb=F(53,100);aa=F(*rr['alpha']);gg=F(*rr['seed_guard']);ep=F(1,100000)
 assert aa==F(100,47)+ep and gg==1+4*aa*bb+ep<6-F(1,20)
 assert 1-bb-1/aa>=ep*F(2209,10000)/(1+F(85,100)*ep)
 zz=F(*rr['zeta']);nn=F(2000001,1000000)*zz+F(1,1000000)
 assert (bb-F(1,1000000))/2-nn>F(36,100000000) and F(*rr['arity'])>F(2,25)
 result['expanded_rate_uniform_endpoint_margin_cases']=1
 expanded=F(13,25)<F(550791,1000000)/F(1042927,1000000)<F(53,100)
 assert expanded
 result['new_box_configuration_exceeds_previous_rate_cap']=expanded
 assert f.response_row_counts((1,2,1),(1,1,2),3)==(5,1,1)
 assert f.response_row_counts((10,10,10),(1,1,1),3)==(2048,1024,1)
 weighted_cases=0
 for width in (16,64,128):
  pp=dict(n=1<<60,H=3,ell=4096,P=width,active=(0,1,3),products=(0,10,20,30),arities=(10,10,10,99),blocks=(1,1,1),thresholds=(1024,)*3,response_rows=(2048,1024,1))
  old=pp|dict(response_rows=(4096,)*3)
  assert f.response_rank_budget(pp)<f.response_rank_budget(old)
  assert f.response_resiliency_budget(pp)<f.response_resiliency_budget(old)
  dn=f.minimum_geometry_dimension(pp);do=f.minimum_geometry_dimension(old)
  assert 0<dn<do<=pp['n'] and f.geometry(pp,dn)is not None and f.geometry(pp,dn-1)is None
  weighted_cases+=1
 result['actual_public_response_row_count_fixtures']=2
 result['weighted_response_strict_dimension_cases']=weighted_cases
 root_cases=0
 for e in range(1,101):
  for d in range(1,16):
   for n in ((1<<e)-1,1<<e,(1<<e)+1):
    r=f.root(n,d);assert r**d<=n<(r+1)**d;root_cases+=1
 result['exact_power_root_and_neighbor_cases']=root_cases
 arity_cases=0
 for b in range(1,13):
  for e in range(1,11):
   for L in range(2,17):
    for upper in range(1,8):
     choices=[A for A in range(max(1,f.ceil_log2(b)),upper+1)if f.comb(1<<A,b)<=(1<<e)//(L+1)]
     assert f.binomial_arity(1<<e,L,b,upper)==(max(choices)if choices else None)
     arity_cases+=1
 result['exact_binomial_arity_maximum_cases']=arity_cases
 result['all_passed']=True
 result['scope']='Finite BASE component predicates, conditional error arithmetic, scoped symbolic geometry/probability fixtures and public-grid fallback only; no complete accepted tree or universal extraction/circuit acceptance.'
 return result

if __name__=='__main__':
 request=json.loads(input())
 if request!={}:raise ValueError('empty input object required')
 print(json.dumps(run(),sort_keys=True))
