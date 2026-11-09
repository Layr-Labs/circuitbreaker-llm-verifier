"""Exact exhaustive terminal cover and rational potential witnesses only."""
import json,itertools,time,importlib.util,sys
from math import lcm
from fractions import Fraction as F
from pathlib import Path
sys.dont_write_bytecode=True
family_path=Path('/snapshot/files/family.py')
if not family_path.exists():family_path=Path(__file__).with_name('family.py')
spec=importlib.util.spec_from_file_location('candidate_family',family_path)
f=importlib.util.module_from_spec(spec);spec.loader.exec_module(f)

READONCE_COVER={"version":1,"common":[4911,5008,627,3146,4426,563,1621,1992,1558,3791,4823,5042,939,98,3598,3651,4347,4648,4835,4468,4376,4667,4671,4446,7,624,2813,4083,1164,4622,4604,2005,4190,4643,372,4771,4763,3003,4472,5026,2437,4370,4328,4722,130,4741,4697,887,1494,4292,3751,2772,2153,4813,2801,1120,733,82,63,4349,5033,4942,4591,904,1998,3894,4485,4422,4869,729,4755,4653,2098,4545,1939,3515,1846,1580,1068,3081,4566,633,3931,4280,4344,4520,2812,795,2701,905,4814,3603,4213,3025,1578,4550,4410,893,3140,4884,4037,3647,1655,4570,2799,1615,158,4618,1188,2725,3249,4502,2476,738,4843,686,4989,139,4811,2705,4153,4279,1557,2988,4384,3259,240,1789],"covers":[[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[4483],[],[84],[],[1330,847,3035],[997],[],[],[],[],[],[],[],[],[],[],[],[3035,847],[1607,3037],[],[],[],[],[],[],[],[],[],[],[],[3870],[1330,1358,1761],[],[1766],[847,848],[],[199],[],[],[],[],[],[],[],[],[],[],[],[2856],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[3806],[282,3813],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[4575],[],[1529],[],[],[23,1510],[4528],[],[],[3975,857],[],[941,2204],[],[],[],[],[],[4400],[4396,1269],[],[8,1790,2478],[],[4727,4775],[],[4400],[4009,3607,3849],[2478,1404,427,1332],[],[],[],[2146],[],[397,1970,2428,1528,1358],[2787],[1495],[],[42,3493,4460,2232],[4488,3634,397,1852,858,1890,2436,4616,2003],[2003,732,4616],[],[],[210,1083,895,1772,3260,355],[855,4770,3517,2121,1282,2767,4452],[855,3997,858],[],[],[2428],[1279],[1130,387,170,1065],[],[],[],[],[],[],[],[2489,3834],[],[405,1086,3045,4821,99],[],[],[],[],[],[],[],[3790],[4400],[],[4366,4536],[4252,5035],[1441,3790],[],[2478,1529,1404],[],[1330,1800,3264,4747,4804],[387,3030,2490],[],[],[744,99,3680,42,4088],[],[2186,4094,3489,3010,2979,4440,3492,3963,3542,42,1973,4929],[99,3996,42,3988],[3514],[],[],[],[4455,3260],[],[630],[],[],[],[],[1,46,251,3312,4425],[4396,4323],[157,4237,3840,1380,3480,4323,2540,1394,3834,513,3842],[3744,189,86,3514],[4774,3523,3480],[],[],[3870],[4036,3910],[],[],[3810,4278],[153],[4278,4036,4008],[4542,1199,2482,2003,4578,3634,1726,2656,4720],[2325,4368,266],[1128,4720,1130,4332,2767,223,1282,3608],[2547,1497,4038,975,4330,2931],[5037],[],[635,4512,3045,4278,153,2704],[4278,3806,2704,296],[2736,1611,3271,2358,4330,4008,3586,1940],[],[3959],[],[2666,4352],[],[],[],[],[],[],[],[],[],[4536],[],[],[4342,4775],[],[],[],[4575],[],[4727,4352],[],[],[],[],[],[],[],[3144],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[4325,2429],[],[1510,941],[],[],[],[1441],[4251,3789],[189,2681,3498],[],[],[4391],[4117],[1097,207,387,1112],[],[],[],[],[],[],[3250,3702],[4400],[3250],[4483,4325],[2429,1790,858,3255,1501],[],[],[],[4774,3523],[],[4798,3744,4323,5022,2549],[4038],[3789],[5035],[4251,3789,4009],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[]]}

def affine_flat_indices(legal,n,d):
    legal_set=set(legal);index={v:i for i,v in enumerate(legal)};planes=[]
    def extend(offset,points,pivots,low):
        if len(pivots)==d:
            planes.append(tuple(index[v] for v in points));return
        k=len(pivots)
        for pivot in range(low,n-(d-k)+1):
            if offset>>pivot&1:continue
            free=[j for j in range(pivot) if j not in pivots]
            for bits in range(1<<len(free)):
                v=(1<<pivot)|sum(((bits>>j)&1)<<b for j,b in enumerate(free))
                more=[x^v for x in points]
                if all(x in legal_set for x in more):
                    extend(offset,points+more,pivots+[pivot],pivot+1)
    for a in legal:extend(a,[a],[],0)
    return planes


def transform_coloring(value,permutation,flip):
    return sum((((value>>i)&1)^flip)<<j for i,j in enumerate(permutation))


def canonical_colorings(group):
    remaining=set(range(4096));representatives=[]
    while remaining:
        value=min(remaining)
        orbit={transform_coloring(value,p,f) for p,f in group}
        assert value==min(orbit)
        representatives.append(value);remaining.difference_update(orbit)
    return representatives


def ramsey84_normalized_roots():
    legal=[a for a in range(256) if not(a&1 and a&4)
           and not(a&2 and a&8) and not(a&16 and a&32)]
    planes=affine_flat_indices(legal,8,3)
    assert len(legal)==108 and len(planes)==14391
    ids={};values=[]
    for a in legal:
        key=(int(not(a&1) and not(a&2)),a>>2)
        if key not in ids:ids[key]=len(ids)+1
        values.append(ids[key])
    assert len(ids)==84
    edges=sorted({tuple(sorted({values[j] for j in plane})) for plane in planes})
    assert len(edges)==11349
    base=[]
    for edge in edges:base.extend([list(edge),[-v for v in edge]])
    small=[a for a in range(16) if not(a&1 and a&2)]
    small_planes=[sum(1<<small.index(x) for x in points)
        for points in itertools.combinations(small,4)
        if points[0]^points[1]^points[2]^points[3]==0]
    bad={v for v in range(4096) if not any(v&m in (0,m) for m in small_planes)}
    assert len(small_planes)==39 and len(bad)==48
    core=238;orbit=set();stabilizer=[]
    for pair in itertools.permutations((0,1,2)):
        for r1,r2 in itertools.permutations((1,2,3),2):
            for shifts in itertools.product(range(4),repeat=3):
                permutation=[]
                for a in small:
                    b=a>>2
                    free=((b&r1).bit_count()&1)|(((b&r2).bit_count()&1)<<1)
                    free ^=shifts[a&3]
                    permutation.append(small.index(pair[a&3]|(free<<2)))
                for flip in (0,1):
                    result=transform_coloring(core,permutation,flip);orbit.add(result)
                    if result==core:stabilizer.append((permutation,flip))
    assert orbit==bad and len(stabilizer)==96
    zeros=[i for i in range(12) if not(core>>i&1)]
    ones=[i for i in range(12) if core>>i&1]
    allowed=set()
    for iz,io in itertools.product([-1]+zeros,[-1]+ones):
        code=4095^core
        for j in (iz,io):
            if j>=0:code ^=1<<j
        allowed.add(code)
    assert len(allowed)==49
    cases=[];remaining=set(allowed)
    while remaining:
        code=min(remaining)
        current={transform_coloring(code,p,f) for p,f in stabilizer}
        assert current<=allowed
        cases.append(code);remaining.difference_update(current)
    assert cases==[1809,1811,1813,3857]
    core0=[ids[(0,t<<2)] for t in small]
    core1=[ids[(1,t<<2)] for t in small]
    A=[ids[(0,1|(t<<2))] for t in small]
    B=[ids[(0,2|(t<<2))] for t in small]
    translations=[[small.index(a^(shift<<2)) for a in small] for shift in range(4)]
    canonical_B=canonical_colorings([(p,0) for p in translations])
    assert len(canonical_B)==1072
    outputs=[]
    for code in cases:
        current=[(p,f) for p,f in stabilizer if transform_coloring(code,p,f)==code]
        group=[]
        for p,f in current:
            for shift in range(4):
                group.append(([small.index(small[j]^(shift<<2)) for j in p],f))
        canonical_A=canonical_colorings(group)
        clauses=[list(q) for q in base]
        for j,v in enumerate(core0):clauses.append([v if core>>j&1 else -v])
        for j,v in enumerate(core1):clauses.append([v if code>>j&1 else -v])
        def prefix_exclusions(codes,variables):
            def visit(xs,bit,prefix):
                if not xs:
                    clauses.append([-variables[j] if b else variables[j] for j,b in prefix]);return
                if len(xs)==1<<(bit+1):return
                assert bit>=0
                for b in (0,1):
                    visit([x for x in xs if (x>>bit&1)==b],bit-1,prefix+[(bit,b)])
            visit(codes,11,[])
        prefix_exclusions(canonical_A,A);prefix_exclusions(canonical_B,B)
        previous=None;nv=84
        for j in range(11,-1,-1):
            a,b=A[j],B[j];prefix=[] if previous is None else [-previous]
            clauses.append(prefix+[-a,b])
            if j:
                nv+=1
                clauses.extend([prefix+[a,b,nv],prefix+[-a,-b,nv]])
                previous=nv
        assert nv==95
        units={};work=[tuple(q) for q in clauses]
        while True:
            new=[];found={}
            for clause in work:
                if any(abs(x) in units and units[abs(x)]==(x>0) for x in clause):continue
                clause=tuple(x for x in clause if abs(x) not in units)
                assert clause
                if len(clause)==1:
                    variable=abs(clause[0]);value=clause[0]>0
                    assert variable not in found or found[variable]==value
                    found[variable]=value
                else:new.append(clause)
            work=new
            if not found:break
            for variable,value in found.items():
                assert variable not in units or units[variable]==value
            units.update(found)
        roots=sorted(set(tuple(sorted(q,key=lambda x:(abs(x),x))) for q in work),
                     key=lambda q:(len(q),q))
        outputs.append((code,roots,len(canonical_A),len(units)))
    assert [len(x[1]) for x in outputs]==[5130,5049,4773,4408]
    assert [x[2] for x in outputs]==[270,364,262,49]
    return outputs


def resolution_refutation(roots,nv,cap=100000):
 clauses=[];watch={v:[] for v in range(-nv,nv+1)if v}
 values=[0]*(nv+1);levels=[0]*(nv+1);reason=[-1]*(nv+1)
 activity=[0.]*(nv+1);phase=[-1]*(nv+1);increment=1.
 trail=[];limits=[];head=0;conflicts=0;resolutions=0;decisions=0
 def enqueue(lit,why):
  v=abs(lit);sign=1 if lit>0 else -1
  if values[v]:return values[v]==sign
  values[v]=sign;phase[v]=sign;levels[v]=len(limits);reason[v]=why;trail.append(lit);return True
 def add(clause):
  assert all(1<=abs(l)<=nv for l in clause)
  assert len({abs(l)for l in clause})==len(clause)
  ix=len(clauses);clauses.append(list(clause))
  if len(clause)>1:
   watch[clause[0]].append(ix);watch[clause[1]].append(ix)
  return ix
 for row in roots:
  if not row:return {'sat':False,'conflicts':0,'resolutions':0,'decisions':0}
  ix=add(row)
  for l in row:activity[abs(l)]+=2.**-len(row)
  if len(row)==1:
   if not enqueue(row[0],ix):
    # Two contradictory root units resolve to the empty clause.
    assert reason[abs(row[0])]>=0 and clauses[reason[abs(row[0])]]==[-row[0]]
    return {'sat':False,'conflicts':1,'resolutions':1,'decisions':0}
 def propagate():
  nonlocal head
  while head<len(trail):
   false=-trail[head];head+=1;ws=watch[false];i=0
   while i<len(ws):
    ix=ws[i];row=clauses[ix]
    if row[0]==false:row[0],row[1]=row[1],row[0]
    assert row[1]==false
    other=row[0];value=values[abs(other)]*(1 if other>0 else -1)
    if value==1:i+=1;continue
    for j in range(2,len(row)):
     lit=row[j]
     if values[abs(lit)]*(1 if lit>0 else -1)!=-1:break
    else:
     if value==-1:return ix
     assert enqueue(other,ix);i+=1;continue
    row[1],row[j]=row[j],row[1];ws[i]=ws[-1];ws.pop();watch[row[1]].append(ix)
  return -1
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
   learnt=set(clauses[bad]);dl=len(limits)
   assert all(values[abs(l)]*(1 if l>0 else -1)==-1 for l in learnt)
   def resolve(lit):
    nonlocal resolutions,learnt
    ri=reason[abs(lit)];assert ri>=0 and lit in clauses[ri] and -lit in learnt
    # This exact set equality is the resolution inference rule.
    learnt=(learnt-{-lit})|(set(clauses[ri])-{lit})
    assert not any(-x in learnt for x in learnt)
    assert all(values[abs(x)]*(1 if x>0 else -1)==-1 for x in learnt)
    resolutions+=1
   if dl==0:
    for lit in reversed(trail):
     if -lit in learnt:resolve(lit)
    assert not learnt
    return {'sat':False,'conflicts':conflicts,'resolutions':resolutions,'decisions':decisions,'clauses':len(clauses)}
   cursor=len(trail)-1
   while sum(levels[abs(l)]==dl for l in learnt)>1:
    while -trail[cursor]not in learnt:cursor-=1
    lit=trail[cursor];cursor-=1;resolve(lit)
   current=[l for l in learnt if levels[abs(l)]==dl];assert len(current)==1
   first=current[0];others=sorted(learnt-{first},key=lambda l:(-levels[abs(l)],abs(l),l))
   jump=max((levels[abs(l)]for l in others),default=0)
   assert jump<dl
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
    assert all(any(values[abs(l)]*(1 if l>0 else -1)==1 for l in row)for row in roots)
    return {'sat':True,'model':values[1:],'conflicts':conflicts,'resolutions':resolutions,'decisions':decisions}
   variable=max(free,key=lambda v:(activity[v],-v));limits.append(len(trail));decisions+=1
   assert enqueue(variable*phase[variable],-1)

def readonce_seven_audit(envelope):
    started=time.monotonic()
    legal=[a for a in range(128) if not(a&1 and a&2) and not(a&4 and a&8)]
    planes=affine_flat_indices(legal,7,3)
    assert len(legal)==72 and len(planes)==5049
    def trees(labels):
        if len(labels)==1:return [labels[0]]
        out=[];first=labels[0];rest=labels[1:]
        for k in range(len(rest)):
            for more in itertools.combinations(rest,k):
                left=(first,)+more;right=tuple(x for x in labels if x not in left)
                for a in trees(left):
                    for b in trees(right):out.append((a,b))
        return out
    def serial(tree,permutation):
        if isinstance(tree,int):return str(permutation[tree])
        a,b=serial(tree[0],permutation),serial(tree[1],permutation)
        return '('+','.join(sorted((a,b)))+')'
    permutations=[]
    for swap,a,b in itertools.product((0,1),repeat=3):
        for free in itertools.permutations((4,5,6)):
            p=[0,1,2,3]+list(free)
            if a:p[0],p[1]=p[1],p[0]
            if b:p[2],p[3]=p[3],p[2]
            if swap:p[:4]=p[2:4]+p[:2]
            permutations.append(p)
    all_trees=trees(tuple(range(7)));representatives=[];seen=set()
    for tree in all_trees:
        key=min(serial(tree,p) for p in permutations)
        if key not in seen:seen.add(key);representatives.append(tree)
    assert len(permutations)==48 and len(all_trees)==10395 and len(representatives)==350
    assert set(envelope)=={'version','common','covers'} and envelope['version']==1
    common=envelope['common'];covers=envelope['covers']
    assert len(common)==128 and len(set(common))==128 and len(covers)==350
    assert all(type(i)is int and 0<=i<len(planes)for i in common)
    for row in covers:
        assert isinstance(row,list)and all(type(i)is int and 0<=i<len(planes)for i in row)
    assignments=16*3**6;full=(1<<assignments)-1;repeat=full//65535
    signs=[sum(((v>>i)&1)<<v for v in range(16))*repeat for i in range(4)]
    op_masks=[]
    for i in range(6):
        stride=16*3**i
        a=full//((1<<(3*stride))-1)*((1<<stride)-1)
        op_masks.append((a,a<<stride,a<<(2*stride)))
    for tree,ids in zip(representatives,covers):
        operations=[]
        def compile_tree(t):
            if isinstance(t,int):return t
            a,b=compile_tree(t[0]),compile_tree(t[1])
            operations.append((a,b));return 6+len(operations)
        root=compile_tree(tree);assert root==12 and len(operations)==6
        values={}
        for j in {j for i in common+ids for j in planes[i]}:
            a=legal[j]
            v=[signs[i] if not(a>>i&1) else full^signs[i] for i in range(4)]
            v.extend(full if a>>i&1 else 0 for i in (4,5,6))
            for (left,right),(aa,oo,xx) in zip(operations,op_masks):
                x,y=v[left],v[right]
                v.append(((x&y)&aa)|((x|y)&oo)|((x^y)&xx))
            values[j]=v[root]
        covered=0
        for i in common+ids:
            if covered==full:break
            zero=one=full
            for j in planes[i]:one &=values[j];zero &=full^values[j]
            covered |=zero|one
        assert covered==full
    return {'passed':True,'shape_orbits':350,'assignments_each':assignments,
            'assignments_verified':350*assignments,
            'plane_witnesses':sum(map(len,covers)),'shared_planes':len(common),
            'elapsed_seconds':time.monotonic()-started}



def ramsey_direct_audit():
 started=time.monotonic();cases=[]
 for code,roots,canonical_A,units in ramsey84_normalized_roots():
  result=resolution_refutation(roots,95)
  assert result['sat']is False
  cases.append(dict(code=code,roots=len(roots),canonical_A=canonical_A,units=units,**result))
 assert len(cases)==4
 return dict(passed=True,cases=cases,elapsed_seconds=time.monotonic()-started)

def run():
 result={'ramsey84_resolution':ramsey_direct_audit(),'readonce_seven':readonce_seven_audit(READONCE_COVER)}
 # New exact convex-midpoint potential certificate; calculus is supplied in proof.md.
 a=F(550791690219,10**12);X=F(104292796651,10**11);wm=X+F(1,10**12);w0=F(1683067911,500000000)
 pc=a/4-F(102,10**8)*w0-F(1,10**13);h=pc*X/F(2000001,10**6)-F(1,10**13)
 assert f.analytic_prefilter((h,a,wm,w0,X,pc))
 assert w0+X==F(440906378851,10**11)<F(8818127579,2000000000)
 M=65536;step=(w0-wm)/M;scale=10**15
 V=lcm(wm.denominator,(step/2).denominator);U=int((wm+step/2)*V);inc=int(step*V)
 assert F(U,V)==wm+step/2 and F(inc,V)==step
 ap,aq=a.numerator,a.denominator;total=0;cross=0
 for j in range(M):
  z=U*aq-ap*V
  num=scale*100000*U*z*z
  den=V*(100000*aq*aq*U*U+300000*ap*aq*U*V+4*ap*V*z)
  assert den>0
  integer=num//den;assert integer*den<=num<(integer+1)*den
  if j in (0,1,M//2,M-1):
   w=wm+(F(2*j+1,2))*step;beta=a/w;alpha=1/(1-beta)+F(1,100000)
   value=scale*w*(1-beta)/(1+4*alpha*beta)
   assert integer==value.numerator//value.denominator;cross+=1
  total+=integer;U+=inc
 assert total==49109014265147760617
 ratio=8*h*step*F(total,scale)
 assert ratio==F(40960020509752098178821891163443417586661219884120373286607,40960020480000000000000000000000000000000000000000000000000)>1
 W=w0-wm;mu=ratio-1-8*h*W/10**12;assert mu>0
 target=176*h*W*W/mu;J=1+-(-target.numerator//target.denominator)
 assert J>0 and ratio-8*h*W/10**12-88*h*W*W/J>1+mu/2>1
 result['convex_midpoint_certificate_terms']=M
 result['convex_midpoint_certificate_sum']=str(total)
 result['convex_midpoint_rational_crosschecks']=cross
 result['original_left_mesh_existence_integer']=str(J)
 result['all_passed']=True
 result['scope']='Exact four normalized Ramsey CNF refutations by checked resolution, complete read-once seven-variable cover and convex-midpoint potential witness only; no full accepted tree or universal theorem execution.'
 return result

if __name__=='__main__':
 request=json.loads(input())
 if request!={}:raise ValueError('empty object required')
 print(json.dumps(run(),sort_keys=True))
