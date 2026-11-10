pass;import json,itertools,time,importlib.util,sys;from math import lcm
from fractions import Fraction as F;from pathlib import Path;sys.dont_write_bytecode=True;family_path=Path('/snapshot/files/family.py')
if not family_path.exists():family_path=Path(__file__).with_name('family.py')
spec=importlib.util.spec_from_file_location('candidate_family',family_path);f=importlib.util.module_from_spec(spec);spec.loader.exec_module(f)
from functools import lru_cache
COMMON={"near_cross":[4436,617,1,5040,4323,3680,3697,5023,3710,147,4763,3716,5022,4863,629,4797,619,4971,603,2701,4825,3249,3477,0,4736,4762,2013,4418,3651,5030,4406,3476,5008,4761,3629,5028,2097,4447,4964,2975,4843,3,2652,3815,3647,4787,631,3565,3607,4428,621,3475,4556,5043,3835,4320,3073,3474,3603,3649,3259,4775,4130,4773],"near_mixed":[2652,4801,4856,5040,389,4897,2813,4618,4436,4406,7,4326,2674,2655,1210,3024,4323,1066,4515,3,4433,2680,1735,2675,4520,1720,4328,42,45,4176,1,729,568,4797,4517,4332,3514,2656,1365,4964,732,5001,4900,5028,4564,397,2658,1098,3522,1613,4899,0,15,562,1093,4091,2725,2673,387,4771,1088,4769,4265,4518],"readonce8":[1504,2094,1581,1295,1323,244,1927,383,230,987,981,1078,2079,1378,233,660,2052,372,586,1381,1819,355,1377,136,147,1006,865,2100,7,1379,696,1859,183,868,207,923,1324,1341,2031,1353,156,989,1401,674,870,1038,2086,850,109,353,1345,1428,632,1207,127,284,864,1326,335,188,945,139,135,2090],"reverse_A":[266,515,269,732,1641,518,1890,2652,1188,735,4485,4545,4662,4692,1185,1437,3024,3249,872,2722,729,1164,2658,4542,2322,4689,4482,1644,2673,1992,4659,74,4351,1995,3045,2725,1893,1434,12,3484,245,4911,1167,4881,13,2701,4722,95,4365,2655,782,4775,4767,251,4986,2325,4971,45,4352,4328,76,544,302,46],"reverse_C":[732,1185,4482,3894,735,2132,3059,8,2592,3414,1188,4485,872,1437,518,1995,1164,4545,2652,71,2465,5024,4328,729,4662,2135,4773,869,4613,4692,4814,5028,248,4398,1434,3020,2095,4083,1988,2425,4326,3522,3988,5036,1551,3520,4937,3954,4801,5034,1495,4927,3480,3484,12,4659,3354,4986,4911,4401,1307,1167,3284,3255],"extension":[6,8,4,1,0,82,11,9,10,45,49,2,37,85,158,28,50,32,18,51,61,15,53,294,14,52,81,25,42,270,33,24,27,16,278,93,272,273,23,196,17,13,87,121,59,281,26,54,195,40,174,55,66,323,43,186,22,316,114,21,83,75,7,73],"readonce7":[4911,5008,627,3146,4426,563,1621,1992,1558,3791,4823,5042,939,98,3598,3651,4347,4648,4835,4468,4376,4667,4671,4446,7,624,2813,4083,1164,4622,4604,2005,4190,4643,372,4771,4763,3003,4472,5026,2437,4370,4328,4722,130,4741,4697,887,1494,4292,3751,2772,2153,4813,2801,1120,733,82,63,4349,5033,4942,4591,904,1998,3894,4485,4422,4869,729,4755,4653,2098,4545,1939,3515,1846,1580,1068,3081,4566,633,3931,4280,4344,4520,2812,795,2701,905,4814,3603,4213,3025,1578,4550,4410,893,3140,4884,4037,3647,1655,4570,2799,1615,158,4618,1188,2725,3249,4502,2476,738,4843,686,4989,139,4811,2705,4153,4279,1557,2988,4384,3259,240,1789]}
EXTENSION_FLATS=[[0,6,18,64,256],[0,6,18,64,512],[0,6,18,66,512],[1,3,28,76,512],[0,6,18,128,256],[0,1,4,8,16],[0,6,18,130,256],[1,3,28,64,256],[0,6,18,64,514],[0,6,18,66,514],[0,6,18,128,258],[0,6,18,130,258],[1,3,28,64,512],[1,3,28,128,268],[1,3,28,128,256],[0,6,18,64,552],[0,9,17,65,513],[0,9,17,64,513],[0,9,17,65,512],[0,18,40,68,512],[0,18,40,68,516],[0,18,40,64,512],[0,9,17,64,516],[0,9,17,64,548],[0,9,17,64,512],[0,6,18,64,520],[0,6,18,66,520],[0,6,34,66,512],[0,6,34,64,512],[0,6,18,128,512],[128,1,4,8,16],[0,6,18,128,552],[1,4,24,64,512],[1,4,24,72,512],[1,4,24,72,520],[0,9,17,128,516],[0,9,17,128,548],[1,4,24,64,544],[2,8,20,68,512],[2,8,20,68,516],[16,6,32,64,514],[16,6,32,66,514],[16,6,32,64,512],[16,6,32,66,512],[1,3,28,128,512],[0,6,18,168,256],[0,18,40,128,512],[0,6,18,168,552],[1,4,24,136,264],[1,4,24,128,256],[1,4,24,128,264],[0,9,17,129,256],[0,9,17,129,257],[0,9,17,128,257],[0,18,40,128,256],[0,9,17,132,256],[0,18,40,128,260],[0,18,40,132,260],[0,9,17,132,516],[0,9,17,164,256],[0,17,36,128,512],[1,4,24,160,256],[0,9,17,164,548],[1,4,24,160,544],[2,8,20,132,260],[2,8,20,128,260],[0,9,17,128,256],[0,6,18,136,256],[0,6,18,136,258],[16,6,32,130,256],[16,6,32,130,258],[16,6,32,128,256],[16,6,32,128,258],[0,6,34,128,256],[0,6,34,128,258],[0,6,18,64,544],[0,6,18,66,544],[64,6,18,192,256],[0,6,18,128,514],[0,6,18,130,514],[0,6,18,104,512],[0,1,4,24,128],[1,4,24,64,256],[0,1,4,24,136],[0,1,4,24,512],[0,18,40,64,256],[0,1,4,24,520],[0,6,18,128,296],[0,6,18,66,552],[0,1,4,24,544],[0,1,4,40,512],[0,1,4,24,640],[1,4,24,72,547],[0,6,18,168,258],[1,4,24,163,264],[1,4,24,128,544],[0,6,18,104,514],[0,1,4,128,512],[0,1,4,24,160],[0,1,4,40,128],[1,3,12,128,512],[0,4,24,136,257],[0,4,18,40,640],[1,3,28,96,512],[0,4,24,65,520],[0,4,18,40,512],[1,3,28,128,288],[1,9,17,128,257],[0,12,20,65,512],[0,4,18,40,576],[0,9,17,65,516],[1,9,17,65,516],[1,9,17,65,548],[0,4,24,65,544],[0,9,17,65,548],[0,4,24,65,512],[0,4,18,40,514],[1,3,28,96,544],[0,4,18,40,520],[0,4,40,65,512],[0,9,33,65,512],[0,4,18,40,128],[1,4,24,65,520],[1,4,24,65,512],[1,4,24,65,544],[0,4,18,130,512],[1,4,24,136,257],[1,4,24,128,257],[1,12,20,65,512],[0,12,20,128,257],[1,12,20,128,257],[0,4,18,40,384],[1,23,40,136,264],[1,23,40,128,264],[4,1,24,64,256],[0,9,17,132,257],[1,9,17,132,257],[0,9,17,164,257],[1,4,24,160,257],[1,9,17,164,257],[1,3,28,160,288],[0,4,18,40,130],[0,4,24,128,257],[0,4,18,160,512],[0,4,24,160,257],[0,4,18,128,512],[0,4,18,40,136],[0,4,24,34,640],[0,4,34,128,512],[0,9,33,128,257],[0,9,17,65,544],[0,4,9,17,128],[0,1,24,68,520],[0,1,24,68,512],[0,1,24,68,544],[0,1,40,68,512],[0,1,24,128,260],[0,4,9,17,129],[0,9,17,64,256],[0,4,9,17,512],[0,1,24,136,260],[0,4,9,17,640],[0,9,17,68,548],[0,9,17,68,512],[0,4,9,33,512],[0,9,17,68,513],[0,9,17,129,260],[0,1,24,128,512],[0,1,24,132,260],[0,1,24,132,516],[0,9,17,164,260],[0,9,17,128,260],[0,4,17,128,512],[0,4,9,33,128],[0,2,8,16,32],[4,18,40,128,260],[9,4,24,64,256],[6,8,20,160,259],[0,9,17,132,260],[0,18,40,68,514],[6,8,20,132,259],[6,8,20,128,259],[0,18,40,130,260],[0,9,33,68,512],[0,9,17,68,544],[0,9,17,68,516],[0,2,8,20,128],[1,4,8,128,512],[0,2,8,20,132],[0,9,17,100,512],[0,2,8,20,512],[0,2,8,20,516],[0,9,17,128,292],[1,4,8,128,528],[0,2,8,20,544],[0,17,36,64,256],[2,8,20,64,256],[0,2,8,36,512],[0,2,8,20,640],[0,2,8,128,512],[0,2,8,20,160],[1,4,24,96,512],[1,4,24,72,544],[0,9,17,100,513],[1,4,24,128,288],[0,2,8,132,512],[1,4,24,160,264],[0,2,8,20,384],[1,4,24,96,544],[0,2,8,36,128],[0,2,12,20,128],[0,2,12,20,132],[2,6,18,128,258],[0,18,40,66,512],[0,18,40,66,516],[0,8,20,66,516],[0,8,17,36,512],[0,12,20,66,512],[0,8,17,36,576],[0,8,17,36,640],[0,2,24,128,512],[1,4,8,144,544],[0,2,12,20,512],[0,2,12,20,516],[1,4,8,144,512],[0,2,12,20,640],[2,6,18,168,258],[0,2,12,36,512],[0,8,17,128,512],[0,8,17,128,513],[2,18,40,66,512],[2,18,40,66,516],[0,8,17,132,513],[0,18,40,128,258],[2,12,20,66,512],[0,18,40,132,258],[2,18,40,128,258],[2,18,40,132,258],[0,8,20,132,258],[0,8,17,36,128],[2,8,20,132,258],[2,8,20,128,258],[0,12,20,128,258],[2,12,20,128,258],[0,8,17,36,384],[1,4,24,160,288],[8,2,20,64,256],[0,2,20,128,512],[0,2,12,36,128],[0,6,18,72,512],[0,2,20,72,516],[0,6,18,128,264],[0,6,18,72,552],[0,2,20,72,544],[0,2,20,72,512],[0,6,18,72,514],[0,2,20,72,520],[0,2,20,132,264],[0,2,36,72,512],[0,2,20,128,264],[0,6,18,168,264],[0,6,18,130,264],[0,18,40,72,512],[1,4,8,144,528],[0,2,20,160,264],[0,6,34,72,512],[0,18,40,128,264],[0,2,20,136,264],[0,6,18,72,520],[0,2,20,136,520],[0,6,16,32,128],[1,3,16,32,128],[0,6,16,32,130],[0,6,16,32,512],[0,6,16,32,520],[0,9,16,32,512],[0,6,16,32,514],[0,9,16,32,516],[16,6,32,64,256],[0,6,16,32,640],[0,16,32,128,512],[0,6,16,32,136],[1,3,16,32,136],[0,9,16,32,128],[0,9,16,32,132],[1,3,16,32,132],[16,6,32,128,512],[80,6,32,192,256],[16,6,32,128,514],[80,6,32,194,256],[16,6,32,130,514],[0,6,16,32,384],[0,6,16,34,128],[0,6,16,34,130],[0,6,34,64,256],[0,6,32,80,514],[0,6,32,80,512],[0,6,34,80,512],[1,3,32,80,520],[0,6,16,34,512],[0,6,16,34,640],[1,3,32,80,516],[1,3,32,136,272],[1,3,32,132,272],[0,6,16,34,136],[16,6,32,80,514],[16,6,32,80,512],[0,6,32,130,272],[0,6,32,128,272],[16,6,34,80,512],[0,6,34,128,272],[16,6,34,128,272],[0,6,32,128,512],[32,6,16,64,256],[16,6,34,128,512],[0,6,16,128,512],[0,6,18,160,256],[0,6,18,160,258],[1,3,28,64,544],[65,3,28,192,544],[0,6,18,128,544],[1,3,28,128,544],[64,6,18,192,544],[1,3,28,160,256],[128,1,4,8,272],[1,3,28,160,544],[0,6,18,160,544],[17,3,32,128,512],[0,1,4,8,128],[16,1,4,8,32],[0,1,4,8,528],[1,3,28,76,544],[0,1,4,8,512],[0,6,18,128,288],[0,1,4,8,640],[0,6,18,96,512],[0,1,4,8,144],[1,3,28,160,268],[1,4,24,96,520],[0,9,17,96,512],[0,6,18,96,514],[1,4,24,136,288],[9,4,24,136,288],[9,4,24,128,288],[2,6,18,128,288],[1,4,24,129,288],[0,9,17,129,288],[0,9,17,128,288],[0,6,18,130,288],[0,6,18,96,544],[0,1,4,8,544],[1,3,12,128,528],[1,3,12,144,544],[0,4,18,160,514],[1,9,17,160,257],[0,9,17,160,257],[1,3,28,68,544],[0,1,24,160,260],[1,3,28,160,260],[0,4,9,128,512]]
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
def ramsey6_roots():
 pass;flats6=affine_flat_indices(list(range(64)),6,3);assert len(flats6)==11160;roots=[]
 for flat in flats6:
  for omit in range(8):
   clause=tuple(v+1 for j,v in enumerate(flat)if j!=omit);roots.extend((clause,tuple(-v for v in clause)))
 roots.extend((-(a+1),)if a<6 else(a+1,)for a in range(8));roots.extend(((-9,),(-17,),(-33,)));assert len(roots)==178571;yield 'maximum-six',64,roots;del roots
 flats5=affine_flat_indices(list(range(32)),5,3);assert len(flats5)==620;assignment_count=1<<16;full=(1<<assignment_count)-1;masks=[]
 for j in range(16):
  stride=1<<j
  masks.append((((1<<stride)-1)<<stride)
  *(full//((1<<(2*stride))-1)))
 case_zero_sets=({0},{0,1},{0,1,2},{0,1,2,4});expected_patterns=(3136,5760,6232,5296);expected_groups=(168,48,24,24);expected_orbits=(6,18,28,30)
 for case,zeros in enumerate(case_zero_sets,1):
  core=[0]*7+[1]+[int(j not in zeros)for j in range(8)];allowed=full
  for flat in flats5:
   ones=zeroes=full
   for j in flat:
    if j<16:
     ones&=full if core[j]else 0;zeroes&=full if not core[j]else 0
    else:
     ones&=masks[j-16];zeroes&=full^masks[j-16]
   allowed&=full^(ones|zeroes)
  patterns=set()
  while allowed:
   bit=(allowed&-allowed).bit_length()-1;allowed&=allowed-1;patterns.add(bit)
  assert len(patterns)==expected_patterns[case-1];group=[]
  for columns in itertools.permutations(range(1,8),3):
   linear=[0]*8
   for v in range(8):
    for j in range(3):
     if v>>j&1:linear[v]^=columns[j]
   if len(set(linear))!=8:continue
   shift=7^linear[7]
   for shear in range(8):
    permutation=[(linear[v&7]^shift^(shear if v&8 else 0))
    |(v&8)for v in range(16)]
    if all(core[permutation[v]]==core[v]for v in range(16)):
     group.append(permutation)
  assert len(group)==expected_groups[case-1];seen=set();representatives=[];orbit_sizes=[]
  for pattern in sorted(patterns):
   if pattern in seen:continue
   orbit=set()
   for permutation in group:
    image=transform_coloring(pattern,permutation,0)
    for shift in range(16):
     orbit.add(sum(((image>>j)&1)<<(j^shift)
     for j in range(16)))
   assert orbit<=patterns;seen|=orbit;orbit_sizes.append(len(orbit));representatives.append(min(v for v in orbit if not(v&1)))
  assert seen==patterns and len(representatives)==expected_orbits[case-1];assert sum(orbit_sizes)==len(patterns)
  for index,pattern in enumerate(representatives):
   fixed=core+[(pattern>>j)&1 for j in range(16)];roots={(-1,)}
   for flat in flats6:
    for color in(0,1):
     if any(j<32 and fixed[j]!=color for j in flat):continue
     roots.add(tuple((1 if color==0 else-1)*(j-31)
     for j in flat if j>=32))
   roots=sorted(roots,key=lambda clause:(len(clause),clause));assert roots and all(clause for clause in roots);yield 'core-seven-%d-%d'%(case,index),32,roots
@lru_cache(maxsize=None)
def binary_trees(labels):
 if len(labels)==1:return(labels[0],)
 output=[];first=labels[0];rest=labels[1:]
 for size in range(len(rest)):
  for selected in itertools.combinations(rest,size):
   left=(first,)+selected;right=tuple(v for v in labels if v not in left)
   for a in binary_trees(left):
    for b in binary_trees(right):output.append((a,b))
 return tuple(output)
@lru_cache(maxsize=100000)
def tree_string(tree,permutation):
 if isinstance(tree,int):return str(permutation[tree])
 return '('+','.join(sorted((tree_string(tree[0],permutation),
 tree_string(tree[1],permutation))))+')'
def clear_tree_caches():
 binary_trees.cache_clear();tree_string.cache_clear()
def repeated_component_shapes(kind):
 permutations=[]
 if kind=='cross':
  for swap in(0,1):
   for free in itertools.permutations((3,4,5)):
    p=[0,1,2,*free,6,6]
    if swap:p[1],p[2]=p[2],p[1]
    permutations.append(tuple(p))
  forbidden=(1,2)
 else:
  assert kind=='mixed';forbidden=(1,2,3)
  for a,b in itertools.product((0,1),repeat=2):
   p=list(range(7))+[6]
   if a:p[2],p[3]=p[3],p[2]
   if b:p[4],p[5]=p[5],p[4]
   permutations.append(tuple(p))
 seen=set();shapes=[]
 for size in range(1,6):
  for subset in itertools.combinations(range(6),size):
   if size==1 and subset[0]in forbidden:continue
   images=[tuple(sorted(p[v]for v in subset))for p in permutations]
   if subset!=min(images):continue
   stabilizer=[p for p,image in zip(permutations,images)if image==subset];outside=tuple(sorted(tuple(i for i in range(6)if i not in subset)+(6,7)))
   for inside in binary_trees(subset):
    for outer in binary_trees(outside):
     key=min(tree_string(inside,p)+'|'+tree_string(outer,p)
     for p in stabilizer)
     if key in seen:continue
     seen.add(key);shapes.append((inside,outer,size))
 if kind=='mixed':
  for outer in binary_trees(tuple(range(7))):
   key='extra|'+min(tree_string(outer,p)for p in permutations)
   if key not in seen:seen.add(key);shapes.append((-1,outer,0))
 assert len(shapes)==(4208 if kind=='cross' else 12468);clear_tree_caches();return shapes
def near_readonce_audit(envelope,kind):
 started=time.monotonic();legal=[a for a in range(128)if not(a&1 and a&4)and not(a&2 and a&8)];flats=affine_flat_indices(legal,7,3);assert len(flats)==5049
 if kind=='mixed':
  permutation=(1,3,2,4,0,5,6);legal=[sum(((a>>j)&1)<<permutation[j]for j in range(7))for a in legal]
 common=COMMON['near_'+kind];shapes=repeated_component_shapes(kind);prepared={};base_arrays={};bottom_cache={};previous=None;cached=None
 total_assignments=total_hints=total_steps=0
 for inside,outside,size in shapes:
  signbits=3 if kind=='cross' else(5 if size==0 else 4)
  if signbits not in prepared:
   signs=1<<signbits;count=signs*3**6;full=(1<<count)-1
   sign_masks=[sum(((v>>i)&1)<<v for v in range(signs))
   *(full//((1<<signs)-1))for i in range(signbits)]
   operation_masks=[]
   for i in range(6):
    stride=signs*3**i;mask=full//((1<<(3*stride))-1)*((1<<stride)-1)
    operation_masks.append((mask|(mask<<stride),
    (mask<<stride)|(mask<<(2*stride))))
   prepared[signbits]=(count,full,sign_masks,operation_masks)
  count,full,sm,om=prepared[signbits];total_assignments+=count;operations=[]
  def compile_tree(tree):
   if isinstance(tree,int):return tree
   a=compile_tree(tree[0]);b=compile_tree(tree[1]);operations.append((a,b));return 7+len(operations)
  inner=compile_tree(inside)if size else-1;cut=len(operations);root=compile_tree(outside);assert len(operations)==6
  if signbits not in base_arrays:
   rows=[]
   for point in legal:
    if kind=='cross':
     values=[full if not(point&1)and not(point&2)else 0];values+=[full if point>>i&1 else 0 for i in(2,3,4,5,6)];values[1]^=sm[0];values[2]^=sm[1];values+=[0,0]
    else:
     g=full if point&1 and not(point&2)else 0
     if not size:g=(sm[4]^(full if point&2 else 0))if point&1 else 0
     values=[g]+[full if point>>i&1 else 0 for i in(2,3,4,5,6)]
     for j in range(3):values[j+1]^=sm[j]
     values+=[0,0]
     if not size:values[6]=(full if point&1 else 0)^sm[3]
    rows.append(values)
   base_arrays[signbits]=[list(column)for column in zip(*rows)]
  def node_value(values,a,b,index):
   key=(signbits,index,min(a,b),max(a,b));limit=7 if not size else 6;simple=a<limit and b<limit
   if simple and key in bottom_cache:return bottom_cache[key]
   ma,mb=om[index];output=[((x&y)&ma)^((x^y)&mb)for x,y in zip(values[a],values[b])]
   if simple:bottom_cache[key]=output
   return output
  if previous!=(inside,signbits):
   cached=base_arrays[signbits][:]
   for index,(a,b)in enumerate(operations[:cut]):
    cached.append(node_value(cached,a,b,index))
   previous=(inside,signbits)
  values=cached[:]
  if size:
   values[6]=values[inner];values[7]=[v^sm[2 if kind=='cross' else 3]for v in values[inner]]
  for index,(a,b)in enumerate(operations[cut:],cut):
   values.append(node_value(values,a,b,index))
  hints=[];total_hints+=len(hints);total_steps+=check_union_cover(values[root],flats,common,hints,full)
 pass
 return{'passed':True,'shapes':len(shapes),'assignments':total_assignments,
 'hinted_flats':total_hints,'cover_steps':total_steps,
 'elapsed_seconds':time.monotonic()-started}
def readonce_eight_audit(envelope):
 started=time.monotonic()
 legal=[a for a in range(256)if not(a&1 and a&4)and not(a&2 and a&8)
 and not(a&16 and a&32)and not(a&64 and a&128)]
 flats=affine_flat_indices(legal,8,3);assert len(flats)==2106;common=COMMON['readonce8'];permutations=[]
 for a,b,e,swap in itertools.product((0,1),repeat=4):
  p=list(range(7))
  if a:p[1],p[2]=p[2],p[1]
  if b:p[3],p[4]=p[4],p[3]
  if e:p[5],p[6]=p[6],p[5]
  if swap:p[3:7]=p[5:7]+p[3:5]
  permutations.append(tuple(p))
 seen=set();representatives=[]
 for tree in binary_trees(tuple(range(7))):
  key=min(tree_string(tree,p)for p in permutations)
  if key not in seen:seen.add(key);representatives.append(tree)
 assert len(representatives)==880;clear_tree_caches();count=64*3**6;full=(1<<count)-1
 signs=[sum(((v>>i)&1)<<v for v in range(64))*(full//((1<<64)-1))
 for i in range(6)]
 masks=[]
 for i in range(6):
  stride=64*3**i;mask=full//((1<<(3*stride))-1)*((1<<stride)-1);masks.append((mask|(mask<<stride),(mask<<stride)|(mask<<(2*stride))))
 rows=[]
 for point in legal:
  row=[full if not(point&1)and not(point&2)else 0];row+=[(full if point>>(i+2)&1 else 0)^signs[i]for i in range(6)];rows.append(row)
 base=[list(column)for column in zip(*rows)];bottom={};total_hints=total_steps=0
 for tree in representatives:
  operations=[]
  def compile_tree(tree):
   if isinstance(tree,int):return tree
   a=compile_tree(tree[0]);b=compile_tree(tree[1]);operations.append((a,b));return 6+len(operations)
  root=compile_tree(tree);assert len(operations)==6;values=base[:]
  for index,(a,b)in enumerate(operations):
   key=(index,min(a,b),max(a,b));simple=a<7 and b<7
   if simple and key in bottom:result=bottom[key]
   else:
    ma,mb=masks[index];result=[((x&y)&ma)^((x^y)&mb)for x,y in zip(values[a],values[b])]
    if simple:bottom[key]=result
   values.append(result)
  hints=[];total_hints+=len(hints);total_steps+=check_common_cover(values[root],flats,common,hints,full)
 pass
 return{'passed':True,'shapes':880,'assignments':880*count,
 'hinted_flats':total_hints,'cover_steps':total_steps,
 'elapsed_seconds':time.monotonic()-started}
def graft_new_leaf(tree,label):
 yield('n',tree,label)
 if isinstance(tree,int):return
 if tree[0]=='o':
  _,table,a,b=tree
  for value in graft_new_leaf(a,label):yield('o',table,value,b)
  for value in graft_new_leaf(b,label):yield('o',table,a,value)
 else:
  _,a,b=tree
  for value in graft_new_leaf(a,label):yield('n',value,b)
  for value in graft_new_leaf(b,label):yield('n',a,value)
def tagged_tree_string(tree,permutation):
 if isinstance(tree,int):return str(permutation[tree])
 label='o'+str(tree[1])if tree[0]=='o' else 'n'
 return label+'('+','.join(sorted((tagged_tree_string(tree[-2],permutation),
 tagged_tree_string(tree[-1],permutation))))+')'
def extension_shapes():
 a=('o',14,0,2);b=('o',14,1,3);c=('o',14,4,5);originals=(('o',6,a,('o',6,b,c)),('o',6,('o',6,a,b),c));group=[]
 for swap,free,partners in itertools.product((0,1),repeat=3):
  p=list(range(9))
  if swap:p[0],p[1]=p[1],p[0];p[2],p[3]=p[3],p[2]
  if free:p[4],p[5]=p[5],p[4]
  if partners:p[7],p[8]=p[8],p[7]
  group.append(tuple(p))
 seen=set();output=[]
 for tree in originals:
  for first in graft_new_leaf(tree,6):
   for second in graft_new_leaf(first,7):
    for third in graft_new_leaf(second,8):
     key=min(tagged_tree_string(third,p)for p in group)
     if key not in seen:seen.add(key);output.append(third)
 assert len(output)==1450;return output
def extension_values(legal):
 prepared={};bases={};bottom={}
 for tree in extension_shapes():
  operations=[]
  def compile_tree(tree):
   if isinstance(tree,int):return tree,tree<6
   a,old_a=compile_tree(tree[-2]);b,old_b=compile_tree(tree[-1])
   if tree[0]=='o':table=tree[1];options=None
   else:
    assert not old_b;table=None;options=(4,6,14)if old_a else(2,4,6,8,14)
   operations.append((a,b,table,options));return 8+len(operations),old_a or old_b
  root,_=compile_tree(tree);assert len(operations)==8;choices=tuple(options for _,_,_,options in operations if options);assert len(choices)==3
  if choices not in prepared:
   count=4
   for options in choices:count*=len(options)
   full=(1<<count)-1;stride=4;coefficients=[]
   for options in choices:
    k=len(options);mask=full//((1<<(k*stride))-1)*((1<<stride)-1);blocks=[mask<<(j*stride)for j in range(k)];m1=m2=m3=0
    for table,block in zip(options,blocks):
     if table&4:m1|=block
     if table&2:m2|=block
     if((table>>1)^(table>>2)^(table>>3))&1:m3|=block
    coefficients.append((m1,m2,m3));stride*=k
   prepared[choices]=(count,full,coefficients)
  count,full,coefficients=prepared[choices]
  if count not in bases:
   repeat=full//15;rows=[]
   for point in legal:
    row=[full if point>>j&1 else 0 for j in range(6)];alpha=0 if point&64 else 1;beta=0 if point&256 else 1;row.append((1<<(2*alpha+beta))*repeat)
    row.extend((full if point&128 else 0,full if point&512 else 0));rows.append(row)
   bases[count]=[list(column)for column in zip(*rows)]
  values=bases[count][:];new_index=0
  for a,b,table,options in operations:
   if table is not None:
    assert table in(6,14);key=(count,table,min(a,b),max(a,b));simple=a<6 and b<6
    if simple and key in bottom:result=bottom[key]
    elif table==6:result=[x^y for x,y in zip(values[a],values[b])]
    else:result=[x|y for x,y in zip(values[a],values[b])]
    if simple:bottom[key]=result
   else:
    m1,m2,m3=coefficients[new_index];new_index+=1
    if m1==full:
     result=[x^(y&m2)^((x&y)&m3)for x,y in zip(values[a],values[b])]
    else:
     result=[(x&m1)^(y&m2)^((x&y)&m3)for x,y in zip(values[a],values[b])]
   values.append(result)
  yield count,full,values[root]
def exception_extension_audit(envelope):
 started=time.monotonic()
 legal=[a for a in range(1024)if not(a&1 and a&2)
 and not(a&64 and a&128)and not(a&256 and a&512)]
 index={v:j for j,v in enumerate(legal)};assert len(legal)==432;flats=[];seen=set()
 for offset,*basis in EXTENSION_FLATS:
  assert len(basis)==4 and type(offset)is int and 0<=offset<1024;assert all(type(v)is int and 0<v<1024 for v in basis);points=[offset]
  for v in basis:points+=[x^v for x in points]
  assert len(set(points))==16 and all(v in index for v in points);key=tuple(sorted(points));assert key not in seen;seen.add(key)
  flats.append(tuple(index[v]for v in points))
 assert len(flats)==360;common=COMMON['extension'];total_assignments=total_hints=steps=cases=0
 for count,full,values in extension_values(legal):
  hints=[];total_assignments+=count;total_hints+=len(hints);steps+=check_common_cover(values,flats,common,hints,full);cases+=1
 assert cases==1450;pass
 return{'passed':True,'shapes':cases,'assignments':total_assignments,
 'verified_flats':len(flats),'hinted_flats':total_hints,'cover_steps':steps,
 'elapsed_seconds':time.monotonic()-started}
def one_pair_counterexamples():
 pass;legal=[a for a in range(64)if not(a&1 and a&2)]
 indexed=affine_flat_indices(legal,6,3);flats=[tuple(legal[j]for j in flat)for flat in indexed];assert len(flats)==1770
 group=[tuple(([1,0]if swap else[0,1])+list(p))
 for swap in(0,1)for p in itertools.permutations(range(2,6))]
 seen=set();trees=[]
 for tree in binary_trees(tuple(range(6))):
  key=min(tree_string(tree,p)for p in group)
  if key not in seen:seen.add(key);trees.append(tree)
 assert len(trees)==41;clear_tree_caches();count=4*3**5;full=(1<<count)-1;signs=[sum(((v>>i)&1)<<v for v in range(4))*(full//15)for i in range(2)];masks=[]
 for i in range(5):
  stride=4*3**i;m=full//((1<<(3*stride))-1)*((1<<stride)-1);masks.append((m|(m<<stride),(m<<stride)|(m<<(2*stride))))
 canonical=lambda a:int(bool(a&1 or a&4)^bool(a&2 or a&8)^bool(a&16 or a&32));orbit=set();truth_full=(1<<64)-1
 for shift in range(16):
  truth=sum(canonical(a^(shift<<2))<<a for a in range(64));orbit.update((truth,truth^truth_full))
 originals=(((0,2),((1,3),(4,5))),(((0,2),(1,3)),(4,5)));exceptions=[]
 for tree in trees:
  operations=[]
  def compile_tree(tree):
   if isinstance(tree,int):return tree
   a=compile_tree(tree[0]);b=compile_tree(tree[1]);operations.append((a,b));return 5+len(operations)
  root=compile_tree(tree);outputs=[]
  for point in range(64):
   values=[full if point>>j&1 else 0 for j in range(6)];values[0]^=signs[0];values[1]^=signs[1]
   for(a,b),(ma,mb)in zip(operations,masks):
    x,y=values[a],values[b];values.append(((x&y)&ma)^((x^y)&mb))
   outputs.append(values[root])
  left=full
  for flat in flats:
   zero=one=full
   for j in flat:one&=outputs[j];zero&=full^outputs[j]
   left&=full^(one|zero)
  while left:
   bit=(left&-left).bit_length()-1;left&=left-1;truth=sum(((v>>bit)&1)<<j for j,v in enumerate(outputs));assert truth in orbit and tree in originals
   value=bit//4;tables=[]
   for _ in range(5):tables.append((8,14,6)[value%3]);value//=3
   assert value==0
   # The three disjoint cherries are AND/OR; both joins are XOR.
   for(a,b),table in zip(operations,tables):
    assert table in((8,14)if a<6 and b<6 else(6,))
   exceptions.append({'tree':tree,'signs':bit%4,'tables':tuple(tables)})
 assert len(exceptions)==16;assert all(sum(row['tree']==tree for row in exceptions)==8 for tree in originals);return exceptions
def cleanup_legal(partner):
 base=[a for a in range(128)if not(a&1 and a&2)and not(a&4 and a&8)];remaining=[j for j in range(2,6)if j!=partner];permutation=(0,1,partner,6,*remaining)
 return[sum(((a>>j)&1)<<permutation[j]for j in range(7))for a in base]
def counter_tag(tree,operations):
 if isinstance(tree,int):return tree
 a=counter_tag(tree[0],operations);b=counter_tag(tree[1],operations);index=len(operations);operations.append((a,b));return('o',index,a,b)
def counter_masks(cases,full):
 repeat=full//65535;sign0=sum((row['signs']&1)<<i for i,row in enumerate(cases))*repeat;sign1=sum(((row['signs']>>1)&1)<<i for i,row in enumerate(cases))*repeat
 # Low three bits select one of eight exceptions, bit3 is the extra sign.
 sign0|=sign0<<8;sign1|=sign1<<8;flip=65280*repeat;old=[]
 for j in range(5):
  tables=[row['tables'][j]for row in cases]
  if set(tables)=={6}:old.append(None)
  else:
   assert set(tables)<={8,14};value=sum(int(t==14)<<i for i,t in enumerate(tables));old.append((value|(value<<8))*repeat)
 return sign0,sign1,flip,old
def zero_table_masks(options,stride,full):
 k=len(options);mask=full//((1<<(k*stride))-1)*((1<<stride)-1);m1=m2=m3=0
 for i,table in enumerate(options):
  block=mask<<(i*stride)
  if table&4:m1|=block
  if table&2:m2|=block
  if((table>>1)^(table>>2)^(table>>3))&1:m3|=block
 return m1,m2,m3
def reverse_A_values(exceptions):
 originals=list(dict.fromkeys(row['tree']for row in exceptions))
 for partner in range(2,6):
  legal=cleanup_legal(partner)
  for original in originals:
   cases=[row for row in exceptions if row['tree']==original];tree=counter_tag(original,[]);seen=set();shapes=[]
   for first in graft_new_leaf(tree,6):
    for second in graft_new_leaf(first,6):
     key=tagged_tree_string(second,tuple(range(7)))
     if key not in seen:seen.add(key);shapes.append(second)
   assert len(shapes)==77;prepared={}
   for tree in shapes:
    operations=[]
    def compile_tree(tree):
     if isinstance(tree,int):return tree,tree<6
     a,old_a=compile_tree(tree[-2]);b,old_b=compile_tree(tree[-1])
     if tree[0]=='o':slot=tree[1];options=None
     else:
      assert not old_b;slot=-1;options=(4,6,12,14)if old_a else(2,4,6,8,14)
     operations.append((a,b,slot,options));return 6+len(operations),old_a or old_b
    root,_=compile_tree(tree);assert len(operations)==7;choices=tuple(options for _,_,_,options in operations if options);assert len(choices)==2
    if choices not in prepared:
     count=16
     for options in choices:count*=len(options)
     full=(1<<count)-1;stride=16;new=[]
     for options in choices:
      new.append(zero_table_masks(options,stride,full));stride*=len(options)
     prepared[choices]=(count,full,new,counter_masks(cases,full))
    count,full,new,(s0,s1,flip,old)=prepared[choices];rows=[]
    for point in legal:
     row=[full if point>>j&1 else 0 for j in range(7)];row[0]^=s0;row[1]^=s1;row[partner]^=flip;rows.append(row)
    values=[list(column)for column in zip(*rows)];new_index=0
    for a,b,slot,options in operations:
     if slot>=0:
      mask=old[slot]
      if mask is None:result=[x^y for x,y in zip(values[a],values[b])]
      else:result=[(x&y)^((x^y)&mask)for x,y in zip(values[a],values[b])]
     else:
      m1,m2,m3=new[new_index];new_index+=1
      if m1==full:result=[x^(y&m2)^((x&y)&m3)for x,y in zip(values[a],values[b])]
      else:result=[(x&m1)^(y&m2)^((x&y)&m3)for x,y in zip(values[a],values[b])]
     values.append(result)
    yield partner,count,full,values[root]
def reverse_C_values(exceptions):
 originals=list(dict.fromkeys(row['tree']for row in exceptions))
 for partner in range(2,6):
  legal=cleanup_legal(partner);others=[j for j in range(6)if j!=partner]
  for original in originals:
   cases=[row for row in exceptions if row['tree']==original];tree=counter_tag(original,[])
   for coupled in(False,True):
    choices=others if coupled else[-1]+others;count=16*len(choices)*(1 if coupled else 4);full=(1<<count)-1;s0,s1,flip,old=counter_masks(cases,full)
    choice_masks=[];stride=16;k=len(choices);mask=full//((1<<(k*stride))-1)*((1<<stride)-1)
    for i in range(k):choice_masks.append(mask<<(i*stride))
    rows=[]
    for point in legal:
     row=[full if point>>j&1 else 0 for j in range(7)];selected=0
     for choice,mask in zip(choices,choice_masks):
      if choice>=0 and point>>choice&1:selected|=mask
      if coupled and point&64:row[choice]^=mask
     row[partner]=(full if point>>partner&1 else 0)^selected^flip
     if coupled and point&64:row[partner]^=full
     row[0]^=s0;row[1]^=s1;rows.append(row)
    base=[list(column)for column in zip(*rows)];new=zero_table_masks((4,6,12,14),16*k,full)if not coupled else None
    shapes=(tree,)if coupled else tuple(graft_new_leaf(tree,6));assert len(shapes)==(1 if coupled else 11)
    for shape in shapes:
     operations=[]
     def compile_tree(tree):
      if isinstance(tree,int):return tree
      a=compile_tree(tree[-2]);b=compile_tree(tree[-1]);operations.append((a,b,tree[1]if tree[0]=='o' else-1));return 6+len(operations)
     root=compile_tree(shape);values=base[:]
     for a,b,slot in operations:
      if slot>=0:
       mask=old[slot]
       if mask is None:result=[x^y for x,y in zip(values[a],values[b])]
       else:result=[(x&y)^((x^y)&mask)for x,y in zip(values[a],values[b])]
      else:
       m1,m2,m3=new;assert m1==full;result=[x^(y&m2)^((x&y)&m3)for x,y in zip(values[a],values[b])]
      values.append(result)
     yield partner,count,full,values[root]
def reverse_cleanup_audit(envelope,mode,exceptions):
 started=time.monotonic();assert mode in('A','C');legal=[a for a in range(128)if not(a&1 and a&2)and not(a&4 and a&8)]
 flats=affine_flat_indices(legal,7,3);assert len(flats)==5049;common=COMMON['reverse_'+mode];expected=616 if mode=='A' else 96
 generator=reverse_A_values if mode=='A' else reverse_C_values;count_cases=total_assignments=total_hints=steps=0
 for partner,count,full,values in generator(exceptions):
  hints=[];total_hints+=len(hints);total_assignments+=count;steps+=check_common_cover(values,flats,common,hints,full);count_cases+=1
 assert count_cases==expected;pass
 return{'passed':True,'cases':count_cases,'assignments':total_assignments,
 'hinted_flats':total_hints,'cover_steps':steps,
 'elapsed_seconds':time.monotonic()-started}
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
 common=COMMON['readonce7'];assert len(common)==128;assignments=16*3**6;full=(1<<assignments)-1;repeat=full//65535
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
def ramsey_direct_audit():
 started=time.monotonic();cases=[]
 for code,roots,canonical_A,units in ramsey84_normalized_roots():
  result=resolution_refutation(roots,95);assert result['sat']is False
  cases.append(dict(code=code,roots=len(roots),canonical_A=canonical_A,units=units,**result))
 assert len(cases)==4;return dict(passed=True,cases=cases,elapsed_seconds=time.monotonic()-started)
def midpoint_certificate():
 a=F(57893,125000);X=F(876987,1000000);wm=F(219247,250000);w0=F(283057,100000)
 pc=a/2-F(102,10**8)*w0-F(1,10**13);h=pc*X/F(2000001,10**6)-F(1,10**13);assert f.analytic_prefilter((h,a,wm,w0,X,pc))
 assert w0+X==F(3707557,1000000)<F(92689,25000)
 M=65536;step=(w0-wm)/M;scale=10**15
 V=lcm(wm.denominator,(step/2).denominator);U=int((wm+step/2)*V);inc=int(step*V);assert F(U,V)==wm+step/2 and F(inc,V)==step
 ap,aq=a.numerator,a.denominator;total=0;cross=0
 for j in range(M):
  z=U*aq-ap*V;num=scale*100000*U*z*z;den=V*(100000*aq*aq*U*U+300000*ap*aq*U*V+4*ap*V*z);assert den>0;integer=num//den;assert integer*den<=num<(integer+1)*den
  if j in(0,1,M//2,M-1):
   w=wm+(F(2*j+1,2))*step;beta=a/w;alpha=1/(1-beta)+F(1,100000);value=scale*w*(1-beta)/(1+4*alpha*beta)
   assert integer==value.numerator//value.denominator;cross+=1
  total+=integer;U+=inc
 assert total==41296754879821793427;ratio=8*h*step*F(total,scale)
 assert ratio==F(13653377949173357636608938738496614772013943,13653340160000000000000000000000000000000000)>1
 W=w0-wm;mu=ratio-1-8*h*W/10**12;assert mu>0;target=176*h*W*W/mu;J=1+-(-target.numerator//target.denominator)
 assert J>0 and ratio-8*h*W/10**12-88*h*W*W/J>1+mu/2>1
 return dict(terms=M,total=str(total),crosschecks=cross,left_mesh=str(J))

def run():
 result={'ramsey84_resolution':ramsey_direct_audit()};cases=[]
 for name,nv,roots in ramsey6_roots():
  outcome=resolution_refutation(roots,nv);assert outcome['sat']is False;cases.append(dict(case=name,roots=len(roots),**outcome))
 assert len(cases)==83;result['ramsey6_resolution']={'passed':True,'cases':cases};result['readonce7']=readonce_seven_audit(None)
 result['near_cross']=near_readonce_audit(None,'cross');result['near_mixed']=near_readonce_audit(None,'mixed');result['readonce8']=readonce_eight_audit(None)
 result['extension']=exception_extension_audit(None);exceptions=one_pair_counterexamples();assert len(exceptions)==16
 result['classified_exceptions']=len(exceptions);result['reverse_A']=reverse_cleanup_audit(None,'A',exceptions)
 result['reverse_C']=reverse_cleanup_audit(None,'C',exceptions)
 # New exact convex-midpoint potential certificate; calculus is supplied in proof.md.
 result['accuracy_midpoint']=midpoint_certificate()
 assert result['accuracy_midpoint']['left_mesh']=='24642872'
 result['all_passed']=True
 result['scope']='Direct checked resolution and exhaustive finite formula/cleanup obligations, plus exact midpoint witness; no universal theorem or accepted large-tree execution.'
 return result
if __name__=='__main__':
 request=json.loads(input())
 if request!={}:raise ValueError('empty object required')
 print(json.dumps(run(),sort_keys=True))
