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

# Plain exhaustive-cover witness; all plane IDs and truth coverage are checked.
TERMINAL_COVER={"version":1,"common":[12138,7520,133,421,14337,12139,13666,12934,13987,13718,12289,10822,14120,14376,2051,12313,5145,13153,11890,8628,12141,7508,9725,1454,13597,9371,12334,13526,13149,6724,13166,5133,5121,6814,13614,18,13492,6934,14107,5624,3799,13269,12487,12639,13623,382,12423,11030,143,6554,13782,10945,11174,13495,3929,11274,13339,11991,12920,4178,9964,13151,7232,13096,10075,6846,13885,8095,14365,12531,13627,601,114,9978,14132,183,251,9963,6606,1640,58,290,4650,3904,1398,12326,3785,12149,564,12381,6336,9339,14316,14179,9143,14343,10889,7428,13574,1971,5638,6498,13752,2495,5787,13670,10888,101,12391,5110,12631,11343,2268,11220,5212,1080,7483,8,13609,13570,10341,585,13665,10971,12333,13703,5134,6933,12603,7189,7459,4729,12392,12619,13085,13875,9137,533,11186,9498,13047,6560,12962,13207,4402,10097,452,12454,12655,13653,5880,5104,10819,71,13363,339,12354,1220,4960,12924,6860,4875,14069,6381,11985,758,4964,5101,6881,12190,8076,3946,14174,14004,2242,12482,9266,5125,10909,1194,12987,12481,13347,12578,314,12771,4245,13667,12948,11222,3820,12394,5393,6618,14104,930,142,7434,2132,4892,2826,3846,12504,4280,13799,1184,11148,3720,6097,10995,6138,6507,3761,11901,13583,3796,9136,14301,1098,12198,13215,629,9590,7512,1742,3873,14279,1499,1406,6324,13990,12773,13748,13995,6385,12923,185,11059,6335,4180,9709,2254,62,1157,13765,8216,4954,3749,12638,13143,13360,12709,14384,2036,6732,4559],"covers":[[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[342],[12396],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[1016,12042],[14002,938],[],[],[11673],[],[],[9204],[],[6937,6973],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[5284],[],[],[7751],[],[267],[282],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[11002],[],[],[6160],[],[],[],[1736],[12396],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[6603],[14002,4299],[],[],[],[3933],[608,911],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[1227,12247],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[4685],[],[],[],[],[],[13506],[4014],[],[6973],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[3813],[4014],[],[5092],[],[14003],[14003],[],[],[],[],[],[8155],[7465],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[1016,12409],[13826,14003],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[342],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[4083],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[11483],[],[6937,331],[],[10319],[],[6973],[],[13722],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[9204],[342,7075],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[4884],[331,11192],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[9210,12242],[],[],[12251],[],[267],[267],[],[],[],[],[],[],[],[],[],[13726,11442],[],[],[],[],[],[],[],[8155],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[356],[],[6065,13094],[13107,12222],[],[13094],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[1016,12042],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[1016,4508],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[282],[8155],[],[],[],[2185],[],[],[],[2996],[],[],[],[],[11356],[],[],[],[],[],[],[],[],[],[],[],[928],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[11356],[],[],[3831],[],[],[],[14001],[11002],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[8155],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[4171],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[7075],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[1016,4401],[7699],[],[],[],[3855],[1484,273],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[3813],[4014],[],[608],[],[14003],[14199],[],[],[],[],[],[],[],[],[],[],[],[],[1016,4401],[11469],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[3855],[928],[],[1002],[2015],[],[6603],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[1004],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[5626],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[1016,4401],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[4401,1016],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[968],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[5626],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[1016,4401],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[]]}

def terminal_cover_audit(certificate):
    started=time.monotonic()
    n,d=8,3
    legal=[a for a in range(256) if not(a&1 and a&4)
           and not(a&2 and a&8) and not(a&16 and a&32)]
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
    assert len(legal)==108 and len(planes)==14391
    def trees(labels):
        if len(labels)==1:return [labels[0]]
        out=[];first=labels[0];rest=labels[1:]
        for k in range(len(rest)):
            for more in itertools.combinations(rest,k):
                left=(first,)+more;right=tuple(x for x in labels if x not in left)
                for a in trees(left):
                    for b in trees(right):out.append((a,b))
        return out
    def serial(tree,perm):
        if isinstance(tree,int):return str(perm[tree])
        a,b=serial(tree[0],perm),serial(tree[1],perm)
        return '('+','.join(sorted((a,b)))+')'
    permutations=[]
    for a,b,e in itertools.product((0,1),repeat=3):
        p=list(range(7))
        if a:p[1],p[2]=p[2],p[1]
        if b:p[3],p[4]=p[4],p[3]
        if e:p[5],p[6]=p[6],p[5]
        permutations.append(p)
    all_trees=trees(tuple(range(7)));seen=set();representatives=[]
    for tree in all_trees:
        key=min(serial(tree,p) for p in permutations)
        if key not in seen:seen.add(key);representatives.append(tree)
    assert len(all_trees)==10395 and len(representatives)==1710
    assert set(certificate)=={'version','common','covers'} and certificate['version']==1
    common=certificate['common'];covers=certificate['covers']
    assert isinstance(common,list) and len(common)==256 and len(set(common))==256
    assert len(covers)==len(representatives)
    assert all(type(i) is int and 0<=i<len(planes) for i in common)
    for row in covers:
        assert isinstance(row,list) and all(type(i) is int and 0<=i<len(planes) for i in row)
        assert len(row)==len(set(row))
    assignments=64*3**6;full=(1<<assignments)-1
    repeat=full//((1<<64)-1)
    signs=[sum(((j>>i)&1)<<j for j in range(64))*repeat for i in range(6)]
    op_masks=[]
    for i in range(6):
        stride=64*3**i
        a=(full//((1<<(3*stride))-1))*((1<<stride)-1)
        op_masks.append((a,a<<stride,a<<(2*stride)))
    total_planes=0
    for tree,ids in zip(representatives,covers):
        operations=[]
        def compile_tree(t):
            if isinstance(t,int):return t
            a,b=compile_tree(t[0]),compile_tree(t[1])
            operations.append((a,b));return 6+len(operations)
        root=compile_tree(tree);assert root==12 and len(operations)==6
        points={j for pi in common+ids for j in planes[pi]};values={}
        for j in points:
            a=legal[j]
            leaves=[signs[i] if not(a>>i&1) else full^signs[i] for i in range(6)]
            leaves.extend(full if a>>i&1 else 0 for i in (6,7))
            v=[leaves[0]&leaves[1]]+leaves[2:]
            for (left,right),(aa,oo,xx) in zip(operations,op_masks):
                x,y=v[left],v[right]
                v.append(((x&y)&aa)|((x|y)&oo)|((x^y)&xx))
            values[j]=v[root]
        covered=0
        for pi in common+ids:
            if covered==full:break
            zero=one=full
            for j in planes[pi]:
                one &= values[j];zero &= full^values[j]
            covered |= zero|one
        assert covered==full
        total_planes +=len(ids)
    return {'passed':True,'labeled_component_trees':10395,'shape_orbits':1710,
            'domain_points':108,'affine_three_flats':14391,
            'normalized_assignments_per_shape':assignments,
            'normalized_assignments_verified':assignments*len(representatives),
            'cover_plane_witnesses':total_planes,'shared_planes':len(common),'certificate_data_bytes':len(json.dumps(certificate,separators=(',',':')).encode()),
            'elapsed_seconds':time.monotonic()-started}



def run():
 result={'terminal_read_once_cover':terminal_cover_audit(TERMINAL_COVER)}
 h,a,wm,w0,X,pc,D=F(7180240479,100000000000),F(86061,156250),F(1042926401,1000000000),F(1683069,500000),F(651829,625000),F(68847083257,500000000000),F(881813,200000)
 assert f.analytic_prefilter((h,a,wm,w0,X,pc)) and D>w0+X
 M=4194304;step=(w0-wm)/M
 V=lcm(wm.denominator,step.denominator);U=wm.numerator*(V//wm.denominator);inc=step.numerator*(V//step.denominator)
 assert F(U,V)==wm and F(inc,V)==step
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
 assert total==3142984930683370297
 ratio=8*h*step*F(total,10**12)
 assert ratio==F(52428816687623435373084947320302198537,52428800000000000000000000000000000000)>1
 result['current_exact_cost_grid_explicit_box_certificate_terms']=M
 result['current_exact_cost_grid_explicit_box_certificate_sum']=str(total)
 result['integer_box_independent_rational_crosschecks']=cross
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
 result['scope']='Exact exhaustive read-once terminal lemma, retained left potential sum and new convex-midpoint witness/error bound only; no complete accepted tree or universal circuit/extractor execution.'
 return result

if __name__=='__main__':
 request=json.loads(input())
 if request!={}:raise ValueError('empty object required')
 print(json.dumps(run(),sort_keys=True))
