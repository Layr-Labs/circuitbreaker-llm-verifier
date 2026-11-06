#!/usr/bin/env python3
"""Exact finite witnesses for the printed constants; no universal proof or sampling.
Run: echo '{}' | python3 verify.py. JSON stdin/stdout; integer/Fraction arithmetic.
Relations are explicit below. Source locations refer to the accompanying files.
"""
import json, sys
from fractions import Fraction as F
from math import ceil, floor

count = failures = 0
LINES = []

def check(name, src, lhs, rhs, op='>='):
    global count, failures
    ok = {'>=': lhs >= rhs, '>': lhs > rhs, '<=': lhs <= rhs,
          '<': lhs < rhs, '==': lhs == rhs}[op]
    count += 1
    failures += not ok
    # Full values only on failure; the verifier caps stdout at 64 KiB.
    LINES.append(f"PASS {name} {src} {op}" if ok else f"FAIL {name} {src} lhs={lhs} {op} rhs={rhs}")

def lg(x):
    x = F(x)
    j = max(0, x.numerator.bit_length() - x.denominator.bit_length())
    return j + (2**j < x)

P = 'proof.md:'
W = 'affine-and-standard-wrappers.md:'
V = 'linear-and-bilinear-primitives.md:'
B = 'linear-and-bilinear-primitives.md:'
S = 'sampler-and-disperser.md:'
EPS = F(1, 2**20)
# affine-and-standard-wrappers.md:484 and proof.md:697: selected vectors.
specs = [
 ('I',8,'1124607/500000','7/50',True),
 ('F',8,'2063933/1000000','1/5',True),
 ('G',2,'3444761/1000000','7/20',True),
 ('J',2,'6019061/1000000','19/200',True),
 ('Ext0',22,'11300207/1000000','1',False),
 ('Ext2',5,'19574499/1000000','6/5',False),
 ('L0',13,'836983/500000','7/20',True),
 ('Lr',22,'1149889/500000','171/200',True),
 ('E',5,'2153593/1000000','5/32',True),
 ('Lm',1,'7290999/1000000','1',False)]

def envelope(name,h,R,b,aghp):
    R,b = F(R),F(b)
    rates = [1+b*(R-1)]*(h-1)+[R]
    # linear-and-bilinear-primitives.md:636: descending A,B,C recurrence.
    a, bb, cc = F(2), F(1,2**h), F(1)
    product = F(1)
    for i, ri in reversed(list(enumerate(rates,1))):
        check(f'PXT.{name}.rate{i}',V+'626',ri,1,'>')
        lam = 2*ri/(ri-1)
        a,bb,cc = lam*(a+1+(ri+1)/2),lam*bb,lam*(cc+F(3,2))+1
        # Independently expand this suffix into prefix products and forcing sums.
        prod, sa, sc, sg = F(1), F(0), F(0), F(0)
        for r in rates[i-1:]:
            prev = prod
            prod *= 2*r/(r-1)
            sa += prod*(r+3)/2
            sc += F(3,2)*prod+prev
            sg += prod*((r+3)/2+F(3,2)*EPS)+prev*EPS
        check(f'PXT.{name}.A{i-1}',V+'637',a,2*prod+sa,'==')
        check(f'PXT.{name}.B{i-1}',V+'637',bb,prod/2**h,'==')
        check(f'PXT.{name}.C{i-1}',V+'638',cc,prod+sc,'==')
        product *= ri/(ri-1)
    check(f'PXT.{name}.Lam',V+'641',bb,product,'==')
    g = a+EPS*cc
    check(f'PXT.{name}.g',V+'640',g,(2+EPS)*prod+sg,'==')
    kap,tau = (R+1)/(R-1),(R+1)/(2*R)
    if aghp:
        # linear-and-bilinear-primitives.md:670: PXG terminal seed envelope.
        c = F(2,2**h)
        pi = 6+sum(r+1 for r in rates)+F(1,8192)+(h+2)*EPS
        check(f'PXG.{name}.rate-min',V+'664',R,F(5,3))
        check(f'PXG.{name}.rate-max',V+'664',R,64,'<=')
    else:
        # linear-and-bilinear-primitives.md:642: PXT terminal seed envelope.
        c = kap/2**h
        pi = (2+sum(r+1 for r in rates)+kap*(3+(R+1)/2)+(R+1)/2
              +EPS*(h+F(5,2)*kap+tau+F(3,2)))
    return c,pi,bb,g

v = {name:envelope(name,h,R,b,ag) for name,h,R,b,ag in specs}
# affine-and-standard-wrappers.md:477: rational slots and all charged rows.
a = dict(zip(('a','s','s0','h','d','k'),map(F,[
 '4058372281055/137438953472','30857602586255/1099511627776',
 '16189511144005/549755813888','700165329211955/1099511627776',
 '1140949273153177/274877906944','605411727703975008315/274877906944'])))
ci,pi,li,gi = v['I']; cf,pf,lf,gf = v['F']
cg,pg,lg_,gg = v['G']; cj,pj,lj,gj = v['J']
rows = [
 ('seed-F','a',{'h':cf},pf,0),
 ('seed-J','a',{'s':cj},pj,146797863),
 ('seed-G','s',{'a':cg},pg,0),
 ('seed-I','s0',{'h':ci},pi,0),
 ('entropy-J','h',{'s':15+lj},gj+1,0),
 ('entropy-G','d',{'s0':15,'a':120+lg_},gg+1,0),
 ('entropy-F','k',{'h':134+lf},gf+1,1609194937627445464019),
 ('entropy-I','k',{'h':li},gi,0),
 ('input-d','k',{'d':F(1)},0,2421642347018806371976)]
for row,(name,x,weights,u,cert) in enumerate(rows):
    rhs = sum(c*(a[y]+EPS) for y,c in weights.items())+u
    check('F#.'+name,W+str(505+row),a[x],rhs)
    check('F#.certificate-'+name,W+str(505+row),floor(2**40*(a[x]-rhs)),cert,'==')
for x,y in [('h','s'),('d','s0'),('k','h'),('k','d')]:
    check('F#.order-'+x+'-'+y,W+'517',a[x],a[y])
check('F#.capacity',W+'516',a['h'],128)
check('F#.error',W+'494',66,2**7,'<')
check('F#.unit-height',W+'523',7+9+6,32,'<=')
# proof.md:348: B9 exact rational substitution, with its distinct epsilon.
nb = F(71,137438953472)
eta = F(256000000,2109543000001)
theta = eta/8
c0,p0,l0,g0 = v['Ext0']; c2,p2,l2,g2 = v['Ext2']
df,kf = a['d']+nb,a['k']+nb
M0 = c2*df+p2*(1+nb)+nb
K2 = l2*df+g2*(1+nb)+nb
check('B9.denominator',P+'363',theta,7*c0,'>')
Rnew = (M0+p0*(1+nb)+7+7*nb)/(theta-7*c0)+nb
d0 = 7*c0*Rnew+p0*(1+nb)+nb
d1 = 2*(M0+d0+3)+1
K00 = 7*l0*Rnew+g0*(1+nb)
Kother = kf+14*Rnew+84+nb
check('B9a.K00',P+'369',K00,2874390323,'<')
check('B9a.Kother',P+'370',Kother,2874390323,'<')
check('B9a.dominant',P+'372',K00,Kother,'>')
check('B9b.seed',P+'378',2*d0+2*d1+df+K2+1,13792,'<')
check('B9.log-constant',P+'374',lg(ceil(Rnew*M0)+1)+1,240,'<')
check('B9.log-endpoint',P+'349',8*41+240,nb*2**40,'<=')
check('B9.log-induction',P+'349',8*42+240,2*(8*41+240),'<=')
check('B9.quarter-endpoint',P+'375',2*41+240,F(2**40,4),'<')
check('B9.quarter-induction',P+'375',2*42+240,2*(2*41+240),'<')
check('B9.unit-additive',P+'352',14,nb*2**40,'<=')
check('B9.unit-height',P+'353',max(22,5),38,'<=')
check('B9.capped-output',P+'355',82,128,'<=')
check('B9.numerator',P+'366',F(41,8)+1+F(50,2**40),7,'<=')
check('B9.prime-threshold',W+'596',8*2**20/eta,3200000,'>')
check('B9.pi0',W+'596',p0,1)
check('CZ.expansion',P+'322',103*eta,F(256,20481),'<')
check('CZ.eta',P+'323',eta,F(1,8000),'<')
check('CZ.rate',W+'591',eta*(1+eta),1,'<')
check('CZ.error',P+'330',14,2**4,'<')
# proof.md:674: graph, advice and error constants.
B0,P3,Q = 392,215,144
p = 3*P3*B0
q = (P3-1)*B0//2
C = 4*(q+1)
t = 4*C-1
for name,line,x,y in [('p',674,p,252840),('q',674,q,41944),('C',674,C,167780),
                  ('t',397,t,671119),('lgC',678,lg(C),18),('lg26t',678,lg(26*t),25),
                  ('26t',678,26*t,17449094),('j',678,lg(t),20),
                  ('e-intercept',680,33+lg(C)+lg(26*t),76)]:
    check('outer.'+name,P+str(line),x,y,'==')
cO = F(6*p+1,Q)+2+EPS
check('outer.cO',P+'688',cO,F(99439673353,9437184),'==')
delta = ceil(12156847042*cO)
kappa = ceil(59628486398890*cO)
check('outer.delta',P+'702',delta,128096781715697,'==')
check('outer.kappa',P+'704',kappa,628305775328681190,'==')
ce,pe,le,ge = v['E']; cr,pr,lr,gr = v['Lr']
cm,pm,lm,gm = v['Lm']; cz,pz,lz,gz = v['L0']
TR = (t+1)*(39+13792+20*(46+164))
TL = (t+1)*(13792+40*23)
check('outer.TR',P+'709',TR,12100964720,'==')
check('outer.TL',P+'710',TL,9873517440,'==')
for name,line,x,y in [
 ('seed-E',713,23,ce*82+pe),('seed-Lr',713,82,cr*15435805+pr),
 ('seed-Lm',714,82,cm*23+pm),('seed-L0',714,39,cz*13792+pz),
 ('r',715,15435805,(t+1)*23+lm*23+gm+1),
 ('d-NB',716,12156847042,(t+1)*(39+13792)+2874390323+1),
 ('d-E',717,12156847042,TR+le*82+ge+t*82+1),
 ('k-Lr',718,59628486398890,TL+lr*15435805+gr+t*15435805+1),
 ('k-L0',719,59628486398890,lz*13792+gz)]:
    check('A.'+name,P+str(line),x,y)
for name,line,x,y in [('s0-d',726,39,12156847042),('dy-r',726,23,15435805),
                  ('m-r',726,1,15435805),('r-k',726,15435805,59628486398890),
                  ('d-k',726,12156847042,59628486398890)]:
    check('A.order-'+name,P+str(line),x,y,'<=')
outputs = dict(I=a['h']+EPS,F=a['h']+EPS,G=a['a']+EPS,J=a['s']+EPS,
               L0=13792,Lr=15435805,E=82)
for name,h,R,b,ag in specs:
    if ag: check('PXG.'+name+'.output',V+'664',outputs[name],256*2**h,'<=')
check('outer.error',P+'735',13*2**20-8,26*t,'<')
check('outer.unit-height',P+'731',23+6+6,40,'<=')
# linear-and-bilinear-primitives.md:522: base 7f envelopes.
alpha = F(384290345147091,93824902758400)
Ro = F(50981815897738182537983470459,6296488643820466995200000000)
Rl = F(6386477,1000000)
R0 = 1099511627776
lo,ll = 4*Ro/(Ro-3),4*Rl/(Rl-3)
spo,spl = (Ro+1)*(1+EPS)+EPS,(Rl+1)*(1+EPS)+EPS
kl = ll*(alpha+EPS+3*(1+EPS)+F(3,4)*spl+EPS)+EPS
f = lo/(lo-1)
kfac = lo**25*(kl+f*((1+EPS)+F(3,4)*spo+EPS)+EPS)
ceilk = ceil(kfac)
check('base.Ro',B+'522',Ro,3,'>'); check('base.Rl',B+'522',Rl,3,'>')
check('base.last-block',B+'553',15*Rl,kl)
check('base.leaf',B+'555',alpha+EPS,16,'<=')
check('base.last-dimension',B+'555',F(15,4)*spl+EPS,64,'<=')
check('base.kfac-ceiling',B+'564',ceilk,12484765614829834828034,'==')
check('base.kfac-lower',B+'564',ceilk-1,kfac,'<')
check('base.kfac-upper',B+'564',kfac,ceilk,'<=')
check('base.kfac-bits',B+'564',ceilk,2**75,'<')
check('base.reverse-ceilings',B+'568',F(1,(lo-1)*2**40),EPS,'<=')
i0 = F(R0+1,R0)*(145+EPS)+EPS
seed = F(R0+1,R0)*(145+EPS)+25*spo+spl+2*alpha+4+10*EPS
check('base.seed',S+'320',seed,B0,'==')
check('base.padding',S+'324',(R0+1)*ceilk+i0,2**128,'<')
check('base.output',S+'329',2**52*alpha*(1-EPS),144*delta,'==')
check('base.initial-rate',S+'325',R0,2**40+64,'<=')
check('base.failure-bits',S+'332',2**7,100)
check('base.error',B+'545',3*4**26,2**64,'<')
# proof.md:754: two branches, all-length guards, K.
T = ceil(F(12*p+2,144)+EPS)
raw = 4*C*delta+kappa+T
base = ceil(kfac/144)+7
H = max(raw,base)
for name,line,x,y in [('T',755,T,21071),('H.raw',757,raw,86596617920367272901),
                  ('H.base',757,base,86699761214096075202),
                  ('H',757,H,86699761214096075202),('K',899,4*H,346799044856384300808)]:
    check(name,P+str(line),x,y,'==')
Smin = 2**40
Lmin = 144*Smin
for name,src,x,y,op in [
 ('B9.half-input',P+'346',H,2*delta,'>='),('guard.Om',P+'343',(6*p+1)*Smin+76+2*Lmin+40,2**40,'>='),
 ('guard.primitive-unit',V+'635',Smin,2**20,'>='),('guard.outer-intercept',P+'689',EPS*Lmin,116,'>='),
 ('guard.T-intercept',P+'759',EPS*Lmin,110,'>='),('guard.tag-slope',P+'763',p,129*144,'>'),
 ('guard.tag-intercept',P+'763',Lmin,p+128,'>='),('guard.tag-length',P+'764',p*Smin,128*(Lmin+1),'>='),
 ('guard.tag-blocks',P+'766',128,2,'>='),('guard.root-length',P+'765',Lmin,2048,'>='),
 ('guard.N-exponent',P+'779',Lmin,256,'>='),('guard.H',P+'780',H,2**128,'<'),
 ('guard.all-n-anchor',P+'780',2**256,2**128*256,'>='),('guard.all-n-doubling',P+'781',2*256,257,'>='),
 ('guard.base-endpoint',S+'310',600+2*40,2**20,'<='),
 ('guard.base-induction',S+'310',600+2*41,2*(600+2*40),'<='),
 ('guard.rho1',S+'312',322+2*41,600+2*40,'<='),
 ('guard.rho0',S+'312',287+41,600+2*40,'<='),
 ('guard.floor-output',S+'328',1,EPS*(Smin+1),'<='),
 ('guard.base-log',B+'557',13+2*40,2**20,'<='),
 ('guard.base-log-induction',B+'557',13+2*41,2*(13+2*40),'<='),
 ('guard.PXG-log',V+'667',22+2*20,2**6,'<='),
 ('guard.PXG-log-induction',V+'667',22+2*21,2*(22+2*20),'<='),
 ('guard.PXG-dimension',V+'666',1199*2**20+9,2048*2**20,'<='),
 ('guard.PXG-product',V+'666',2048*256,2**19,'==')]:
    check(name,src,x,y,op)
# sampler-and-disperser.md:145: walk and support thresholds.
xi = F(1,427)
mu = (P3-1)*xi-F(1,2)
for name,line,x,y,op in [
 ('walk.even',148,2*q,(P3-1)*B0,'=='),('walk.mu',160,mu,F(1,854),'=='),
 ('walk.positive',160,mu,0,'>'),('walk.log',167,B0+2,512,'<='),
 ('walk.endpoint',168,(P3-1)*42,mu*2**32,'<'),
 ('walk.induction',168,(P3-1)*43,2*(P3-1)*42,'<='),
 ('walk.floor',201,(mu*B0/2-xi)*1024,1,'>='),
 ('walk.support',222,P3*B0*4096,300000000,'>'),
 ('walk.guard',228,Smin,2**32,'>='),('walk.local-guard',215,2**32,4096,'>='),
 ('walk.Lambda-integer',221,p%2,0,'=='),('walk.Lambda-exponent',221,p*4096//2,15,'>=')]:
    check(name,S+str(line),x,y,op)
# sampler-and-disperser.md:806, RG0 and S1' numerical margins.
er = F(10,711)
rr = (1-2*xi)/(1+2*xi)
for name,line,x,y,op in [
 ('RG.small',806,800000*(1-71*er),252,'>'),
 ('RG.large',806,800000*(rr-34*er),122,'>'),
 ('RG.nontrivial',847,800000*(1-3*er),7,'>'),
 ('RG.small-constant',851,64*33*2**240,2**252,'<'),
 ('RG.large-constant',853,3*2**120,2**122,'<'),
 ('S1.log',886,F(300000000,2)-6,800000,'>='),
 ('S1.slope',900,er/2-3*xi,F(2,303597),'=='),
 ('S1.margin',900,F(600000000,303597),18,'>'),
 ('S1.intercept',900,17+6*er,18,'<'),
 ('S1.zero-mass',898,800000*(1-er),2,'>')]:
    check(name,S+str(line),x,y,op)
# proof.md:622: finite coefficients in the tag/quartic argument.
# D3: conservative endpoint A=2^128, S=0; 2^-S<=1 for actual S>=2^40.
# After multiplying by A^2 the remaining term decreases with A and S.
bias_A, bias_S = 2**128, 0
bias_theta = F(1,1024*bias_A**2)
bias_epsilon = F(1,2**(bias_S+33)*C*bias_A**6)
for name,line,x,y,op in [
 ('D2.union',564,64*2**16*2*1024,2**33,'=='),
 ('D1.collision',474,F(2,2048)+F(1,4),F(257,1024),'=='),
 ('D3.good-fraction',573,F(1,2**14),F(1,256),'<='),
 ('D3.bias',581,2*(bias_theta+bias_epsilon),F(3,1024*bias_A**2),'<='),
 ('D4.collision-ratio',622,F(257,1024)*F(256,255),F(257,1020),'=='),
 ('D4.r',623,F(257,1020)+F(3,1024*2**256),F(13,50),'<'),
 ('D4.monotone',629,6-8*F(13,50),0,'>'),
 ('D4.quartic',630,6*F(13,50)-4*F(13,50)**2-2,F(-444,625),'=='),
 ('D4.g',630,F(4*256,255*2**128),F(1,1000),'<'),
 ('D4.margin',626,F(-444,625)+F(1,1000)+F(12,1024),F(-1,2),'<='),
 ('D4.delta-squared',637,F(1,2**27),F(1,10000**2),'<'),
 ('D5.perturbation',651,F(432,10000),F(1,8),'<'),
 ('D4.fourth-moment',617,3+F(3,1024),4,'<='),
 ('D3.success',665,1-F(2,100),F(98,100),'==')]:
    check(name,P+str(line),x,y,op)
RESULT = {"all_passed": not failures, "checks": count, "failures": failures,
          "H_raw": str(raw), "H_base": str(base), "H": str(H), "K": str(4*H), "lines": LINES,
          "scope": "Exact finite check of the displayed constants only; universal statements are proved in the Markdown files."}

if __name__ == "__main__":
    if json.load(sys.stdin) != {}:
        raise ValueError("verify.py takes an empty input object")
    print(json.dumps(RESULT, sort_keys=True))
    raise SystemExit(0 if RESULT["all_passed"] else 1)
