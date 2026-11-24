"""Directed fair-cut certificate with a reconstructed intermediate table."""
from fractions import Fraction as F
from math import isqrt
S=1<<50;H=4096
def up(n,d): return -(-n//d)
def pi():
    def atan(x):
        v=sum((-1)**j*x**(2*j+1)/(2*j+1) for j in range(20))
        return v,v+x**41/41
    a,b=atan(F(1,5));c,d=atan(F(1,239));return 16*a-4*d,16*b-4*c
PI=pi()
class Gaussian:
    __slots__=('q','dg','log','p','cum','dens','err','tail')
    def __init__(self,q,logs=26):
        if q not in (1,2) or not 12<=logs<=29:raise ValueError('CDF parameters')
        self.log=logs;self.p=1<<logs;self.q=q;self.dg=S//self.p;self.tail=1
        def rt(x):return isqrt(S*S*x.numerator//x.denominator)
        cc,cd=rt(1/(q*PI[1])),1+rt(1/(q*PI[0]));D=4*H*H*256*q
        c=[(S//2,S//2)];density=[]
        for k in range(8*H):
            z=(2*k+1)**2*S;a,b=z//D,up(z,D)
            # Enclose consecutive odd/even exp(-x²/(256q)) partial sums.
            m=n=S;r=s=S
            for j in range(1,12):
                m=m*a//(S*j);n=up(n*b,S*j)
                if j&1:r-=n;s-=m
                else:r+=m;s+=n
            s+=up(n*b,12*S)
            for _ in range(8):r=r*r//S;s=up(s*s,S)
            f,g=r*cc//S,up(s*cd,S)
            h,i=c[-1];c.append((h+f//H,i+up(g,H)));density.append((f,g))
        self.cum,self.dens=c,density
        self.err=up(8*S,12*H*H)+up(S,H*H)
    def __call__(self,n,d):
        if d<=0:raise ValueError('CDF positive denominator')
        if n<0:
            a,b=self(-n,d);return self.p-b,self.p-a
        if n>=8*d:return self.p-1,self.p
        k,l=divmod(n*H,d);a,b=self.cum[k];f,g=self.dens[k]
        a+=f*l//(d*H);b+=up(g*l,d*H)
        return max(0,(a-self.err)//self.dg),min(self.p,up(b+self.err,self.dg))

"""All-positive finite Gaussian propagation; exact nonnegative integer products.
Consecutive identical transitions are coalesced; ±1 states have analytic merges.
"""
from operator import itemgetter,mul
from math import sumprod as dot
def ce(n,d):return -((-n)//d)
def get(ids):
    # itemgetter(one) returns a scalar, so retain a zero-weight dummy entry
    return itemgetter(*ids) if len(ids)>1 else itemgetter(ids[0],ids[0])
def reducer(xs,ys,W,T,N):
    reg={};neg={};pos={};other={};posQ=0
    for x,y,w in zip(xs,ys,W):
        if not w:continue
        if y==N:
            if x!=T:neg[x]=neg.get(x,0)+w
        elif x==T and y==T:posQ+=w
        elif y==T:pos[x]=pos.get(x,0)+w
        elif x==T:other[y]=other.get(y,0)+w
        else:k=x,y;reg[k]=reg.get(k,0)+w
    def single(d):return (get(list(d)),tuple(d.values())+((0,) if len(d)==1 else ())) if d else None
    bil=(get([a for a,b in reg]),get([b for a,b in reg]),tuple(reg.values())+((0,) if len(reg)==1 else ())) if reg else None
    return bil,single(neg),single(pos),single(other),posQ
def compute_plan(qshape,m,p,ns,ps,rest,Q):
    bits = 18
    scale = 1 << bits
    bucket, remainder = divmod(Q, scale)
    if not bucket or remainder:
        raise ValueError("incompatible probability scales")
    m = [x // bucket for x in m]
    p = [(y + bucket - 1) // bucket for y in p]
    ns = [x // bucket for x in ns]
    ps = [(y + bucket - 1) // bucket for y in ps]
    qp = [scale - y for y in p]
    qm = [scale - x for x in m]
    on = [scale - x for x in ns]
    z = [scale - x - y for x, y in zip(ns, ps)]
    K = 2 * bits
    M=[];P=[];ma=M.append;pa=P.append;ceil=(1<<K)-1
    for B,N,L,R,sm in qshape:
        u=0;v=(rest+sm)*scale*scale
        if B:
            f,g,w=B;zz=g(z)
            u=dot(map(mul,w,f(m)),zz)+dot(map(mul,w,f(qp)),g(ns))
            v+=dot(map(mul,w,f(p)),zz)+dot(map(mul,w,f(qm)),g(ps))
        if N:f,w=N;u+=scale*dot(w,f(qp))
        if L:f,w=L;v+=scale*dot(w,f(qm))
        if R:f,w=R;v+=scale*dot(w,f(on))
        ma(u>>K);pa((v+ceil)>>K)
    if any(x+y>Q or x<0 or y<0 for x,y in zip(M,P)):raise ArithmeticError('tail simplex')
    return M,P
def run_layers(Q,half_mass,T=-11200,dx=800,dw=1000,dep=10,dm=35,pm=60,*,debug=False,initial=None):
    if debug:import time;t=time.process_time()
    G=40000; nd=2*dm+1;nj=2*pm+1;size=nd*nj;nw=ce(3*G,dw);N=2*nw
    wg=tuple(max(0,half_mass((g+1)*dw,G)[0]-half_mass(g*dw,G)[1]) for g in range(-nw,nw));aw=tuple(2*w for w in wg[nw:]);rest=Q-sum(wg)
    if rest<0 or wg!=wg[::-1] or min(wg)<0:raise ArithmeticError('mass')
    ids=list(range(size));top=size;neg=top+1;K=2*(Q.bit_length()-1)
    if (1<<K)!=Q*Q or dep<1:raise ValueError('scale/depth')
    def ix(i,j):return top if i>=nd or j>=nj else ids[(j if j>=0 else 0)*nd+(i if i>=0 else 0)]
    ds=[[max(0,dm+ce(2337*(i-dm)*dx+1568*(g+1)*dw,3249*dx)) for g in range(-nw,nw)] for i in range(nd)]
    dh=[[max(0,dm+ce(1596*(i-dm)*dx-1148*g*dw,3249*dx)) for g in range(-nw,nw)] for i in range(nd)]
    ph=[[ce(57*r*dx-2045*g*dw,3249*dx) for g in range(-nw,nw)] for r in range(57)]
    ps=[[ce(57*r*dx+1568*(g+1)*dw,3249*dx) for g in range(-nw,nw)] for r in range(57)]
    au=[ce(28*(v+1)*dw,57*dx) for v in range(nw)];al=[ce(-28*v*dw,57*dx) for v in range(nw)]
    Am=[];Cm=[];ac=Am.append;cc=Cm.append;m=[];p=[]
    for j in range(nj):
        cut_num = 57*((j-pm)*dx-T)
        bin_index, remainder = divmod(cut_num,13*dw)
        cut=max(0,min(N,nw+ce(cut_num,13*dw)))
        weights = wg
        extra = 0
        split = bool(remainder and -nw <= bin_index < nw)
        if split:
            at = bin_index + nw
            cut_lo, cut_hi = half_mass(cut_num,13*G)
            lower = half_mass(bin_index*dw,G)[1]
            upper = half_mass((bin_index+1)*dw,G)[0]
            first = min(wg[at],max(0,cut_lo-lower))
            second = min(wg[at]-first,max(0,upper-cut_hi))
            extra = wg[at] - first - second
            weights = wg[:at] + (first,second) + wg[at+1:]
        m.extend([0]*nd)
        p.extend([rest+extra+sum(weights[:cut])]*nd)
        for i in range(nd):
            y,r=divmod(41*(j-pm)+28*(i-dm),57); yp,rp=divmod(41*(i-dm)+28*(j-pm),57)
            ac(reducer([ix(i+u,j+u) for u in au],[ix(i+l,j+l) for l in al],aw,top,neg))
            xs = [ix(d,jp+y+pm) for d,jp in zip(ds[i],ps[r])]
            ys = [ix(d,jp+yp+pm) for d,jp in zip(dh[i][:cut],ph[rp][:cut])]
            if split:
                split_d = dm+ce(2337*(i-dm)*dx*13+1568*cut_num,3249*dx*13)
                split_p = ce(57*r*dx*13+1568*cut_num,3249*dx*13)
                xs.insert(at,ix(split_d,split_p+y+pm))
            ys.extend([neg]*(len(weights)-cut))
            B,N0,L,R,sm = reducer(xs,ys,weights,top,neg)
            cc((B,N0,L,R,sm+extra))
    if debug:print('setup',time.process_time()-t,flush=True)
    if initial is not None:
        m, p = initial
        if len(m) != size or len(p) != size:
            raise ValueError("prefix table extent")
        if any(x < 0 or y < 0 or x + y > Q for x, y in zip(m, p, strict=True)):
            raise ValueError("prefix table simplex")
    for it in range(int(initial is None), dep):
        mm,pp=compute_plan(Am,m,p,m,p,rest,Q)
        nm,nl=compute_plan(Cm,m,p,mm,pp,rest,Q)
        m=list(map(max,m,nm));p=list(map(min,p,nl))
        if any(x+y>Q or x<0 or y<0 for x,y in zip(m,p)):raise ArithmeticError('sup simplex')
        if debug:print('layer',it+1,m[pm*nd+dm],p[pm*nd+dm],time.process_time()-t,flush=True)
    return m,p

from functools import lru_cache
from itertools import combinations_with_replacement,product
def root(m,p,dm,pm,logs,standard,half,dx=800,T=-11200,dr=800,du=800,kd=1,kp=1):
    # Q=message probability scale, dx and all innovation widths in 1/40000.
    Q=1<<logs;G=40000;nd=2*dm+1;np=2*pm+1;limit=nd*np
    if len(m)!=len(p) or len(m)!=limit:raise ValueError('root table extent')
    if not all(0<=x<=Q and 0<=y<=Q-x for x,y in zip(m,p)):raise ValueError('root table probabilities')
    if standard.p!=Q or half.p!=Q:raise ValueError('CDF scale')
    u_end=(4*G)//du;r_end=(4*G)//dr
    def mass(c,a,b):return max(0,c(b,G)[0]-c(a,G)[1])
    # Every omitted strip mass is explicitly assigned root/line value +1.
    R=[]
    for j in range(-r_end,r_end):
        a,b=j*dr,(j+1)*dr
        if a<T<b:
            R.extend(((T,mass(standard,a,T)),(b,mass(standard,T,b))))
        else:
            R.append((b,mass(standard,a,b)))
    missR=Q-sum(w for _,w in R)
    if missR<0:raise ValueError('mass')
    step=114*dx
    # All V bins sharing both upward-rounded child-state indices can be combined exactly.
    @lru_cache(maxsize=4096)
    def offsets(d,p):
        factor = 56
        end = factor * 3 * G
        points = {0, end}
        for center in (d,p):
            points.update(range(center % step,end,step))
            points.update(range((-center) % step,end,step))
        order = sorted(points)
        rows = []
        for left,right in zip(order,order[1:]):
            w = 2*max(0,half(right,factor*G)[0]-half(left,factor*G)[1])
            if w:
                rows.append((up(d+right,step),up(p+right,step),
                             up(d-left,step),up(p-left,step),w))
        return tuple(rows)
    sh=logs;Q2=Q*Q
    # Replace C++'s 6-category line table by cumulative tails on exactly 4 payoff groups:
    # --,-0 | -+,00 | 0+ | ++ . Missing mass is put in ++.
    if kd not in (1,2,4) or kp not in (1,2,4):raise ValueError('root lattice factors')
    @lru_cache(maxsize=65536)
    def shared(xx,yy,positive):
        x,d=divmod(xx,kd);y,e=divmod(yy,kp)
        d*=step//kd;e*=step//kp
        x+=dm;y+=pm
        a=b=c=0
        for i,j,k,l,w in offsets(d,e):
            i+=x;j+=y;k+=x;l+=y
            if i>=nd or j>=np:d0=0;p0=Q
            else:
                x0=(j if j>0 else 0)*nd+(i if i>0 else 0);d0=m[x0];p0=p[x0]
            if k>=nd or l>=np:d1=0;p1=Q
            else:
                x1=(l if l>0 else 0)*nd+(k if k>0 else 0);d1=m[x1];p1=p[x1]
            a+=w*((Q2-p0*p1)>>sh)
            if positive:
                b+=w*((Q2-(p0*(Q-d1)+p1*(Q-d0)-p0*p1))>>sh)
                c+=w*((d0*(Q-p1)+d1*(Q-p0)-d0*d1)>>sh)
        # Compute all three upper tails independently with directed rounding; clip
        # valid upper tails against each other in the stochastic order.
        # For A/B, lower complement computation eliminates a dim-dependent correction:
        A=Q-(a>>sh)
        
        if not positive:return A
        B=Q-(b>>sh)
        C=Q-(c>>sh)
        return min(A,B,C),min(B,C),C
    @lru_cache(maxsize=16384)
    def u_mass(start,stop):
        return 2*mass(standard,start*du,stop*du)
    whole=missR*Q
    dd,pd=step//kd,step//kp
    x_stride,y_stride=28*du,41*du
    for z,r_weight in R:
        if not r_weight:
            continue
        positive=z>T
        x_center,y_center=28*z,97*z
        total_weight=0
        value=0
        index=0
        while index<u_end:
            xp=up(x_center+x_stride*(index+1),dd)
            yp=up(y_center+y_stride*(index+1),pd)
            xm=up(x_center-x_stride*index,dd)
            ym=up(y_center-y_stride*index,pd)
            stop=1+min(
                u_end-1,
                (xp*dd-x_center)//x_stride-1,
                (yp*pd-y_center)//y_stride-1,
                (x_center-(xm-1)*dd-1)//x_stride,
                (y_center-(ym-1)*pd-1)//y_stride,
            )
            if stop<=index:
                raise ArithmeticError('empty root interval')
            weight=u_mass(index,stop)
            index=stop
            if not weight:
                continue
            total_weight+=weight
            left=shared(xp,yp,positive)
            right=shared(xm,ym,positive)
            if positive:
                A,B,C=left
                D,E,F=right
                pay=2*(C*E+B*F-B*E)+A*(Q-F)+(Q-C)*D
            else:
                pay=Q*(left+right)
            value+=weight*up(pay,2*Q)
        if total_weight>Q:
            raise ArithmeticError('root angular mass')
        value+=(Q-total_weight)*Q
        whole+=r_weight*up(value,Q)
    return min(Q,up(whole,Q))
def truth_controls():
    v=list(combinations_with_replacement((-1,0,1),2));g=[0,0,1,1,2,3];sg=lambda x:(x>0)-(x<0)
    for b,i,j in product((0,1),range(6),range(6)):
        z=v[i]+v[j];M=bool(b and sg(z[0]+z[1])+sg(z[2]+z[3])>0)
        e=[z[k]+sg(z[k^1]+(sg(z[k^2]+z[k^3]) if b else -1))>0 for k in range(4)]
        u,vv=g[i],g[j];pay=(int(u==3)+int(vv==3)) if not b else 2*int(u>=1 and vv>=1 and max(u,vv)>=2)+int(min(u,vv)==0 and max(u,vv)==3)
        if pay!=int(M+e[0]+e[1]>=2)+int(M+e[2]+e[3]>=2):raise ValueError('wrong root quotient')

def literal_minimum_controls():
    from itertools import product
    tr=(-1,0,1);bb=(0,1);sg=lambda x:(x>0)-(x<0)
    split=lambda a,b,c:int(not a==b==c)
    for x,y in product(tr,repeat=2):
        c=[min(split(k,i,j)-x*i-y*j for i,j in product(bb,repeat=2)) for k in bb]
        if c[0]-c[1]!=sg(x+y):raise ValueError('literal NAE recurrence')
    for b in bb:
        cand=[x for x in product(bb,repeat=5) if b or not x[4]]
        for z in product(tr,repeat=4):
            vals=[split(x[0],x[1],x[4])+split(x[2],x[3],x[4])-sum(v*w for v,w in zip(z,x)) for x in cand]
            lit=tuple(min(x[j] for x,c in zip(cand,vals) if c==min(vals)) for j in range(5))
            e=[int(z[i]+sg(z[i^1]+(sg(z[i^2]+z[i^3]) if b else -1))>0) for i in range(4)]
            m=int(bool(b and sg(z[0]+z[1])+sg(z[2]+z[3])>0))
            if lit!=(*e,m):raise ValueError('literal root minimizer')
# Width-proof parameters. The stored table is reconstructed in the energy phase.
STATE_D = 35
STATE_P = 60
PROBABILITY_LOG = 26
PREFIX_LAYERS = 24
SUFFIX_LAYERS = 18
INNOVATION_BIN = 300
PREFIX_TABLE = (
    'c-jn?hdbA86voeU?o^6WiIz}Q8rp>xN=rj~Y9~s<E0u_(Ni<ZXGE-'
    '7U6ip*#B@#l@mc6oC@9z)zj_X|a_}pU@eq(pSFTC3O3;nzPg0WvG&X4OvP;>_pCU>AMpdIJ>v}2%I8*1ZP(L$||)'
    'A@;Bds^_-'
    'z8Tj%nxOC02p!7?#OT&Txmz6uKdHg^h1EEo`vb?uS7OeYas)(|A||aEg>S#(ocC9x?D>qu$^snS{}G<|bFk|Ddo2'
    '3+4ofekW4cia5Rim9zvEF~{~Uwf#3J818rm`8*o9E!KfaCU4t|h&dmUWt3B4EYXdiqALrxq=!Z17ZEw_f>@4cY*9'
    'Z0ym5#w*Iz=rD!u}D!Def_3Ef5}*=eNlx`g#x@T^@OJ$;YoEf8{}4S+Ur6le#csF@w_AI7c|!WLXWmibU)jP^@BQ'
    'b{YD4=4d}oj*LG}cXoId%8@hP6Vs6n-'
    '6!!fI`B^QPv#uElo0~9WV<VPNZ$M4odK`LN3txvC_$gGw_recE{;Pn>ta4mGT8e|`iV=GHJDT=<#k#?t(VAF*X{s'
    'O5+a(7zf$y>J-a9-yn+|i06zqALgo_pN_%ZMWavNh&78nh$aS@obCKUKD7y-'
    'TgF(SqXhd*7x+&mA|Upb3`RnGW6${yL#HrTa$KOD?<Lv^jO<eyffFx?PQMf!L)RR>WC6Hsk90@I2H!&0#iTHi{+z'
    '@(GymFn1da50}-'
    'p3C(PDV(kLjE5}xC3&tB^+!6f5S^ei9hkr!Xt!v`tb#V28rueY$5uE!{fPmMEf}QOf|VN0I6AEf>cbncvv&h-'
    '7u3PuyB2Ha)?i9b6|y(~z{H0Y=qX)}rtzi7TU3k#2H$a7{wqEef5IV?0*q<Q!?!6pSZ($m&u!k}z}$4i^iIJ7yCe'
    '*$ibwC?&+)S@7OUc-'
    'v3pVkM%)X<;v2!Z;Omd_;WzNZ^(w6IUBK;~=WuL_3)*KoV7#?0a+@unYrO|AowwrHuyvS`whVKx7$75ZCU&<@#=n'
    '`P;VGwr???JWzPJapjsCG?w-'
    '!ELTEQ1PKk>eWnH;bsk#8+~#9cmgV&kq(v^I2L)PWB8RJJ2#Z9DeGw&6>^HdJqHh2f2#xRcz1#M)-0Q!_S-CJgFq'
    'fL2~TEbi7JX=5$Ebg4o4l`1^!^#fZDSD+}S3`1*6Ftw!!n;X7Cq<+B)hff$jyZ|!pc?hihfD_%`qu0Q9nBSQOkNe'
    '5ky)X&iuEoPI>p3!GW8vo!jVg@@m^+7J&6{BO2Kd8g?+vKNU&T<mh{*1)kiGAMkGjXupnnJ!b1d;+j44#QZ%0bE4'
    'Vd6zgiiH^aQ!<QpL3_;WAAa;x^y_gl?GzsnO-=S1spbO=P#pc`Hn&n`zz#f*5MTP4|>LfzK8Mj)t&fW(19<TI$-'
    'm@9o=WP<C%XO9yYYXYECP<9{h={J1w~MrWtlQP3S0WL~}v|itg4Ud~Y4{)NA3GQjL8psxTqD5@o|H5N2D36P_io4'
    'lTm*8{gov_X~PV{Ddx1`EZ|+2hYPFa6Tdn`yXe(%_9x-'
    'HYKBX<|}+MjEB?i=eT4Oi+ziuG4NA3OiV&C+BO(;!yn0WZeT~@RZJao5%-'
    '6=;*Z)%6xtp`OZy=x<yxXR$_y1sJCHnnBl^!=iM3vfkT5_GV=||sz(*6eveY4;HyC<ZeXz8=EB@{O&GGdOJn~f;|'
    '9DZr8a^3pwm*Tp9*btRH+T4JN(a>DccAb|JEDiTqn~pd8uD8~gIiI({wF$|TX4s>8PcImFu2o*-'
    'QEpYWKoYKojUBTu0hJ#YB=<%g5~i_Y%M8=^3XEOTvmcXYm4wn{~KQQ{(>V3g&43kAI7h8Q9SqqWc9P)Z;$~}PD5w'
    '*8&q4o!l+kq_^-'
    '=z9F~qnYJL>n91X|8%6sUX5saUc{ZT1@1DivxVpIDC98q$`iB=aVJUs^Wu7|N}ffbCD_hNPE4h;Xf5ue*uBD-'
    '=i{!Y?|(Q$2@G@JmRf)SYiPYIL0_QSkvX<X~_kK=c?aQ2`doPOgA|4Vz%-'
    'Va}Mzw1wU<CXh7rPY_8PwaqHU_0b`wqv4s8y3X0qE)owz@(oTzr6*{hnq3Pxd}<f8}VmH1B&NLGWT~a@*mcqi(xg'
    '2-'
    'u*zVK_zt0mt*dyQb;uvqdOO&pzJG*B0eK|Zy~(PKVr3ME*4$ShHX|RRMXz#xPK~?cfNu0$3zTS6^BOiXISGLgOD{'
    '*(ES-E>70AGcKbH63jA=@%?HhMuc9H~0!ot4;b4Iaa(5j=Oy9$>8fb-'
    '0)Ayn(ZwE%}Y{D+DRj^xX2=%dZadfs0o*bTtQk_xAzdjVXZ~8-'
    'LS5IuW5Ky|(%6q<6^Xl{8`9PmsmOYZnnL01{t!5;*&kN!YZ?3VCdpiQ#+F-'
    'G)4dH&Rxb^iX1`Yg)=QCSya7{BRwl?9?ibnLA*#OlZ^;nx#i>`+y*&<hs@M}L%)mnke%gPb&Sc(TB#Ta_$J2c(Dq'
    'Qv+!4tEydu<J)GtjIyVRyIOcW@6cjw`iS^icfD}W6+jFYzTje)VimTtA31$2@i3>E)0)9hQMjeZQRlDgTi1Rq(@('
    'Y7<K`j2hO2-p$p1Xj-'
    'k~35Z?c^#GnQ<bUfOD>TR1~>$nOk!G@5TKNqp(I`EXAg!x@HpqZwO__F?(S|)=9p$wOk?W_@7%L-'
    '}5tP=W>ceST;$G&*JJw2LVYTxDMSAE#8)14c?wBdzz8x|gK#ru?>=-$!-'
    'SLGH=8{dq&DNS$~*N7L=4TvqN!@tY5FqvP2n(Qjn8~s4avkI8}DT9-'
    'KDOx05RJ`^(GRA&|Ldz#q`xanq|BtYC&B29-@8MYT4!_gV@%csyPAz*4?UxDIH}NI<&3OvnwU0rA9-'
    '`vq1EkIi!FRJDB>DToYxQ+B4!HuWBOb7dI17!`6F6Wl>7#LnpkZW*n`g~%Y3UBUF5ZZMt*cP6+z{PO=fW^u2O&n2'
    'kZq)a;(jW)^Fjf#*)mv_4wSFy;2#S0EXJ3z-'
    '1GttnfQ*Qwk7f&sTgj)caP7!xy7ZUF7Y3`ll)k@4KplTvGvYRlznNzp<m5-'
    'E7y!+@=Z9(jnJ#EhyH^)M3_tNB2$A)zE!9h^aEM{RUkX943n!$@Tf;IHgtSLVc{2e-'
    'TH(LQwwnAbsqNj$brGC_jt199Y(KB$7Qt?bdj*jf#nHMJ@EqiPoF?N_z|8Tjl>)E2S|6ii^-'
    '*d=#zE}t}ndtVU;IlrMV+P?kpNRow2gT9%k!oaU<gZ4&<BS(#!2IJ-iXSl2@W+gdtSE>0?5q4$^y1g13PN+6`1NE'
    '>!`ion;}h9Uczo<k`3DxvN?kM>rI6`Mga2wlay^EMqyj@jj1e^y7vkPu^j5hIM_8viXWu6nXx{&Sx#CEN(_lbrW7'
    'SG-Cbt23(D;hqh}Sj!&sYTTwO4cU8gRYb9pSs6hLnGOT)7f(PM6I2I-8=%Zh-bj~NBFdtb9^N=3-'
    '0rr_$xY3#ckBl@p_#|VIRuYc+O1enu1%kRffsgznguIS`llgtD&A)?TodMWA<`$~*z0ewP8M=+<;WXn6F6cTVqUk'
    '6?qYmPJp#_ePGR3X>ZFn5D0Y#%%LQi)ww5|2w*Q$-'
    '=yAxsZb`*vdD<kcx0<QFz#fLCJ*0h7K{i^44)@9t|OCeXNXR%ju5{JKg!nz3K2i^hvWdBtjx8xib@G-'
    '8_wqlb*lDmYqU`awVo@6&+b5<jQo;6_e&3ZI$sKYqLTEvG|W4?A3UWQb{_jfrQ7nfn%_7arZ7vb>6Z?Ks91?wvdk'
    '#{H`JAdV(kJbmQkif$+n+)jeNQ2SPWDI}s3hI;MF?+*vY(F21gsst-En(fpgix3&-'
    'oZ+<09;vc6YKQ6p!)0*tfriY{Jzt;Y<wKiOOE1AxDBecEU<3U9;_R?4MD~0F>me)Yzbe4XZCuiw$qld+yt0wkHW*'
    '9%Gf$s0ad?cuxv2Ud8&i!^y@kPNExs1`H5}DW^sVaYi=I=lui4EbM~4*KJ0gmqa$27bb}+`w)v0mx|=XuTfkw>NR'
    '4QMW>h0q`843Abv;sN)Zu7R4bsk3!<nlPZ(E70Ddq6&Qif#{N|2^jgzBN+;PB})=3gtsbCrDjaLYw(SvF>MW}?4B'
    '1~R{-'
    'VpqT$<Y>Oa`73d#D|!a6>=?95Iy(DcI9mFJ;;(5iHed9|g`gW4J^va`OJ9QFJU5hII)(2pPH<mnhnDL95Z!k_UYn'
    'R;huT(XHmt*$Ek;;*V<8Hj>tfkjE!<U_fOpGB;^*0+(6Z@|Qe_#edCYkHu${9y>e%vfDJz5(@>%&TzA@`Hf7tbue'
    '^-QYYhEBn{kq0E^WC`50VnpEW5ef1ney*tn^@JO85Xyi5bN29WA+XBwYDBpde^}rr3MSvR-'
    '^jc4@{DDbcAO){(UUP=Za#i>--'
    'L#w_h<j@H5UB6apC^5i}<kUv09{8IcKV@3$D|l#1ILZx9=qh_32!P|<ycg?nP)Jt_(*4Pg*_@1f<<ZNxSEK_}S<h'
    '`x$-'
    '11}=|k}C#AoWv7<M_Ah*L9egY7`WXWu133|vT_S78rQ;U<8m0jTmTubIf%YF9a%IUi;s>#>@6jXbdX1ySx?*w7jT'
    '{1&YK3+@xLLZTq$8YEr(3*@h^!d=04$!p5feiCXnw;_Tq9+H=gM0#7}ULqfeRg6#vcK{r7SXaA?95`$n8MZb0{m^'
    '|;(zi?i2jaBX@uF5UZqr2{K*U}rf-UoXX*XT>-'
    'a{vB>UU(s*rXWXYktUUS={$F!&S0fwlvobN(;4L=%O+noA*SM^eh;dh5;%(Mb{7QU`agh(PL?;Xp86oiPaT_Pb_#'
    'sil2k!5$V2t$zytsZ2)<rJJdv*-TM-Ss%4{Pi++=nuUU63AYjK;P#Xk5Jv3%uqdHf%N;Os3&%e@$$gtBxsGhQQIT'
    'AO33fz}|0!&_Qi{)wq_2^($et<nxRD@0bQ8vBHigthp+TgY*J<SJ*YqljucPau=D`HvFo^l()xkX8kB5_82mcWw$'
    'nh=QUucVm+)SJ9@fP4NAIIBf;$l@_$r7R;wI?txDl|v{=Fg-'
    '{H0VD;_F;hGldC0;YaMqgM{xlHa3i>pQ&7OovWb3X&JShId*549C7i@XV)h+wd53dOSp2$^$fNgkZQu5T0J~h3%H'
    '>@X)>juX7$4^7|}8{+&QawFAzE9)farE2uf{Md7xcSh;00{3=$XZ0%BPylntu`&l@BbSg@6$D+_!4aJs&@n%b3ES'
    '@cmN#6f>u68S5QmNr`{bEjR_{duwGq_uiSDd^)mJOpGu<g14J~`_eJHB+~cdAZ2{)P=dbv9-'
    'FZxX%O8*$m`c|4%E7Q3o7;6h6s-'
    'bdDA=js}$*Hq!y)*o<wRe_GK<rt=2inJBQ(3$fc>qmS=T>d8n*cISeMIPb}auDYH9@6LEfzPC)kDi2Gs*@ypApzq'
    ';Ucf5n30B5Ef^I-03^gA>A@D9fqy%Dg$t`qA^M>gzPsFvjWADhb*t_5aW_EW#_H$bd)3U@<3p3by?m)-'
    'tO(<?#g^lx;pu~9|&Yqfyo?25-P(B7y6NaOG(;(3JK6umC4fn46WyP;8+%USDy?Pb#g5o?rz2z-'
    ';DNf+F$ua!V@jmaq<j?DuUFD{xbL_F*k!_Rz<3eRqUOZwmtFAHP6Pokbq<K2OzBPu+KGs3TtrogEH5l`}3TI~iK<'
    '1?iO#NDhQ@u*jBwvg&#qa3&`~^C9K7l6`pfEHKJCTFTG4HW***lyXoQ@lf$&j{BLUC(63g$^VddU-'
    'f+5QOq5Q)Q2?qheCyNELgMB3_GD3JG-yzeq*4|c~umowO9;fw$ydr2)E#OuWe&=zhAb;s>EWVI1XD^_Ca8bbuS&B'
    'd;*GZ17n8Opw+(V;#J2UQ1RTYotyeea4r^Z&5^*Jf5u{K1}0-}t>(E=SBy=b>-'
    'oxuf<Gw=NB3q3_2p7GL23MQ3^Bf5&*f_>W5m?cr(HH*uQPa&Eadmkq6_^Ss$(`1xpMj@?^}aXo7A^+px$NK10pya'
    'FYWWhl!nL2+3Tq_e)ESJW5eZvKSno%vASn}^iLAD~s91*IPuFfB+!`{`s{)kwnVYw-wv_Z*GovG|rAjmW?VWDdQL'
    '8|Uvp@nrxaAKk>UlU_J0e_7(G^Prg1c$$73{jMIxwjl?R;cS8K_xE60{5CAzw*lXzSE6I@VwfG)N9i6N#7US=Ek*'
    '-ZYE|GcZ~&Hd=?xu!DMXI>&6lj3SoV1(pSt>0!n8RYshY+;9O8I&N;H>Lgz$wfzC5m*CrhQAVU;5eyd~X=wI`eK{'
    '=FMH=E*YN`9z;%VyE%z=c9S&p`jdL*q3ukC3|s46|BEiV)cXy3^-'
    'DTZhj@$9ascc|8Lm1{R>8n`Xt#Y`53L72M_ZP=yx#-'
    '_IEPyY+o9@XC*^6_7zme#N*oT=dj!ti|&@u_#qX6^U<Mb_z{feu>rWHb`zJsUBl8#mtdFT2A>h9ak$M13$pC+ZM_'
    'X@LiQu9%>?%!Z^doj^=STWgvOnV5PDn>UQ4tQHGU#e{u>43PG$VfQ-DRPEb@+ZK}2aMd+%%H^aU0C_4pTlEcb!`C'
    'a18}jF;@XFN$?j@AA2eH~DwiCB9&KnpFqcbNO9MF8r~Ji)=RV);&vkK(Zcp`c37C2O9j}eI>Rt>BIkyNb#QWRk(P'
    'q5>q?MF>$(tbvKkCc3%-'
    '*F8T)55nnJPwGbXV^KmvU7mF1?;G$s`Zfj*=;m|bv_vQ^wZGVOAmvLBF{|ugTvB>)vg_C~acqAJN#Y4eZeAORaTy'
    'LOu$2C;Hy$D%7H{4!z3Qvrk@U6EUQZD|7n^OCc{d+fRe{VtB<#nhZXM`XBErd9ti$nHWkZT!_W9}odENUqBed>>w'
    'p)$xf1f&}}_~*ZR-'
    'l|p3MKYh+vEe;OIwx~h>2v<H?jeUPzQfMDZt#qW7dhhJNmePa<AVziu=J0eY}T=k4c0H=W;0zraAgW#3K+#}PYvM'
    '*cjb7>iY{EWxLwpQuf%{G<#<z43TIJ*RVc#IUtiHJ_A_*j6vDXvBc#o9apQhAu4ZLo$fvhxzMhKyrf)FeV<L<e#G'
    '#kPGniY);N<Qogtvsj|LQ$_i@OaA;g67SKDgv{6>Fzn#Ls=MD1LSliboyMzWj(Jr>&tTZ4T?%yV14L7!%&DMO6Rg'
    '@U&k5|ATYzT5CE!NRNl~+!2V3SHjbi^5}WCC-'
    'y5cRy(xw)~Y(rds)hVnhSY~Wfo^oeZ$r#o^j`k2oC8P%rib*=QqtByzM_1_R>GXdL{e0+teMr@8w#akhqw8&YZ&m'
    '11Iy))g!sqU@$LK@6BBnf))R^3cnH6V${@f^xR*Hob$!#>HQt5hrVLX>d)|~D!^=~kLX*R0~5_`d{~!>CG*~5!uV'
    '9oee)VSHzgu1^d%}jK80QFV}yr2#HYhyXh{#jp$WIKceWpL$NAt&!xhwSxPX0m=Mb0b0;iN?ShVLbK2=!Z#JqhNJ'
    '!uy*wT<y0WDV-ImZ53ad<@?;8@-'
    'lHgHPNzsG6z6YvT}HcIbymGkQStD`Ci<Hm<X(<?#3tp1P}mH|uBevlmJ1F!U)09uMau)<L}Ofj74qxbvW)CwSY_L'
    '%eR1IUlUw#=Uy2;Y6QBeE-}mj!vG)7qZm(XXZeTYL(@NIKf?7ev06zAEGV#i%88WfmL-'
    '8)Ed7@cI+4EUjKwgs|(OevZIgc<)G;Bdj$EugRf6IZf;9~ulj3LoKL{2ju(jT_Y`SLj}h`K5(`g0!1zHSu$~nJwL'
    '`v8pL-o6{(9nriw7n>Ig3f1CvYIl0ft8p!QitcUXR|3d9!!oUd1N3Jzb4AbxTmT!~jc2&%(@KQ;-'
    ')p7KS;)A+J9eN=<z*rL#M(?*GSMDt@y5s2ZO4uZX9Q$meYrGq}y@6(2Yi%O5^HVE=ysZ2IaNR~NXk(#qq^>b5-'
    'V%wAsfYYVrXUd2XN7VxnCGx<gBcvhI9#%AgRxMqe7(_CT`=@yYUuTn_Md=|2^vP6wS5z>czLvx43Ls^Ar-'
    '<}VJj$9ZTeSrG)EX;V40h^#S^mj|fER`gfKZ?ipw&(b*{si)xkC0FkftjKAv9aY2f;R^uQPR=f=6GXs-DMnI?~cm'
    'XXD~U!88%_|xS(Z=vX}#?Br|-i*p8B~8?nh?6-'
    'quDV(7QIn3O*Qts#@qtT+Z!9t?xa`GFYuL=HMDyFvNqAD+6Sh3{Kb@rj4u`S0R9&YAd@3tSTTW#wbuKjl8#4fp4z'
    '%dhg8f^+Pi>d2O}ZCLe%Dd)^G=E&_UxT@ED)@aw}{_8Y3MsFC~tySRRHa&R%qQ64bsYwJZC>Jjd7m5;-'
    'cjC(4M1iI+Xvq77NY4UX)Bgyg$Q(Q#m5rSe4|&<W#S8OPMCiQ1vxY?6TN8(C&z|9&WG_xFjDkUEIO=AEVs~6H6h8'
    'T*L*^#drd`7zr%N#Hbc5ye)40F+I6}G}g=)49daGOD#r{25W4R4;Mr^>GuPg8pi_x@RA3ula;9%WEm^f&_>ZJ-'
    '?#Q;1D?u9eKQb?oUJnL{1i(8d^#pf%RdFJqe)>O_2i{sc2(QG~99#>5E<?6nkZ1U#}N0>YC;W^fv=Vrp&X%pveT+'
    'WuRbNN~qE&iQ5hL=oHVVyqmd}wlaK0fic7_QMEp3W~7d*t%PWTUsDb#}a{JrFIH%r8Vl-+XBL=0bXsL@Q2ND0-'
    '6t_k2kYy-G%|e-fTdNPt`13skp1!R&F5v24^s3@Lj6Bd-u-eGdZb`yqOa4<^Z9MN`BD{2k(oEp{j2waF2xD~>=h+'
    'ZriG=4f@<jm|w=@O|w%SUp;f(RvG^JV_UYW3<pzHXb*PMq=8#p|A_+53f`ij9LuX)^+f9?FRN+SkBFdK68dlHjk)'
    'I=I>)(aNv!H?4EXqM?98vQRYQ%TYZYJ>mOx*-'
    'vhj8<SzbwW<5JUTEaX*kK_7GWu1!}tZb#kYCe5<ik}p>J9mhq(Y0cRQ;|6LELYUrOcljXUx>9ekwUrcw$S^Vi_N|'
    'tP$c&r8VBEDu~#|*qf(%L<29;wCc>~+91?7wLGx`4T0TV~`)fEZ2ZdtaggaOj5rE2^n^3Lyg5k}}xVh3DJzktaaf'
    '>s2dOF}yye+;Qw8W!WGwlAg17B-4K|gpkELSgu=2-'
    ')%yqpCEi>Ywd8wc4KHP|%`hL>tToGX#W##qA2&Q|`jx0dA`O8Do&0<PMW$#tKTc)ICRzUCaxqu&Q{Q>r(&`n$8z$'
    'rIctf0&QOnzNB4bL~ynaAx=-'
    'uKhBLM~s}rpNEWK73)DP=JaAebLNX1TZR1ADzW?US7Ci1TS%XMBW{{M73L4Z#FI|}LazTcai%Q`eeB<1^V4(`H>Q'
    'A_USmdHA}RyokTl{s!tTW4T3$4kNZ4Pw?mo(b?;_kV2<o1`D6hT_q+EgGy$cZ9uIL(Y5^=STuzhm`(vSW_)X4o3e'
    'l)?^xUIP4w;qP4RzUUZBD~wFj}Sc__|BUI|56Po8mVGn^Z*Px-WzrPU2(kEAJz(M=82zwu(E{zmv+zN&z|Y5u{(h'
    '`ct2**^F9YjvSnYdtE}1Al_y+r<e^7xI4{wZmn}8s1=m-'
    'w+FJweH%W(sWi<Jw=`hZjrXb;$9_(lHSG?WcB>uXUi&+~B#apj;qThi;(R%%{==0*9F!^vx{G4}L_-'
    'r~Yo{mpPH~SQ<x$qiEXA^N?a~w3~p5tIZEH1P~<I$W*yj=PK-'
    '6n)!OHq)7Vf;|`(Fa58ufpWtMW{S<!}jH;u_*93oU@PONZvtAJ9z+G#+yOfVFx5z1@8{7M#cE0u+K2S#nxGn9W@P'
    'q8^__4zdE%33_<l2dFZX^i7m1WtB>uRs#(uxEy~!z;uBxo_MW#)NajzspR@mjhy35(JDh#u25&aK$iXvDG4(&nbn'
    'O85n6``0pI*<JAC_<k>2c3-Q+bq$2D`0O;sLIGI7Py_e=c;0?&h_k;Yg9tJeMnWI;ILepBJJ_Sfr3nxGjDac?$y-'
    'cX9fRvv4T075$DS;iROK>-M~a_V{Pe%Z)+1Ry0z4B48DMAGHa0G0-&#&YFI>^wbB-nyx}+$|aOavQ{bcG+^kAeRu'
    '68d2Wlh0hX{iwHMjNJJD`r4E2?3;WTMEj{jW%ucf+>Nz%f=^AiwvZ4_jbRG{dp2$%Q0&@V^|x1xXZ=w8j774(B|#'
    '(d*L*}42EU9uHdByf*@F`TY*pI-(0^VQo|`K!At-'
    '*R!{QmKP1Q*OpH<hQWfjg=gfJ)g~_W^mQA@jOaxIRCI$<bo+ZIdk_vF?vR`C|^<`mhJu|u1(1lp_Z@2SnU|0_Ft%'
    '|KjAC71zr{;bw+d@W-tC^T8dxtyTz2s=THlO0>^caQJ47;`$mLg=CV-iF}nkU1%WVuFK&8YM|1TR)J(gGcVpa;sd'
    'O4^B92QoyFE^yvBme(mdM((7mm`qkat~jleD$?n7177pDjegTs<6lr46Hki75N60ryB%EHxg8gGq9b$?S&x`F}a^'
    '+)th`w}ua%DCQRl`MfFb9V=dx<ch&l4hjutsY$o_@#^cm-qeF#O<Y*I#}V$Dy<ftzJJ_>mEx%JT<ZR11+%0x8|G7'
    'JomnIM9$<e)e^cQe_RGX*{sTP4j-$Xyp52EsNviR@BGhu!!T-'
    '1aF3cFX=#ARtWu}0oW6d$q?Ntvc1*uq$xX<s2+Y@@NQIs&@74`8}J1k<f<W9b}!l=Z%eG7m4%pUW7z&;x<X&LLst'
    'N%Smtg!L6WWE9%qmHYv$>1&2Z&v)RJ-DVtrxCY*0Id11JfM1Yg3(eI=(y@uSm#KlO6ROa&9|(nfIlP?F9nKs7@%x'
    'FbylQwYuiRI{Cz}fRbWJAD`}LX^Wj|wsrICECEtq%ozQO)B7nzrzV&{oR`FI~oKIpcK9X4)Yje|?MLRO!zDNke5R'
    'ijzmekhwQ@5>h6U3s~Fr|9uto%pR%EK(KoL__a1F{$N+m^<X52raxV#tyzNL?3q%J;hl}K5Z*PL-'
    'vXf%eRVquB*hP4GTp?`ApGi6pArx?;w9(AbS71g~IUb*gW<swjRF-gQsptdU6`^kDYO5sRN!?9>UJ))@VCn4r3=1'
    'm@V3duDlT{Ygb|4y(Nh8HNX<P+1R2v9rHcLV|2?%97$4!<|{??UfLULKX%1A*}r^U`X@We)o}iqVop#h;7_wNd7b'
    '5J7F(b3hpv%q86V8^>V5dkhKt;?<rKT@ILZODEcx=6U2OJf1M5{S<p(zUeAagw_ii4|#}$>iSJ!_0K)V~CXy_EPy'
    '4H*S$;IM+OrGd+D^0vGeJKp>9*X#R!D7$q>%w}GyI8fwS-d-ED?VJ>E57vGDz<5@5*_asi25ZnMR~7@qT`3Uc<1L'
    '2UFDlF4)j9&cS#Scd0?a5In4g+g2!==P*`Gzv+r$?*=7Ns$!5rq*@;sXn=!?DEmS4GNt>_;{|(Z|v;H%Xa%nOSj2'
    '(+2Yc<sF8v+eWd9*b4g!2O+HLa69yEL&(X(bQ${>ELG<?(dew|p%tk?-7z<xRd}95^G0mAt(9v$6;OQE=g$O-'
    'ER^rv)cg?O^#U>-gwDLsq||%le^HI95x8e~nP$Cp-J_jdfD&r`aJA*VT%ijo-'
    'y%*<4|FD@E7^J{JS`M~M5QgG9$vFR{6+o6ypB5*~eQMCuGv5xZ-'
    '%aB?&fdU6Kh+YD`Shct!R;bCH?siM$ZA|u=%UxV?(OYrhN5B;TQaJkMIr?xrZ@8d&Q^4AK=wfk_Vb~nWNt#~zK1H'
    'K2WL}Za6&h?ONpQ>4ykv0ul+T(Hm?no>vQie;pB6{WZhVkNVc>Cus4;|3T`<B=8)zPKAx2%wNU3<?Tvy$0m(o0tC'
    '8O=AP?(y1aUyht`g+~XUW0!}H{PV93tH+x0_RCxN(yUcHG;RUs?VQPX4o+aD4mIv58^9Zf$+D56;EK)`aY*Ke=pF'
    'r8Tz#D-vX&(Yxe2jCxhYgAMEZ&Vzsusc)@c#C=BTh5dqAj9*(v7et`jpC8Vae1IU+A^vMAUwN-V!WM9eemBi4-'
    'QDt`XyBw4ivY8UQEUUCK>E;}RRoIMs@wZ-'
    '7AmXPVO7n8hp!twiNj9j`FUG$gZ%=(2`Gei$cY1+7EJPGH#N2B=bFl4_Pg!U(W@oj4lD1H`rKBt2_W;gI`s|sE#^'
    '_2^<awL6`&W4@|ygfOFLoPpHwabD0OWvEuM!R$S?-'
    'M+1#1U>^ZoyVkJK0HUJu8_lVY{JveEjfKj!x0w4FiU<V|yR=T-'
    '}uoWI9EXTAj%JS|ol;c=nJ_s<<%Zg(&V5Df*rc60R0r;$}}bu`}IKDEa&+Jhb<SjVm?@2bJaG-'
    'SoNQ(Y)zGX5SdmB~4k_m-G{6LEXiQ?tiGPsgb_CEGNhJpXlJg_vD-Jg5tf$U{Z4kujg1n-'
    'D)q^p4o}uZJTjp)*5)cT8f^#=41E0*@*Zt4Q|rokyAbbxAKM}!%zVan`Kd{FNF&SezUf7Ge^v<;!5=*-'
    'gNFGPwM%O!^S4@iIq>;d{YFMQZU=3`|yZ?7kNg|Db`6m$|dhCxifD!-'
    '*?=|$967f^_aO#GqiYv^H}zOtim2`@_g=fcdikCM2mleI389icF)fjMdGb^^fpfDERGT-'
    'DtAPp-F2~gp1bf&J}%a7J1B+?G7|%?Zx(UxM&hLEJaJ>7mI&xMR<v|g5hIt#i`wPgg-qFR3Le%--'
    'fd+xEVz&Yy)r5BV-'
    'j6acuJc~!)dOmy@WvyV*J1T=;>|(zm_eyFnt{^Npx~icLC%ZXCv?aG(7LCiA~4U!3jg4n9~o11wEkVB5-'
    'X*JFDc?@%52qyzl%cHb4EIXGAA+`l=T^YGM@E`Q7C>^;@iMb(ys;pW(x?4t)HGHK&{J;Y_iaWv>|V1r-'
    'BM*r3hxo{i&C=T+Gwvp-'
    'J^kY>G8f5h2UjiRZgRBTSp7v~4R6;Yq##K;*@V&#WmF>>v75fpu1Ofonw{2gsX##mF~cy5yzqp)1WN9l{R3#N%mh'
    '0$V2fRcE7w~r{fFC}u`bx_%`I$Ca3Ov9#pB+Xapq}CEo)6^eRNO&mC2=%8MkFHVcI(N#?wZUMy{WvpVH$GQyhVOy'
    'Z$SN_!K$?pHc^zmpPk^iMNW@GWibpQ;INI0)8s&tVmR2r_ui=^Q#r#J-'
    'pBKAk@WNiNxIi_Q<JLalmAV1!HtZUk|8eDLTPJ=r^&oFQX2#>&jQM`Qm3%{YK0jKa!&xIV*?r0|4v6c|2Dhbo#Na'
    '>TOhkjwnpP@uOg{>>-sz(Mz?Y)sMx>Y*79{4wOLSY{DsFE-CZ;D_3GMV<;^W-'
    ')qI$ofnDBPC_&H;e@bpj@FINo|r(ep5j4%>XHZ9cTT1n303uzd>r#G|T(9>?u$mv)FMau+};~yW=I&+DdyPl<j)l'
    'T&9-XV&+d4Q(KSzy@x-B22A3|;Njkd`%sN`xLX>!-s?RTKMrs9{*tK%CL;4RY-'
    'Q^^6X_`nisMLrVCzZ9aeRlfjqc6ZjuL<~g?a_{Y9mJm}FS{&MydUv;tL^ZFKC-npG?qSx@L$%}Z&{h1ud6ZqHm;e'
    '7Xz0^j%`&4H?a!~@v|F<Q4oM0Cp&Q+lO})Ai3paB7$+{opV3&v=R)%hO_zrk#jV+b@>NY!g#!SBmIv^Mw_(g~OP!'
    'B2HUbocHW2?E7>TFFiZx{GnPJvZRPwJ#uMEMJi1We@QzLO(q*c=<ga|8h-'
    'W)4N`R@E3Fe0bJ31Yowug&`+Lc7;Z9m)zL}Z=*3gT$d*C4P(UJ5u=x@9PIxhN%{!a@Z8^$3=cQ_RE2jFw943u02_'
    'C@?;w~8P9>hTwz==h#1<zKU@Wi0Rb63T-*Z*kkOOT3`ZNe-'
    '1c%=!6y+11gQt5+ED?`VD2HJ{2ZPeyX~{6U<PD8t$GPb|LJD9&Ff5yq-'
    '{V*8vFA%E$q_?>=V6gA!w7x!Hhqb*K|`i6rdzMHAAy0}raNG}z$j?58#XHODZscM1>6h%@@53#!T53SkJNK1-'
    'JX{KBOnf-lBa<>x5+wn0Sd>u-'
    '|qy0(s{58^ccc<3|CuvQwJ?*|{L*{knq<nWbeZ8`UZ2qpJZ#gSyewHDfHlI&cGI})r&33#!xgK#%%Mj5#55}%KaC'
    'ko+Bk!nT_|$<oxn36jlLcZ^TUf8Tf<H-p;_>Mj95^AKHLgG8Q`dv|VeVBP;C6<89Jb?wQs%62)|fM!ma(RVE?-'
    'EU$g%CixN*HaPYaOZgoIYn`)`HtR4Wk49%({8>X}%r_&{iy+!B6K9-'
    '?0PxLD_9E#_<O65a`Ggyq!*V)8+4(f{xm@z+~PtZwQp(xX{SJN%Q{9#@j>l+UzICX0@2N}_dAPiXh(FdCE=NO^a?'
    '$#?Gs8Yy#{aQheq7}(PH=N5F^Y!5XY*+%nQHjquiDpD$2N`EaDQlZ*h(%n9r&V_1I$*!psxojfM`@RW#X0Jq!@<J'
    'p=&4OO)B)A<NiN}$HarRJeB*+3jzkagj{Yu`IT*x`^-'
    'm>A5mpt)hIGdOHv4YD*b~1G4F}~J3rez0L%dTSg3VqHUGMR@}4r9Z_zFcO-'
    'T(qi5Tv0C;VRi3?No9hd_((DJkH5I^<)S#G>?{tOS&J7RcZl<sR*4-W=7}QxsUmji2r(u{QP{gki=;!p$V-'
    'x;`yYL$Ba?E-'
    'HYb@1#yzKdvXS&^MlfB{xj`zMF44NtXDR=b6FnGunC>01B(*wIG9I{tG@UllKGQXH({DL#9le<9I_J@OC0)w1pFy'
    '{Drjw@pWNKeLp2kIup=GN^k=qz`DlS}uh%`g|n5_r1)zk2<X*AY#Q-)b&Uw9-'
    ')p|5K@FAl5XZuXyf=gSN>iF(ONbHn-hG+*8=?ZKP39OE7l`#IEgGn1wv2j^*Xw(KZ2u~lHBBVE`%v`M_F{Vw|S$P'
    '~(ZUx<H~L&dS-*TuE?Gs4H|knk`u5yt1&h~xbXghu!j@%*`(=q>t-K1EW3CbrR>hDy@i`iVSL-_h*E1d99eh^C#n'
    'N7`Ne$p7gTa$M?03nrW(l`uQH|J<4)r1z0qmtC}Tkuk0AxsC>DtRSarhScOZpImS0(VO10=tH~~N#CDBJ~<O;^we'
    '=QF-L<A`;DMA*N2nbM-_T9VkqU@9!#dL%i%0zfTCVHNDkJ7T-'
    'Y#JDJj4}LmDc_e{rZ!Ef3H6%62a@S@lmGKYSC;=k$EJy!kwrx!CiJ!+ZJ4wDo-UgaHq^Fp2Z+mAO<_hLio;#jn;<'
    'p?KuINNsx|G|lb_+X>gi@yZioldP3^H*1Rs<|V?VL|eQ%JWABqDvJAmyNaMWZIsuglI~nDpp=qyN=tZ2idGNlXks'
    'v>=J?R<8yD$%;b}TtbBty!vZX&37Gz##LW5hj(&XJ6Xwt$}^vr1~%^JUu97oJ0xr4LGxu*_!{+voqlP8g-'
    'lO}Z^G=`3T8%dp+YSgz=m3X-_-OU(6cTES;iOB=Vf2ac8SCyw9B?~d5aTYcmoP^;!N8o$)Knw_%LAPH2m^~ZV<xv'
    's$o1V==Hj(@9h~$h4UpAfP&OK6&^6hF<_CK<Q-<;56JM(dTscHa^k?+F)wbY4GK?Rb{pCn|)g^S}=-'
    'oo*Xi#VuZDSo<d5+*kmh}el!MB+|W5x%vz_|))+x_Q*nrKm6DR{oBxzr@ptTT%2-'
    '<qo~@^PyD#3$!5c6gB;|rwMCq=xld$+GMkfe771?zsGB-'
    '%H4==CoZN^O9Oh~s7vjsGf4f^bo%HxnPQ%dC*A2|DWQH8r8cP3?9O4dYmo|_j#r|XRfFkspMi9wQ-'
    'O~Dk*6q?ew4AP59uWIrokh6kzAf0Zf%=}6U#@#^12c(*~ww`c}CHY77i{b<rPD6d1y!yf3uBZ@zS4-'
    '$9u5lVtbZ3V#-'
    'zvR`aW^b9lezXx?ls&jri=ip1awQPS(ZD4zRN=&lPCo6oz6i9Kw^=^op}qkzTY_o=DE_K2$JuGC9}y!=IO8ddcAV'
    'gcpUrPA1i&uH$k2lOq(k9OYir1!_qQiZl7^}l(LX1_D1HjiEO0mjrka4ijYT241z7g0*vTnfycP0ETo^gT?Xji8B'
    'SlRJ(!t{6>yM~tMtebwm60#(YtK9nMKhtTaG1L<9yBF#wbPwj>MNTY9GQnZkx==@&v$V`^%|Mn!u$31AsL1|Jp(L'
    'wnlO$-fJMbIjFjPUA;<d1Fq$fJ_Cx#siwPs!|l@-ffQ3gVap7ui3_fsc5bv2@^SuA(`dQ80>iRQvL)p1;J2btNJ-'
    'BwcV)l;~6IBYI795p&|q#eljsLeFcC_|tu~n5*7joL(t}qFEie`hTYCtaMVxbJ}3|fUY|GQfvGrGOsvA-'
    '+S6qYLgXRjy0u4<G0bIo$G1t#1-'
    '`C`eJ(JHIHoCX49_k+SEQ_D*Y~>Ks7(cQj38GDL1K8>bGH}{YROq_bE|I(;%w5J%FyBRG`GW^5l@!myQhRLuXI)C'
    'I{(WWOPf0_OI?ql?pxRa#MFo%j-s-@4HfRz7);#o`hA6>gf7v0G7t}fd1s)Y@AZV-Clm?A79dW$egE~mVbwbM|!f'
    'uYbUOgG3QlbYuV0Omlc0%@a=%U-'
    '0fGVc>keTythpgi31;s#lhZU*V5x++)^{~D}IGoTBIX3T~iZVHue_Vo_5leqZO3+FNb8>UeRrjD4K9Oh+HpRrQ7~'
    'zX+VMlrAAtlLbxeq9N$V|UDi=Qv*k3>WFbB8txs`dW>QPVH0mKYi4ORWqj4b`)W{<!_1-'
    'YbyroPQpNEj)hCwvFQ;{y__9s<|HXi=#Lpv79(W-a7=;K;hnk&)9tS3E4;es?x`L8<-'
    '+TV?g4|S!Xu2QriybG=Q3OYRmbmRVLm{}>oGqg88lo94fG_uyIBA)O#i*^6J<P^*MtfTAA{-'
    'aLv`VE$Bm%fqfoaXV>JL7n}x&psl@K=O2m5Bt43{eprEz&;vh^9d&#JuNwMe6&N!hXdJG0{&=OnKHze0<+Qo1@C;'
    '^@nU4RTEEv>msPv3qOj}zeK^mT}amNFeO^-'
    'r^RwRsWo~7c~4zI*?Sj}?E!r{D>sv>R!pPweJ0W*<FTYUY!tc9QzO+h74m<eMDkgKXtLe_8d=q!=4ADwH3fZWTW>'
    'iUcCZ(vH^`9ng`Tv1dk=atRhk|S?@r&9y3vP8T`6O=6s0(Ip<ypTodRUMn#n#)P^Y}0;Lu@M|6U$##$8d_w~hads'
    '^CAX@>nT4iA$C|<cTW2yl|W=k2SXCx!K$Jb--etojsZRZ5zTvth(@R=^6>oe-JKbpNjR){^F?T8Ie+UKqy466Sm*'
    'xh+f{K#9O_-V%ELi^vbZ3WCL?(``ScG6OrWC?nf>i7s;uo<Sz3MQLNiOD*myJa)zxV$N9@>*~<Ck-'
    'ER)Ld1z5g=p?#sI*yLq9Ys%etI<<m6$;W+q8CF4(fb{WBv&m@T95kD-bZrOrJ@&|(3hpXuX<9*-'
    'X7F6PMQYDcBg#<y3ud7t`xaIid<~FP`@Xjol>CijZCke3+gvY(CS-'
    '6`Xh*j3?WMO90b*|GPtP!hu3VZWsNVNd9+zNYwAAX%0<Ea_3$M=RVT@q)26)6b0wpE2Cu!U#y4YR*k(hEnEm>b2q'
    '{byKA%HGbHW8-'
    'S9VC8+PqbKpSD0$JsK~<k_L!|Vivyfwe;;v0cC!8Lvk(AWUU`WIgy^!d+8~P`*E1s!^~;p$?ddl<vJ?txr}0i45&'
    'n9HmNL~PGMsvl3K(VGJQLOM(rO)5hsRHMbE)>w@Hy=CibU=iG4}sq8!b<*^62}%h2OxJ;|k1nsj`-'
    ')7!({D9*AgNjXZ<ndmOGyarS^4HW6d^k1zYj}3wrW)a<8Nn|~n$U~DzQe{*%kZ6&AU({p(3${pl>0b%!muGXvpLm'
    'v99LDV*uJd-+Q@s0-'
    '6&uEE;$>3|IB19_KX*{zT#1*4CYA~54{t@|w<s|(>$+H%cwDH$Ozgj4B$BhWL}8e!xHi70h>Q72Dsn|M;b{h8$}{'
    'TO971=}yvZT;9F3^3C*|Ij^nT+m`Vp{!RLYi<kK6+KD$$Ee-'
    'E=Cook+dH$58*{Bj{FxDybAGQP$K!wD_L_$qtt%!+<`t+oU(?n#)q>-'
    'JaAium@Sjbf=q=n{40Gm7I;GNZ+mteGdio2k6#%rtU8VEgmnZR|t`#6OqDBqEbVm71~4vqlpYvh-wvyPL=n-'
    '!wbJSIH;NrWfgMB`c$_0`<SDX0=e+52ag|kl=t4*%~}r2dB}M!HXEtRQgfu)<$HrzDwi+T_K6c)q=Q9<+<7r@n2o'
    'q_a+9dHo-'
    '49y$A~|3`im&le>7o26)hi{N3yeD(fJLLbl=mDKGQ|oa@v`uj5<h%KATb(FJr3Mwu)wSTTE*u{5t5bHZ_fzOar@*'
    'qZW^mG~~fBDqlC0R@)Dv-~ANnft)<$n)M-'
    '*zP;%z$<nm368~J1rqycQNwcIY9ZZs<MX6m#lz{>ULE=NEzPW;Qbp@SkBf9;HXs$QWTU#RAO+?=eh(1ptdN@*Y7b'
    'PNH1){>qfVX2ahYT&@xef1mRBIeh8~%WI9rxw~<4<zNFiT#dwULKEn#-'
    '#XjpYqF^1SihZ(%H7D&|i~7q><{6jG|*;?f)^k#D<43|C(+wr`#;_UWjISFO?_en$(f`tgl!t$9o5#8c{f?=G#g^'
    'CCI(v$WFwD0zBWklpJY^z_d<dak*Q-1Q7-'
    '$&^_%wstBdSd6C=S2QSMhZ@lfWqKbln7aH`q~_Q1G^V}}4KVLbf3;<)d~r|8yedtJirp#Xbysq~DMdL~yHHs;=xz'
    '~{$2g`u$%MPNPmtO$L9*pUX^)7Ox)XJ?APQeil%h-YYl7q+l27FhCVDPURCE3pFR`lTtk6P!s*uX#q91dFLm-DA_'
    'hA1xJ2qLho458~&Nmlm@qtVgeqP?4mE0S|q4JO7OYKXcTz*?rm%53Enl|EF%qCHEeXg*XHdg4G_7~%V{?U??RkY)'
    'L9`S`&H1n}!7q<G*8NG`X_{5oBEICLDrKaS3)|fVoUPU<>izrh|kKVax(@l>_ba>lXlFA-IV+&PD+DeJGo*GC$`z'
    'g@LZvAN78aaw-'
    'm8D%pJtf`MgV?4!#maP}X<1Sf_pA$jcmr}$0C}BY8Ys<lAV`qyDnX}u3A&v}G+*MU1ujHdCPX)v5lx#(<TH+Fmn3'
    'WcO0veHAJNIWX5Qyi!WMJ0`JzEQn_PUr`&zwumG4Quv)_`<A8zFQsq?sG%sBQn>(4d%e?<RxrDA9Qw<1_AN^BZ>U'
    'Dz!;F8+IDD)Q`%#EMug@%O%}NZQg<=(qo*N|Pc<uw_umfoJqx?H;Lq@uu*n=jcy_J*`PPKr1tLl4|jKs;^x}s<j4'
    '^OqxYkHcq42G2^LUng$KFQ=|KF$~67_V2aCCBqcX_T6nt;H7oQc<jauv*B-'
    'P<U7FtccBAP_x>7%FDcWiz(abr}*D|KUO-#;}f=-'
    '^5bn+}g3V$RXN+J4jlW5~nBJHh2MvI6hY7t$RWR0dIYs>}^&FV{pMKv5_^NBy?r*i1=82(}!#HSBmU`{&9Q43ADZ'
    'k7@6J*>^$M-Jnv$vyby-'
    '$pUYxInB8ju%@)?ufV~cTuBgE8Jp@#ezZuao$){s3s|ju{lCKE3cuhAM$C?r`Pmfa}*g33?%JMm&qpD<$w0xIw<R'
    'Z{rX+k=f3Ig4oN{pkrXT}5ET#;#0CryP_YHEMNDk$R_w$8yRf@kvAYBHEa%Ia+4Icri8*J^oM+DY{_*-'
    '_?|lRJdws5Jz1Nc1Xy-WyXE*mogjZ+eY;BFs*3FPPv;jUBgur*AA9mNMhTw6Yuo&!uIUgLb{DTcfCRrl5kr~QI8N'
    '=M%5HH{B;eoL(Mh(-'
    '3o}Cuvyw<?27iy@D)JSck!ie1pT!IudIVot+SMbzJ@cb39d@rze2GFKIU>HZO8UeZVns@Y?I*z~xdQA-'
    'z&eqh&%A;<Vl+v_=N?g%K#ku8t<=2;SO7q`+6&K$GrNOht%Io_+%J+CH<w)Z{GVZ|>Q3aosmxFf5g+&V_|KezgPD'
    'qi`fL5ZjBvd+0bQ9}HL$SU24MXD}Bmd}mjELHY(RbD(aLNL>G@Xof_9HM|vmbV8CnNtwJl>s)f%DvOJZV`E{Z9s='
    'RPjZ(x84ZPamV*{&N$M+9__|kqn?KadU=_`KF<g>>KWi!eLZv@t%D@BHkR(xgwI?Jv{|l3pGzuSbyXpJnu7Wk3XC'
    '@iln#P}R)W~Oz|&1Y&xrslD^gl<pXfDReSvB2Kn8QoN(&&!n7W<6tu$?RN_qKbmvZ9iN+rl8U)ec5OL^R+r!qb;R'
    '+(H7sw`;jrp(x`uN)urNnG@9OF_h8xiG6xCSIH_w|s}lulHRf#JH)*Szmc)X(uL;8uGf)JE-g5#$eY{Om-'
    '<klRB#~xBF~Vo;w~k4~L+jRv+wI(-pq&+9JMA3-'
    's#P2v#YzF~lbTz8ihe{Z3^h6uKf%#|e?oZJ}yt1w(yvlmwdK<z7SBPSJ;MA#<*xgVVFL;L~0ct0FWor>h#_3sg{l'
    'Q4rr-LE#sH#cV;_Cfp}2!5wCqUMqpGBY}h@La*sWuNll-ljsUuwWZh4huazPfBSJ+x$b>f**0ObQee4InP8Bs9I?'
    '+-PR{P6=mbY7$n;gN{IXH<t0~IM>(8Z4^aa_ua+f%MTr4-^#!H(=X_EM;t-'
    'Oz^C+=~T<+Yxv_%AJolif3H9eo*#%MPM-m(6JRaxu!MOoOHWD1?2_L<&1(pJ^S@A+#0R+-r&p9m1fytrpCW)Wo-O'
    'RUw}|(Bo?*L``&n-2xl5cDCfRG($`WW4!!rfUPg|P}M;fk-6HqT1^WUqKVu;Y7BE#<8)6IcAir3I9$Px8-'
    'hLAf(3qpZ@jy&905kn2f~I@BYMqwdd)L>&6z4dY9-DOeK?*zJlc>s$LEZ);L9H6m#k3^N6u6%-'
    'j7gr9_g(_Ja4Vob*Za#>*%R8+-'
    ';;Jo&PE?6Yfd7$H!#8??x$FG)w%gvc&&EqIB|VAp`FR$n!`iS=?Gj`gZ+<pz!;!ay^68Uwg3m;aXgIHxGYw^Dx&n'
    '3t25QaHxMGj&F^_w~}b=DsF<=<LaW^GYE%9*1+j8-e_LM9Sf42;o`;a*uxq-Of8V&W{Q!Mji9tJz`d4wQ0D1ikhL'
    '~cXEpI`lLn&qsFCwrg}uQlM66Iy#Z$qS{er%|1a)|KKYI@J-'
    'VMyk2kOylKG176GS?Ib0Lfm!Vkba{J{(OSK4ie%s&!QHKeAQnG-'
    '<JtGa^qpb9bO}DYmPk5#LPNHPv5PmF%FX=V>Uuxvyn_`xPmkUo4e0R!Gr^JXssqUnHc1cpqpec9*M4_dH7})ez~~'
    '_Z8k;xQ-'
    'r9$M7+6JD!&=$A@||k$!$GnzbB^9xqc6li39|_1eOFP;<;X+z@x~g+lv;KTdb`LHP2@oOxHQd+Uh6d$x#Z#C<Y1$'
    'EI)-c$_nY!D@XBI;e|O2OV7Bq=oL8n&{VE1E+_m@n@R~@%k#*j!-a6n>}fzpscmvgoz;eHsDFGX-'
    '2QPL9co9H`j!5mTJ;-=)+%_`3};D+tG*hsLirAC5raALZy-R9Oc~BQA+J&eU#=V?Unl{!j$Y=-'
    'ilF{nUXg9r<e|VB(A?s%8$IQGG1LEW!pxI>C9ddex#)wUmPNn)o#+$-B7Gl-'
    '*NxaW3*aw0oOY0hqv8E40yW`2QE*6|D_x_mGp<Yqz4MtB%s^8STt$d6r;Y>M_i{`(8{Zc8QrU5#%&LJaV6}W>;Pq'
    'r4et3_!qmeI(*_%(p{*g#>*}L#s4h^bjStPWpg3!y*j@uMb=7z^QiaHe3O@4g-uhCIKSeMgjQb>l+;fbOD}eM-'
    'z!>J5b}>K%J5&{BJ{S6M7y9rV`tT6?a5sJC<L*0@-_gqy<D4l<SlJLIrZiERb}~xo-'
    '6B{?tXxUixK~GsKlomnPre~XEDlQUnl%y}HBAhx2aBms7g@0;LLQg;N~c+NQdO-fG0oqj!SFkX-+K~s-'
    'tWREi!~_eG8Yc3a&f_QIQk}~qsOuCC^Cq{jF@Q5O>TmimUXf5cOV>6YM@)TH|%P<qy0c<1lPAmi;>p&UDE<dO-wO'
    '(yAhsb7@*Q%J)A7j!P7e0823&S4KHXw|E3x%G}UMyt3vr61^Yu4xL=?ig9QuOozgx6dF%-'
    'x^qRHwnjT%~DeO?CA;5ZOK0o^K8~X4F`tTz9aETu9DSd;|MOC2m`88I#S-G#GI@&>Tp50iP{n$s*-'
    'e9GKT7x(0Gx=G1PF&)OB-muJ>^wD2B95oYPqTKCd#0X9gqLjXVkSERe&gre=lGO$6@%=L;QPU?Sl)gamTAwx_2Mz'
    '`${U2VMk)C9yEBe0YJ;+%=Fpql5UtOK!tsDVF7@_7+g+70ZN4i;7&*a3!w!|xtgtuC9JBkH;AFWWZe7ra_FG+egz'
    'I4aH7(qpr-=p=G+?uUclQMqMtSi5o~>ZKg986ug5X5%lOy-'
    '4jN0sCW}(;2q}N=J2OPrLwSzbx%zS3_;Z*wYAo}nX`mhJJnb&)z^6|)YC3Sd~qSw2pa`{rMQmuDwCA7d@8R>1PY#'
    '9GVnq=JN7da*a?ro4?lV{1#p;?l*E=k4(Maj`7L2~Gmv&<Z>CoNZh!T#M3v1i{oOg>(Md$-o(z@-'
    'H+{yG_l^>Z*kq(A;d^}zJR1f1#;i$hJC!tG;yR7t6Y3k%p0ct1_j^2AC97aZT~fVA5-'
    '*w@n%=?P{?Dl~>ilp(XdKE9{w;?y;5^k&zNj?jcj0}UMPqQ<!;><GUUjLuZhSrq(SA~0zwh%yrR++aW02-'
    'wqWO!@$;+ps@10Cbu8dQ}GO=)-'
    'r|z3bA4@6&7gaAiFPE>yJga+O`5GnK`wJ1gFgnkrS|{gl}q?UkMv)JlEp*D}5QvM6VY#cbFLc@jHGh9zf8;pGIWU'
    '(i^-Jol02pR8raX_b8I_7-'
    '86Zy~T@DcWo)!V~w^P*0nKBGX*BP8tSh+jLaQ>y8=U+XIcFku|UhUiGRAhsr^yTTlc2ioCHp)g5VDoiQoP9u9}BV'
    'cypQRmYfO*c&6<T4#Xo8}tzPK?mc~v=QT|g$qA4@C_Qc>ZQhnekxSEreJSt1)H7-'
    'Jn{tn>ajmT(C!T9fL=3{UV|RMHF`~OJ<ba=pNR*1H+?vTJ`DQsD|*d79o|`eXDKflj#Sczrz*?-'
    'v{eF)!<3pWyp`WRX3C(cKjqlPM>1vaDd|07o2+U#U!0qcmf<-'
    'm^0Ins`RHFqjy&)X=b9$6Kd~Is;+~;q^(zQd9me)cTVS+fDPHH~V_?rQI97iUZhEF*{QJ)6ysZsR$27;IlMSKLtc'
    '`XEK)XU8G!ZX6dFqPt4o-+|Zih-mRtOqzjsry|m{8jY1B?vdRap-gvvm-zr;XARnsC{wfnA5yaQ&h}=WrD+ZB@`W'
    'P=Vc9!S#XkVjsbm56ruVfQ$5+_UurRU4f|Pz^YK@9A>^t^qP2f@2~Wluk;!#YGX%jCZ|nSZfg!vhG-'
    '=!o1>$Yq`+Whcw8mrn47LrG~k1@U363SO+F-j2NcS-xO|BlI#l}I?Izi$n#pII0MXZWlFY3-'
    '^19Y%Sj>EYk`rfPqAEdht@U`7xd6XbO~%(3BM|7(pP8%&>WoRindz}G7~T|_4Z@ICTniQ-YhuWis)*?7iDm3ho6Q'
    '|Dq^d14c3Gl+p&1r_F~<5ChM1P6kC|(9(ab~##v8QoaIhvc(=_lRm$~Ms3Q>+K9A{66wN)@}r(jff!8>|Q<P**ay'
    '=Ejk)TT^8D*<qeU=L^JTgc2eonBKyuh~Sec~7rtO|KcR&DZoCtL#0Mp=2-'
    'VpqOrKtfZ$_SF{FND~6t;+)8~e&QmUkddhAYI&P^vnw%?Em-<TYlz91^-'
    '%xsaRFi{gR?;R#k!>+=VEF4M^k$tvmVXhFPOd`bj5%=mF#$7|48z8NbPO-'
    ';4he~a?c8Y0KGOsd7wh88q#z{27oDuD;Pqv9ytJ%@S9k0&-'
    'p&S*CoSN8%@l!kj8Xlg0VaRY!>t;+NL;IpZ>_ZOyQ(JkSJQw?D>Y2!sBrPCf<GAwR{asoq1QB`*Z3L>Qf|_R={2k'
    'J*ay-8*LLiZ?3XX;HLvM4%jq?K%zQ)WHCoiBfL`-'
    'Oi#<FkOPM;er_$V~m6B~+N732rq0C)vq(n9OCgnNz#dX~YxxIR`WE`F=E1u-'
    'Y_3pjojBl(=oE|EJ2f9lwEn_iF_=&8UPvNlqGWsn(grf_$;MJU^NScukC1(uoj2MIqsVR8hwhR0%+hWjx=BORm2n'
    '9K{ad}Vx?z>h;%M34s$GhSA2`88p*+HzWQS`_hHio8%ooa+n9SpEJO%Ky(1T(_5VGAuJlxZOHof>tV)cBgBLjOzj'
    ';8+Fe^qPnCn#uJA0~LW8b4?mM6jvJsJ?SH{%z*5dKj<~9=r!l)HOJ^Ro#-{5^x>D>r*vwwDrtb?*rAJ(a5+-hJj-'
    '86FmhBrzS30I4tguwi>^tH`}?Kq*VS^U>NL4Abg*o!)Kzj_BE{mgpV)kJkcGFkrTW-'
    'UXkzmK8Kcj_;9@cCxldiDF2MP>lQGyo2h+#)hryX1*rC-CAIw|fn-'
    'U3)ePI~eGz1T_{1DQq8s^{d#Ljmvcs$V&TNl}4i<1??b<I&4ZGukE4KeYQK9b+-'
    'VoWTv@Cz+0+N_CcOEr+bQ;mLaR9MhJh3lIY><LgXlU|d}Tr;qSApSG3gdOV70$^emJ3=?+mln)8?3Y!mvtPQ=ewq'
    '1FnE9fa`P$QKno*mNn%uSaeUzRx+baWyHc+ZMR#m(*EtK)Kf6I&0Pb5t9oGfloBwM;JmLc25ONMubtkRE@n3fHsU'
    'tm=^z1C8?-xg{5;5Bm3-'
    'h|fd6YyzLg!8f*ulCJBaYioEjfbQE?sP0|nTVRV<4{&32DZJz(Ybp)1XT@2;VfT_DXapkI1ePvuY_xH4k(yrgPF|'
    '2WuazhGT#{cA`MXzrVrIHU7UHPjktU*SSM;CtBVE<a@6oYphAhA3iqcfSW2(SV6M5EB#5pg`27qxx`(-'
    'CCOu>zvt37kR)_dn03!q1y0Bw0^SRK6qv^v_={1|^HCfbVO{XNqUlpV1)(BCaZ*)~^e9~9OtotmQ)$Ypge#b;>=S'
    'ErlZMLNL93l1wJ>`^pOZm1aL`=@PiQPmaaWemj&wZalcgJN|UO9x1_qX7~&!x~eoPiz=W9h4dA*QJ){n`b?&$PwZ'
    'lonX<q!HS?*1<WKK=e9N9kcYkvEr#4hP8A?ogjO(o^B1BRu<Sc&=lFcyR)|&Ab-Ce>gwuZ@I-'
    'CYuB(Ne_L@j@)j)AG?$Znv=6_Mpudf0Z<{IPW><R2p{$|vM9cu4p_NGakjSOHK<HJz)%j$uEAu}KAES}P9-'
    'qLFl=rxz=H4fCqfZDuI?5Ny&6|M|>>Z=@RZKq_6R4e&DuVu;JD{}ZpiS+HbO5SgnA_d`t<cwh#sn(;ZtgYfF59T<'
    '?y*1i$X}~9xR(k;5Q)jWQNeLpiuE#0sh47g=1v5V8U{w7~9GKJ-ulIID%?GXU{b(fijcR}yUqf)s#~)GJKG-'
    'm#GJcMAg{QR>URJWht%X)79BvM$Jti2{)Ce2h4UpSV5A_%6pro2M+P=|5=j$4H^-PT*TQzKYt58O-'
    'F=MVd@Jt|61>YMA!Zif<E-'
    '>e=W^Wt|jO`8FY6To%zx4G3Hn3k#aHJP8^PQpBWYKFXQ=55zy+(sR)US;aak!q6x6@1Mr)#Dt-'
    'sN)r)MFVQaYo`+?~q5<izImZILVJrmpRSa%e$LlGUG-SDV%L7Nv}j68@)k)hg<kzUJ8w8Md-'
    'M4HA;KWMUekQ#JnF4?`au$=ADE#tK%`muq9pwN8oGC`Urgzj5V!mqNGDr>{2{&!OsOPUOQm5tu4H^St5WJY9%cbS'
    'Zy@KvjzINa*Dk>NC!QyXyL>%O$5x*fb|wN8ogBE1-'
    ')h*b4}G?1?$<NCJqyv^b<7y#(g;sbX^SG%Vq{nV&0ABoy2~b$bLDH{qiuqMu$FJM6X#xuNg(Jv7t7M+L$o9g;Ksh'
    'NVzh)l2TBrtBk$zQH~D1Eq8T}%BsW-qW^1_yqlCQRfqPF!RumVZ%T*+m$->fk&#sC{u5KKKE(pxD`+z6FnSekg-'
    '7W!q&}R1Wik$9jE2CcJQW_#x*~W}I~?c`#T{yljlb%^`gtIdC)U8q!`|37%N-'
    'ZrIV1a=Js#TIpwndw%z9~x4;_q=Yh{QHw)*hspvw%T4fBy&*cGP<?Kllg9j?ZPeJVV&SK(>Cf@F56mU{*DdkK!Z@'
    'eY3l%r9YHnZvv{m|3?o5Xye}i2ZUG`(=A(J|||r5_(M$z2*+R<^lI<3HRx_8tA;Vu~Jx6T`5bnQSKKgO8MHCvbpJ'
    'H*?hHF#xz<fQC}vDd(l9de6F(`^KU9;HT>jKiG$?m=}0s_V_)e5tjIivy7v3vv2+6rOcvqk$f-'
    'E=Y9u<x3_wtRFO(hagp=o6qu;(}$n4b++IK^7>W@FHF8iQ3#0x1QZYbXCgm-'
    ')GaKOtNj(^PI=WB{AzTtrr*b^4&q4QfE^c$>=E;Y1}Z?1`&jvDA4sm2p}jrTVN9_&zwY6U%42?E<Q6I%)n-'
    'RJIX=h>$+$MvHpvtQO{zud-t8OY4no|$htGv9b-'
    'zQgnyM{3iO+W1nN(bUHINF8O~Ax{OZO_Y;4Kg2xmk!<ojEn81*mw=`VrLoI6F_Co9Y|>sjmxakwo2pXu-BPYkP-'
    'N4&H;6C0g{KcnarxFRtl6^$-N($siH3RbR%fHdw!X-$pNzbc4lpo`Maz0k;p-'
    'lT8kcHel(!$|`c%WuBc9mx+y%wc9dUbyEj$}qVS=AIY7I9*A6p}o{L#lDPd!wcpo4*q+F0>O6ONZO@V!ip%JiBPd'
    'QIc23RbsM(EXKQalYVMxL~l3;P4fo@jBj16M)!Mz^E-'
    'SxFPRxX1>$Re8ZUeZ0N%Y^x^CDn*Q{f6WpgF?vo$2NqZBZJiFk;jH9DC_x&Ijvu;Va`w>aMxLyp>XGuoeZ0Xsjhv'
    'Ytq5$%!?u@82awGqa$)u$XT`=4Q5m8-'
    'DJIf52vwjul1attw>g@skeV`S~2&<;vNh*5X=U1*OjJ)_a%S`)1PS{J=;1z}9GFSanSKkwjyk!veK$#lTo-'
    '8NX2Xo-S!GsIpn#_Q3B=+B-'
    'IxkDG`?mFmJss)PzP28KI0cDdK_vkeyja1mULqYS}3d*hszK#%d2o@|W2b!K@rdtNA7zrFt24=V9&ahveXXcy5%('
    'sD=uaI8Tm_DqcHe2X5{?w)owRxrjZn*d;H9uP^-n~Rwn(;#7M_-ct{>8HW=L&hfZnBi^9w;}ybe7Z0n~KqHKe?Oc'
    'DCK8#B>4Mhw7l~Wcjui)<J$Yt=hjBVB`ikdxoJob9E}Up24c{)-'
    'dO#)Ge$^TgkNut{Lzh&aHlr5d<%eSQFWwSdSepY@gUO~vl8tQf5IBor&_>uyD9!uHHPhH1N8c%2d4;K=pNTbi@{n'
    'r9jl4HZ8b1`m>RX2Yjm7c*fCo{V;2Q&4+tF71hc#ZW8cv)4)RVa;B2s8j$prR%FLI|%;(O`m%+^UlU~!2na_<rOi'
    'P0~wJGL4edEk5r8Yi^9*U2hv63_XyVRQYPz=|ck{Zd|#VL59%+(tw_iCj}$<p@HWoiR)y;W7t+^~|#nJTIN{4GM9'
    '@37~d!Xmdl$S^I$(5D5^T{8(;9Y<i3Mt^LX)dO0;6R@diD|GA{iC;|`pw`V0boB7YCsQA+URD`1x4U9lloL!k+F{'
    '{6D=a!>j<+zy74`((-UjG6UJoJnb#OUV8|y1;p}CPJYCCA)A-'
    '!h$Y!$BlR9J)*tkB^dwN9Y!DClP=&?*BI_Dg;C%QE)Me(aaA%zWX@eB+t<Iy3X7GV}eS*SOP%JJM@L(rcPin<amJ'
    'nA&_VbW{S`Xe-Lz_u{(urX0O+Sd<m(C3(&)Ig*<#3m^2534>b7o>`${=H?+KolPXI=`S?Vc>(LZYv}j%DEd^}fxa'
    '0lVKaR;ZZ4aE)+>i$SwTAH4M@ZSyLh}@8G|mr!|}aReKdU@jOxiXF*&Cy7KV7jdY}u2`8vY2t1U*r3d8=GVNYul*'
    'nBfY<0twkH`Bw_p*nbKt__X5nus{30f#4QoVHWrc$x|=Z?az|DERSS@Oh44UvqYE1HrEwKso#6ZuU#fbmkZK%S>j'
    'zYRr6g%zRCm`S#IkPSI;>(}!8_@rnC%pZgR;ZCt3$yHA$NkU3E9e0?rw{#=x=*2Pj+R!TpeDRRhZkYq0DB9nSY%I'
    '1Uq;y2MrCYtKW=Q>{z9{w2LO)nz*&;iWrv>AO~F2UW-`6%2q2D9}BL#tUTb|-g5T$^^-=NyF-'
    'I~wD%UR@m33W9b~4SfIYjSCOmVHQ>ikDMH^G2I3?t61V=m>CQ<8e>uqLr4#OoLZraOj{icKd6Nvvo&#jx(1?{Yv!'
    '{<9f?q3#2$qoqoD6CfyG$f-'
    '*p5HMKFr}GIkX`g#9w37f_R#?>sY~J~Q7EX1>A9d@;;?+vqi2=rupNPZ`u^IJHrz%`=W;jTj@v^2IlK8valYHa{i'
    '%hHRIIDGO!v{BhFzdb&It8z+Nv8j7OjBi*fSq{;)eRBiVj4cFgA(Bab<b!0D=ty_l&xeG8ZVhZ|y$U(EInK*9N3+'
    '{tDF$=cFw~fs(XlO%R(yfh{kN~W)u8y=kFGSCD!*^?E9JR4W+IVa3g9X-4GQ|l)V+_5;yZe<MqC#|0v`-'
    't8`f0(0UbCHEb0AxdkB3$C3>5~?SD>z};LcG&)Bb`9KJ?>H?9@lul@{^&GxPoI&bvB_osyZa1v6iFX1)W=eC70-'
    'LVC?~dQC@a^O^g^b;VxJ45K8j@6=Kd`CdX+-4yeZ!;<u7y(IpbCB+dV<X2%&`8GCIEVtE`&LceKW~r(8-'
    '}sFMGhZUL<_)}>cN{nG?!-'
    'RF)v#(h7n|cILbquaI(hcR#cxS)UEBdJjAAjiWm9Oi2}1|F5U8^K@MnBAIJs8F`c|%R|K*604t6kJZ3X8!=CHV9f'
    ';R(<u)K`{RD<=9e_01>lC|;GT?=BYiAd&}JyB|0D^OvOs1TT~Akaub>K4vTH-'
    'U>YyZ3WIqZr7Y&2tU{rgQ>wBLDfG!Ww42Hq3miCfLHvx0+sanqE^%uNh5kPH<+1QkxuVGp!1g<y)SMmCYp?>02yq'
    'hp&|W%cjVDoxw7=psRdX*-Y-'
    '}2Z~z1l32~s7t5F55dEnP#t$wbaoQoi=dG~XunZHuW+Fd(9O~^If)9`T;Q7~XD7)7lMH8d(RVN(R8`Xnby<i+T>x'
    '=7Z_JmU&80+i;$L9{XS=km_&RXKZbu;7#nxMl&LnL0-hmE!#cK6r80uyaS-PXj_^BPEHhwAOD#_0Ykyt=1AubTpm'
    'Z-'
    'P||d2hEA%rO<Tyvy@$<@u%naeaXuaX>*6&ImJ~FEd}qip)2RK3qhvF{U=PsEq}+=}&D!sZA`k*|_?f)T#1N?&P16'
    '7Du;B9s5OM+<m+RJk5}ai{quwqee1uSan&q-'
    '&VSe(UhvLA2D{oeY{<H7FSo4pwjpa*wbkd>|Lf|>yc5o7BLWYclO2xy)N+fYYS_K7Kpmq2+y0;fnItbQhaLQ_cU+'
    'J%XUZf7iSpVw?|PO8$A4Df%HmdSUS@fuUZ?TQ=C3FEzm`{g$^1X&_X!9CTpe!>Tg%$?MD^nwNRn@0R=6>6$CsGoS'
    'G=8)j)8X{qo~QAa5<7|9EDb6u`O-'
    'vq}SC6Eoj+X1)o`e3j_K?dijp=r!Hw!&~V!8@Nxd)CSbXgW9+Sy_Y9NH^sI35jogxgVet^TigP3<j><?^8QFGIZ?'
    'Z=3^VtVUz5zGumEzY|7+|vy9J$DrI>$v7ar)Z#TM5B{0W`}C-)IB{MZlCD|#S2pd&slZiUnA-NUyxK+m?J7|u60b'
    '+He6dwZdA9XHgw>V(1f?Jz9C8b5<9px4(FW4;?9@r(gd9_pd7zb>rzYQv9SW6(+yE#oyXY?K<UN>w=Ir9$>n1u@k'
    'WEIlPSJD9zj{j!k#^4tj^YAJ9ehZ%*LZ%Yi&j+yT&GhaS4UkNi`5i?(X`mjQ6O6WCy)TRZsxkX({sY|&4dS{<Yv('
    '}d+`B<^^c3dU-tES2w?;(=-'
    'tDA&swUBx{g2eBTi&%{@l&mv9@Zt4Sbi8^68x|hH8J3lrYgZu2cs3&XOhCrcVQ8>59XnSf;?$6MY;|afW!VwXUr`'
    '^$Lu(;V@xvLLYA`wDi6jkI%s=FaJ+ExxH_{512bv@GxCyEx8KF|70o>B{kab=Mixagm&{GR%Of>Pck_L9is-eYTW'
    '~8Np{uqUxp}?j{P@F9A`unSR?3XtCfs}ddMneHLGvDOq90g{+YRr6@%zS$E;bi)72z|IFeb}Db?BPCr=FH5YHZ!R'
    'WBQ_=scqn6kostowc1Xv2i^O%(1R3VjPsUeGkmQsmQrO5>-'
    'p;X?TY1{j$Kf;7V;<u6=JQy0d_TOFZbDMeB?z&|#~$9@y=xA}+HtA)exNJdF1Ew|%~2Q>(F6(0>tf{kAS~?S3%?a'
    'ra4O3KLANR)bC&};{kB2iCQI}$F~bZu6a2nxh;FC#ks|Du89MNx*X(4j>BtVX;guSX-PDL2qQZtU1+{uA@a6+fS|'
    '$i+E11lFS^NNaxdS*komw#S88P!sj^J7Ueouj!FO8XRATwVT`fwP1IF(+LORs54Z5B|Q*3^b^3G1m%caK~0a`zDl'
    '3*RW+OXi5__>p26n<Dx{+KA5MdSdy>TWV!lie9lISySKQa_AkzuR8@_?Y-'
    '#cyAJK^FF>WBDaf`TiHCPGF=$9HoO#m;nxSpbv{iG&yEQ_MZMAXkc>oemRfk`YH>^Y45x&P6H+R_Mkc|zLG7C89n'
    'ZaS4G29~!u`*U4C3AE!g(3bby(VU!Cf3c-'
    'fY)v{#(h(vxRnZbkMe$Mt{|g~_tI2BP7}dO_RC*a0jCYX+KJ3BY4lxYzWt5Zs{#N|X1-'
    'C(e4pqwzvwk>>BGhJ8c>_t+$TnD82#-E)FzkO+<S9LLdTZKc<<HXwRxJ{n?F=KO-vM{yHPSNr<T0l?<T9J8jHThF'
    'En*{fl`C(*md|A=JnnQpS!DYIdCp=(<h?)f-'
    'D?f+!u$2C8Mif0%G>WqPu=1IvX`W;n@&$4e*C&10S5E8IJhriX#i1uxqa!(ra5|OC<{!b~MG$*GBky!~h%a>Y=T#'
    'E=KOrMwk9tINDYd9-'
    'TChIbIF7v&=L;DmbiC(4?k<+2;hcvjii81opppC!Z!O=((fmJ3WDivDAv0?{Q6LUS_^c%zW9*eAnnT!{{|e)Mg01'
    'CXL#BrY_g0i<+cRn`*mH%bl(}CB<TiEGfyA&IkHS?YkYtp;LrB^YN42>m6nDI$i0I_7&zQ9;4v*MU4J)5T9Rc!Sd'
    'zHkWhan!j6r@yx^gzF)R%|=5&WkZXCjz#~|x+I3o1wW6IB9#7(V<S|_UFL5?R@pKyWCL`Up8W{bt$tT4N?IU=^3V'
    '0N4lRyQ(0oj!V4eO3qUlC)vxr3EcZO$_tUK>N08xU5j&u$~Hq6BI17QDC-'
    'L@Ph|=R9Rs6o*m#2XJ!G=D2wM~=9|XMS6qj8H#6UQX1<s78h85eX?o2HdQBU8%`5JcA+=dYZI)4&JJcnL)G9bC4p'
    'Ez=nZrDJc72poE>9JW2JIyLW0;IVHSybGEf?Rb<?O5XP*1uCL;bUu8C`+_qc=cv<Rbh^nT8I2qw)LEK%DQD0=q+9'
    'VELsjly@!AWJ6=Dx2TJlfFNxCQUfnrRe@`O2aFe2LW|K3XnDg1?<ZR_V42~xfeB_DGennz`tbjyi@Y=)d@$8U{Rf'
    '(ezp8;DZ`IgSS&j1HD)f25&Yh~@gz)}eEr{wMsKI_&$bOl#8+bjFGcb_Zn3*q&nQs6y-'
    '!x{vnaq4K%zO{&H5=$P$LTfuxlaSB%}(x94{DQ6ZFW(Ynq=mLeNrcCjV!L7FPTQe#l3B^jNBh1L)L^!r6(R@f7w('
    'toBu&x?kkM%e-p+5rTD&o7wiJpB5-H{9<H5)+dD^K<GTL%kk%9C+MVFlw>73OZ-'
    'x%D8=|joZD?c#K&MA_ynW+^<d1HcF~S+gN7&=|duv!9v_Qw#rZ7!4#@PTvL^aUImdUz!sHcPHd$n+5rY7{~X<$yV'
    '8q3O6@M)()x04Dcv}DJ8A&{AZkY<81_RG^ZnQ1ojo}LVR&EVOX`DQWmHDKn`X6BpE%va3Jm&44ri(aEkZG5SXKD9'
    '}tHjD*eHN$M`Vk3WlhN(y*J1&(U33>9Z{Q!A?u(Q<Y5GkGd1&H2DXVEn>keMgHLnq=XmiD=VWdo03f0OO#sl5_uQ'
    ')k0Wm5bcy;b=TI1K0AB5Za*wpL<J0_i75$>0y`~7XqX6e(-'
    '%%4c7UUF}&Cnvyz=~V~ibEm{_Co4|8OOn!@m`5w;Z?;Qlc^taR3eZJ{=%CuyNPQWJ6QH85(l8b42~u*g@1ck2}N2'
    'vVTQi+00E!M9LBUG~e+^X!wWfu-'
    'Y^Ynl1pwq{;o=8I+K+rrG3!^~$+A8tn<zCf>uqt~pZ*DT{cSx}pI)Ww1fq%Qn)_?Ej_77QzpnzP18-'
    'p(}X5E&<X3>ryM>+13hcCsN^OQJh|f|lJw^cs5}i!Sbm(d$j9`FaUz?9Ru}l(ASPLr^}r4=i?egX(B|d|nidIrYM'
    '^dsaQXm>dkJnl<sQPgUFu^MqlZ3;MNj#E<E=co1xb_#kuGO)){g03&3(8bGIo9)h;(AhfwQPMc{lziYzCMgzXh)z'
    'Dt3!YFMO?u}PqZL1)xgqmas5~~SxKLY{Bcwa9DUW@<^C(&ywGM^taUkzryJ<NO?={0re!{6vN8|gJ3)TS=Aq3qa2'
    'U0zX>c7!nrtJkfS7KJlp`jadfTBWCCJ&cvjSL(>%NG~}UWFh^piTu!eix0YY&~@u6>}$FgXIHGl%9ji9=I0c6ejJ'
    'Gt#|Gf;;NI*2oe^HQEp|n>z`9C}F}1J`+Mf@^+(k7IhbplD>JFd8N;ucT0qs`W;C^>Y92{bX7Eg?^VTmF3FVIKwH'
    'C?#0(1FDVEu25BiG{@)IDShF4=XiprKzC4ufQTnK{MtW^$Paxc=|BCCi^M<gZ=W@93W{3`(+pA+=`zg8S4Y+F!QP'
    'C!|}{~f%M_Z)W(R~Y~(&Y;>?VrHlwKx>pVDyxJGT#uPl?-'
    '?<PrF<RBTbt*bnk+gvsj1<8?QSIHb>B&|D_qsOJ^*kg7bt{TTtaFM-x$ZGUf=3#NWJgl3PjREufVRwEH+-'
    '%knrKeh9yJ0h2H*Sber$UkA8~|<0>Ugop3lT+bNNntkq0Q{!wa*&yvn=rJpefciG{!xBL$tTk$J0bzEPSSoM{~8X'
    'bbuzN4b?!;C2B-'
    'HQej1y3ZA<atZJy>;5|Y9B*C!8f~LBH+1KeE^qTMV8l@l4*8zAR&hs(zeW2H%;(NPh^kIMca3Z~C2)!nh+DxK04J'
    '!KOQtIMF4)a5Ttrke_$Z?WdoFSGm9b|1n6WOxGS1wF-ke_8b;;ejyN6}+E4!(p@iH8t7Y%8+5FUN<Pv(W44cucB2'
    '49D`)v2|}Eo@|OoU}8($xf6k@W?`sVUJFKp{qSX0HN0$68Q*eT@w$o=?>Rg4xnqTQC(SWR-xRNx8DUVi0e&si!`?'
    '4C_&in{x5Kp1-AfabgEcTPS&g-'
    'uRminaq1OxrlUx;CI?O#9z+B@i=<uCc<RoCqT$9FJ6WEKnj$U(|UbC29<L^aHnE7tdYwpo&>eGkU(`#OEpH6X~YE'
    'c_aYU4(2Mo<?`vW41Q?>b8!ZXO})8uS+5^43zhetqe;uZoPiW+kPa)iQS1dsuI}hv$9H;^3nae66(+Udf9wqu+Ea'
    '4<Cb~KZB@ED#E^YMXpPGtTBv+WpNYCbgBmrzhErB<%>d>s`&iH1BEdzu(5MQ)%LcS2P^1lnj<g41TI=e$oZrXPgg'
    'xWoS*}54{aQH&Hi*x1LMD{v7(w9eMhL!^`*ituVAf8fwGp~+ljf)kvsJaIJ1{`%Us^y+^3H0P&-'
    '?2HkfM~GuLcqt{F<N*-fwMMz85gAO698>OyV0Q=3=R<rH=KN=^7F@s_$QA-'
    '4`qk+;=`NXJLr#pytlboB|5EoJWVtcIy{`~DjP#=pYCeK#@la48}O?#3&lLc|Q5kLD*Q<Me|ZOuLZ@tHr%Avu<bf'
    'S=9!yrOlDLsuA>p>)_F#K-'
    '6hg1D|(#!)2p8^eR_^nV|#Rd)vUw+Y;Xzo1yL&W3)*(M0u({x)$nUxT6jxozTMN`I<;$PguT7jUQiBm>Q!(%pnCG'
    'n<+?RPw-'
    '<;Xd1z;peIP+K5g4bU*$es;NAV0clR>xQ+Mu@D?3y(?$ciG(<JWG8tzj$_erSDZtl}Z&dd^OlTU5VQkVMV26dT2Q'
    'mRjo?A`rj=gf}Mcy5G@5A~PT&dxGqt-'
    'f?x@g05Ao<isO6&O}IiVW8sc=KZ=Hmsb39>EjQZ&DTt5A?;X1IZY@C;<&>wZfPwkyyO20d6IP!mHFDRSx<<yOkFv'
    'r?}yghBF$l1dq+KMnpRcES_izb3J1;xM_g0r+WC{uZ!}X+St`c3raIhERWGZ_hD*yA5y{9MTJ%c3Z{A~NH`&A$BW'
    'jE0d2=OcI#5+l4ZcYQQRl?gv0C!``Hs_`0>Bz-Cdb?cU|t&Xzr5^_h}RN=@)0_181f^wF#s)%`3ddnH;1pyy}p(d'
    '8Dk4PZj;Rc9N9aK)hS}$lU(65?)11(DM^I_I(H+_X|i{cL0|6Hsi$4rO<ma1ABIkgI(lMjM<%r_;1}&O*bCfZp0w'
    'IYXnM`)Q2*)7V=#EkP%i5Gd_DF)6W%2FCF2bXNRnLRxr&o$B9cONJut9e1rk^C+lJ4Asv{;YNNHe7Ved6pio}}f9'
    'k1meX<JOKP&jyU%@nO1$PSt&K(6`9RzQ9KeggMMa>3|3<mUgKUwg8s#BZUt~xN6_fs9-PuF-'
    'q&FB4;%ll~zXU3G;_*0wX)a5vJc~4FFDe;AxtRk$b*iw6#gzF^BsOl{xKex7+jj1fJubRv7Ga@^~-'
    'okFe9kd^K8jc~w=)7}1c9<<fQs^{HZaf;-Z3p4qu@t-'
    'w>55LH+u{EBDC`Yuf~s5V;^yffWRLbm>(VNyx7>r??}E!$9MIO(7N7Q6qTrAjx;U9&_fbQ3a(&!*ri<`c9US1@-'
    'KRtowYF%W-bpn+YN)ZMg9??Ig}s@Ddp;Gs<vs<ncdz07l*h|<HhcHOJZ83ZAfY|yqA~9(X5lmJ-F_8WcqOy&1I~;'
    'W_bHk8Q%&wuE@y_79`V$M)dQ}CmX7hH{ojB7yPGLf&v%yhw#}q$d7$_lcM;#VM)E$p9Ia};z{j1}G4kGVY<g6Lh{'
    'bE*8CrlUOD5qe1J<AN{s_O{6YZyVf~{T~3~$sNYpOMZa=bQNjRSG}U3FaS=#9Cp-'
    'O>M&GhC0`<DsJs#ysZTU2clH>Bg|~GQ_s(`p8bzh5lo09L(2(b}vl~O4h*332NLuqe5~O75o+}(5|XrJ@@H6@28v'
    '^0tfDs&I$Smjni-'
    '>J5w@H#=N(Vd9Mb0_Z#NDKlETX_U?Roa2P%KE@vj0Gn2}hDWfhIsLL;E(v7^Yn3+C=@}k@LG|5{MCsvk?CE;WZY3'
    'A-'
    'AO&;jTvdFI}_b<b(kC(8&%VDfqx($x=S0HWRYy{ZnV&3xM82Ks#V|A0U;Y|m8T^5VI`jHTq1~~pA1P?m;qhSvp-'
    '27A-vy9vjy}=39_SwO<hBeOWTcBBOQ|KQy!qce+7`sXj-9PFej)7ybKfQ*%JJ(SId68=Dnx#U-'
    '4+YH!D(I`jomwa8&EfCCychTc2w}I3V7I)-Zt2W!X~b@6!EQN|-O_>G@&`Tm7(I9^J-'
    'C=NLpwtc&dhSoOgy#eNo_Vz7Z1W%jj^OH2_#P&_7wlIt;EK)p1jqqBFWRO<i}{W+;IK?RkQnO9(4}$-'
    '|xftzMF7(<q~8Z$%p6av3Qg`1Qu`lAisWhXmpN)WrG+rz7vkSF7?sFsumjTtcfKLs$%m@Py9LOf(26?QLx(<%^F$'
    'ZQ;0b(<(Oc&tr094-'
    'p@PhAtg%(V@$O1<GLnFk80ri6*cPVsiE0b1(P!hBBK=~JrgwKKB4hHze5*$g?IBh`tAh4kq^5X4>i6ay_A`z2{X+'
    'OdT<~;coIGMC1<8RJvfguGlMfD)aGWzzttr1Q+}WpOaI<Yq!Hmd!oa+nT<Y6GY#P>*`|aK3h`Fg$4gG_DMz2va_!'
    'dmIoy67Id$1&S9n6m{fRo!)m^T=Om^uSdR^A&;#&<!~wYJ#(rUjx(8bhZ>UFf$6LVZ(T6em?dpHL4(&aH$Y{T#4m'
    'mkq9WvP5bhGkiW~jMkZkn474NoJG3O*42T-'
    'IxRFDp$XFg8rYYw#?<pF{HdnGmnGaMZw39iPqvwY7Cy{1pXirIdEYJsva^_N5`jBW^c-'
    'fI!}Q>~%rwL4!Sm?B!|1_2^x&VInQ5GvQ`DswXQl^r*-'
    'Bk1SIi985^YETVN^uO`wnt9J6vuI@RJwCPGY!TUqUs1z}Vs$?r*z_*)@(q6}1zeo2|w&Cbqokc_@C7jr&ghvCY0G'
    '8eZ>+7748}Xmc|xTGbFY8`VaYi2;zo)iLOu7o1<XAtlur6Fb<WXNff|vMk`T%oHn4jB)a!0s7q1!vIHJJfEwL>y5'
    'Nv?4gMkUhI}JYCJ1YVO=?QX^?^%nhNeS@5Qn|Ewg4<cm!zdq_0fpZ1m^3*b$~N(;T1&&tgXyM-LuP54NBO^ATVpX'
    'XYGd<^X5Lp4yaA7Y)+0VrE!T!8i>@Y4F+y#|>v6w2_X}!=%>MYBJZ`R<5?zl51t3(9G)*uE|9-'
    '*>e!rY_?)f^m6z#oQ3iF6X3pb7-pDfK&N368V7ek^S3P#JFqFL%nif9t|2J7;D`24t6|xw%9vQ-'
    'iqyJJ@NaI1J;$xke33cE-Za7H&PKT5Zvg9PJ)|zu!A?(Y?0Kw-'
    'c^5QL`<@!#jnp{TT?MVP3RXrb==12G{!|twsIC^cUj#O<W-g&YB=Jy<TG4mw(?{vSp`Oe-'
    '^x!ef!8bTFdpR@h=)uQ1GjFL&HfN^a-'
    '|x*;%nZj5IfS1Aj6`t!=T|;Z`Bt=4+7l|aQdxR!HkYq@iX1)n7FVL~qIah=sGLxOhz}bOnY0-3GpECA-WWXZH5mJ'
    'i`e4McZZKKV9-}8m<7fSFJfBt%DdU1M)u|@t##cqLl_y%Ya=}n5N2u%B;`}8`cpo-'
    'H8*39pZZ*Wax%%)qsEgHZI(SvAg?8gK5kEo$Idj!`b5VsQJ}PYI{j}aw!J<RdB}3rvE%18B=gffSTfpZ#l=l!b&8'
    '23%gXqDj^x&t=!P)fSZ=9JMoEZmtur2rLJZC1D+I*!h`Xu>3W+sVHT2NN}{dE}1iL1|>iq7i*S(fY~kqwPRb-'
    '5gmFT6nXgd3Q7^#qn<HwIo`i^UV>W9_%e2nra9$sq&Kz@axX_H~B2U0d|3*#ZqL8spfkI_S1D5SE!WF#3`=hHi65'
    'o^d5;f3(NoARCmHS>WbtQ}m28hPKKOry1TShU%hhqc-}q*1~dcO;oGS-rY)#;5jPv_^BW}Q$eAkpvrPVNgKg>Gj{'
    'Gf%z9hd5!n$&G1K&D$Ie6pJxUKA&4s8*51zt%dws?C=2mlNesgAiQ5Rl&7*^pk9u+h5Hv{wPN6m_v;i!Qhf);-'
    'o%Fe(VGW)u{_&(H;#_PYL>%B7Uzi=7VS{}jD{oCPwX(fK`n}an;6H(<>7VPZ%;hbL&^wjT&OY>Sm^(YbxzcxVnl2'
    'A-~>W^_}e6Y2R7do_e!-'
    'oq_xO~wL@indCTyBmaYg25WYlPYz43O7L55srsKoz2m$TymByQqPrM{1}|)QIe;!j)159+3*(-'
    '4%Qs!yH^&Q2Lua`80QAIiN-J%;N(;LJwX_4^E&5dou^`uwxG9y*--'
    'ub{uD>6K6(?+SIA|ey$m{=|NpMeppm7GmJ#=_>Y<4%kj;-'
    'u5=mWEr;7#Nyap_beZu11OME|E0^<F_;^2N=Wa%gYfI7X=L~GTKMq6kImEi@aBh=`1zqB?p;Ak@EQmnRo%-'
    '<GTni<Y{Lry_H8|h)MDA}Fl+AU-_Fc9(8eoM+rsfEXFu}<ChB#5A5AEx^C=S*^<C9u=G)WU9voxT?ed^2mDcM5>x'
    'A_YETojxx7DUlBr@9IHyZ~CVBb=PYyJY}9haMb64^E&5$MNEA%nN2djWm|`_E^r$H4Z=^J=jg-pYK^tpe}nWX68`'
    'E%tRB`bvRVa45h?h8K-m)mTynoWbiW+ac}zvW}{!DpT}*K<etLnX?t-'
    'laXo%KS%?ea)3ChvXjE=B2yVYquqUr8ZX9ffvYk<wozMiQPu4}su^?>d=8MLws-'
    'Q=*2bL68LTsV~$`;$eD!>x@)y*(^qA@D@7-EpMKHfIf#m~L$PYGIRTTK%ytN-(NHs0r}@bJBYRek8gKLza<2-'
    '2bjp9}=!*%8v2Y34B4ck9h8%U*eeMsdoQbHIChHSg_ZytlP^$WqSCQXbTcLlDoyYEv8AitjzWtN3?_zca&;13wfT'
    'HSj|~8Nrv~<C>asI@U=(jM10Zjz8c$?ip4lUBmpJ$B<dS2yY|SAW*XadCMlj;@1d-'
    'Ib<UKS5JIA+6ha&+MrSM=16dDgy{!rqu=uYe!l8>Wa^ET'
)


def load_prefix(suffix):
    from base64 import b85decode
    from hashlib import sha256
    from struct import unpack
    from zlib import decompress

    size = (2 * STATE_D + 1) * (2 * STATE_P + 1)
    payload = decompress(b85decode(PREFIX_TABLE + suffix))
    if len(payload) != 8 * size:
        raise ValueError("prefix table byte count")
    values = unpack(f"<{2 * size}I", payload)
    return (values[:size], values[size:]), sha256(payload).hexdigest()


def run(request):
    if request not in ({"phase": "terminal"}, {"phase": "energy"}):
        raise ValueError("expected one of the two declared phases")
    from pathlib import Path
    import json
    import runpy
    import sys

    directory = Path(__file__).resolve().parent
    if (directory / "files" / "moments.py").is_file():
        directory /= "files"
    gain = runpy.run_path(str(directory / "gain.py"))
    energy = runpy.run_path(str(directory / "energy.py"))
    threshold = energy["THRESHOLD"]
    grid_threshold = 40000 * threshold
    if threshold != gain["T"] or grid_threshold.denominator != 1:
        raise ValueError("threshold mismatch")
    grid_threshold = int(grid_threshold)
    ceiling = F(125977975227, 1000000000000)
    coefficient = 1 + 1 / (2 * ceiling)
    bound = json.loads((directory / "bound.json").read_text(encoding="utf-8"))
    for selector in (bound["projection"], bound["guarantee"]):
        c = selector["coefficient"]
        exponent = selector["exponent"]
        if F(int(c["numerator"]), int(c["denominator"])) != coefficient:
            raise ValueError("coefficient mismatch")
        if F(int(exponent["numerator"]), int(exponent["denominator"])) != 1:
            raise ValueError("exponent mismatch")
    result = {
        "phase": request["phase"],
        "threshold": str(threshold),
        "limiting_energy_ceiling": str(ceiling),
        "coefficient": str(coefficient),
    }
    prefix, prefix_hash = load_prefix(energy["PREFIX_TABLE_SUFFIX"])
    logs = PROBABILITY_LOG
    q = 1 << logs
    half = Gaussian(1, logs)
    dm, pm = STATE_D, STATE_P
    result["prefix_hash"] = prefix_hash
    result["prefix_layers"] = PREFIX_LAYERS
    result["suffix_layers"] = SUFFIX_LAYERS
    if request["phase"] == "energy":
        m, p = run_layers(q, half, T=grid_threshold, dx=800,
                          dw=INNOVATION_BIN, dep=PREFIX_LAYERS, dm=dm, pm=pm)
        if (tuple(m), tuple(p)) != prefix:
            raise ValueError("verified prefix table differs")
        print(f"Verified all {len(m)} prefix states.", file=sys.stderr, flush=True)
        moments = runpy.run_path(str(directory / "moments.py"))
        rr, uu = moments["compute"]()
        directed = energy["energy"](rr, uu, moments, gain)
        if F(directed["energy_upper_1e12"], 10**12) != ceiling:
            raise ValueError("energy ceiling mismatch")
        result.update(
            directed_energy=directed,
            directed_initial_majority=energy["initial"](rr, moments, gain),
            required_other_phase="terminal",
        )
        return result

    standard = Gaussian(2, logs)
    depth = PREFIX_LAYERS + SUFFIX_LAYERS
    m, p = run_layers(q, half, T=grid_threshold, dx=800, dw=INNOVATION_BIN,
                      dep=SUFFIX_LAYERS, dm=dm, pm=pm, initial=prefix)
    print(f"Recomputed {SUFFIX_LAYERS} terminal layers after the prefix.", file=sys.stderr, flush=True)
    truth_controls()
    literal_minimum_controls()
    upper = root(m, p, dm, pm, logs, standard, half, dx=800,
                 T=grid_threshold, dr=200, du=200, kd=2, kp=2)
    if 2 * upper >= q:
        raise ValueError("terminal majority is not below one half")
    result.update(
        terminal_upper_num=upper,
        terminal_upper_den=q,
        depth=depth,
        root_truth_controls=72,
        literal_root_minima=162,
        literal_NAE_controls=9,
        required_other_phase="energy",
    )
    return result
if __name__=='__main__':
    import json,sys
    raw=sys.stdin.read() if not sys.stdin.isatty() else ''
    print(json.dumps(run(json.loads(raw) if raw.strip() else {}),sort_keys=True))
