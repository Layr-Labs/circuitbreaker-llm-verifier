"""Full candidate experiment from rational covariance. No precomputed probability tables."""
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
def compute_plan(qshape,m,p,ns,ps,rest,Q,K):
    qp=[Q-y for y in p];qm=[Q-x for x in m];on=[Q-x for x in ns];z=[Q-x-y for x,y in zip(ns,ps)]
    M=[];P=[];ma=M.append;pa=P.append;ceil=(1<<K)-1
    for B,N,L,R,sm in qshape:
        u=0;v=(rest+sm)*Q*Q
        if B:
            f,g,w=B;zz=g(z)
            u=dot(map(mul,w,f(m)),zz)+dot(map(mul,w,f(qp)),g(ns))
            v+=dot(map(mul,w,f(p)),zz)+dot(map(mul,w,f(qm)),g(ps))
        if N:f,w=N;u+=Q*dot(w,f(qp))
        if L:f,w=L;v+=Q*dot(w,f(qm))
        if R:f,w=R;v+=Q*dot(w,f(on))
        ma(u>>K);pa((v+ceil)>>K)
    if any(x+y>Q or x<0 or y<0 for x,y in zip(M,P)):raise ArithmeticError('tail simplex')
    return M,P
def run_layers(Q,half_mass,T=-8800,dx=800,dw=1000,dep=10,dm=35,pm=60,*,debug=False,starting=None):
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

    bins=list(range(-nw,nw))
    for j in range(nj):
        pin=57*((j-pm)*dx-T)
        cut=max(0,min(N,nw+ce(pin,13*dw)))
        split=0<cut<=N and 13*bins[cut-1]*dw<pin<13*(bins[cut-1]+1)*dw
        if split:
            loc=cut-1
            lo,hi=bins[loc]*dw,(bins[loc]+1)*dw
            v,w=half_mass(pin,13*G)
            free=min(wg[loc],max(0,v-half_mass(lo,G)[1]))
            pinned=min(wg[loc]-free,max(0,half_mass(hi,G)[0]-w))
            weights=wg[:loc]+(pinned,)+wg[loc+1:]+(free,)
            omitted=wg[loc]-free-pinned
        else:
            weights=wg
            omitted=0
        m.extend([0]*nd)
        rest_aux=rest+sum(wg[:cut])
        p.extend([rest_aux-(pinned if split else 0)]*nd)
        for i in range(nd):
            y,r=divmod(41*(j-pm)+28*(i-dm),57)
            yp,rp=divmod(41*(i-dm)+28*(j-pm),57)
            ac(reducer([ix(i+u,j+u) for u in au],[ix(i+l,j+l) for l in al],aw,top,neg))
            sx=[ix(d,jp+y+pm) for d,jp in zip(ds[i],ps[r])]
            last=cut-int(split)
            px=[ix(d,jp+yp+pm) for d,jp in zip(dh[i][:last],ph[rp][:last])]+[neg]*(N-last)
            if split:
                sd=dm+ce(2337*(i-dm)*dx*13+1568*pin,3249*dx*13)
                sp=y+pm+ce(57*r*dx*13+1568*pin,3249*dx*13)
                sx.append(ix(sd,sp))
                px.append(ix(dh[i][loc],ph[rp][loc]+yp+pm))
            B,P,L,R,sm=reducer(sx,px,weights,top,neg)
            cc((B,P,L,R,sm+omitted))
    if debug:print('setup',time.process_time()-t,flush=True)
    if starting is not None:
        m,p=starting
        if len(m)!=size or len(p)!=size or any(not 0<=x<=Q or not 0<=y<=Q-x for x,y in zip(m,p)):
            raise ValueError("warm bank dimensions or categorical mass")
    for it in range(0 if starting is not None else 1,dep):
        mm,pp=compute_plan(Am,m,p,m,p,rest,Q,K)
        nm,nl=compute_plan(Cm,m,p,mm,pp,rest,Q,K)
        m=list(map(max,m,nm));p=list(map(min,p,nl))
        if any(x+y>Q or x<0 or y<0 for x,y in zip(m,p)):raise ArithmeticError('sup simplex')
        if debug:print('layer',it+1,m[pm*nd+dm],p[pm*nd+dm],time.process_time()-t,flush=True)
    return m,p

from functools import lru_cache
from itertools import combinations_with_replacement,product
def root(m,p,dm,pm,logs,standard,half,dx=800,T=-8800,dr=800,du=800,dv=800,kd=1,kp=1):
    # Q=message probability scale, dx and all innovation widths in 1/40000.
    Q=1<<logs;G=40000;nd=2*dm+1;np=2*pm+1;limit=nd*np
    if len(m)!=len(p) or len(m)!=limit:raise ValueError('root table extent')
    if not all(0<=x<=Q and 0<=y<=Q-x for x,y in zip(m,p)):raise ValueError('root table probabilities')
    if standard.p!=Q or half.p!=Q:raise ValueError('CDF scale')
    v_end=(3*G)//dv;u_end=(4*G)//du;r_end=(4*G)//dr
    def mass(c,a,b):return max(0,c(b,G)[0]-c(a,G)[1])
    # Every omitted strip mass is explicitly assigned root/line value +1.
    V=[2*mass(half,j*dv,(j+1)*dv) for j in range(v_end)]
    U=[2*mass(standard,j*du,(j+1)*du) for j in range(u_end)]
    root_cells=[]
    for j in range(-r_end,r_end):
        a,b=j*dr,(j+1)*dr
        if a<T<b:
            root_cells.extend(((T,mass(standard,a,T)),(b,mass(standard,T,b))))
        else:
            root_cells.append((b,mass(standard,a,b)))
    missV,missU,missR=Q-sum(V),Q-sum(U),Q-sum(w for _,w in root_cells)
    if min(missV,missU,missR)<0:raise ValueError('mass')
    step=114*dx
    # All V bins sharing both upward-rounded child-state indices can be combined exactly.
    @lru_cache(maxsize=4096)
    def offsets(d,p):
        b=[]
        for j,w in enumerate(V):
            if not w:continue
            e=56*dv*j;hi=e+56*dv
            k=(up(d+hi,step),up(p+hi,step),up(d-e,step),up(p-e,step))
            if b and b[-1][:4]==k:b[-1]=(*k,b[-1][4]+w)
            else:b.append((*k,w))
        return tuple(b)
    sh=logs;Q2=Q*Q
    # Replace C++'s 6-category line table by cumulative tails on exactly 4 payoff groups:
    # --,-0 | -+,00 | 0+ | ++ . Missing mass is put in ++.
    if kd not in (1,2,4) or kp not in (1,2,4):raise ValueError('root lattice factors')
    def line(R,U,positive):
        dd=step//kd;pd=step//kp
        x=up(28*(R+U),dd);y=up(97*R+41*U,pd)
        return shared(x,y,positive)
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
    whole=missR*Q
    for z,r_weight in root_cells:
        if not r_weight:continue
        b=z>T
        val=missU*Q
        for u,w in enumerate(U):
            if not w:continue
            if b:
                A,B,C=line(z,(u+1)*du,True);D,E,F=line(z,-u*du,True)
                pay=2*(C*E+B*F-B*E)+A*(Q-F)+(Q-C)*D
            else:pay=Q*(line(z,(u+1)*du,False)+line(z,-u*du,False))
            val+=w*up(pay,2*Q)
        whole+=r_weight*up(val,Q)
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
# Both fresh processes read the same bounded immutable intermediate.
BANK = (
    '{Wp48S^xk9=GL@E0stWa761SMbT8$j;D$R|8(jbuKfcfKoz?s|pLLQ=Vu=rV+(sQM;3Rn5|Gldc*x4|AEZf3nNXHNm+Dfrh0`RWn'
    'w>68<LjRu|43b<55;0HhV9zCz>jUC79WmB_wF6vHg+N_H35tLCUj|gC#lcoiUe`L?UP~rxDG#Tu??>-2il2}Cm|4;uT(pxcjki@w'
    'k_Sj(3_ue3yz!e{`+jT8XN&PDzxRaGb->*`on}p(br85vhOqfxa=R$aDn<;2ZrC|Z2|5;r!|%yk!R{}Z=;gq6QzAkmJhzZO?WBcD'
    'C_EM)91Ed_X*bO)H<Qy>EMulS!MP4#0(5QKohKvg$-o>8TQ!c-1a=$#BRi~k|B8Qd(V$9t5mW60EdwIKuesazvo!r>lfJYUvR?Il'
    'Pe#qK;L9CLmwjaqw|U9De2FO66*2e0ef+dA0vy?_h|m|&^wM*DG`LFW8J6%_(h!ZAIq&FNyXD#%S&6EgZS*$fRntYBPFi8shu59A'
    'PV=?b!^yHcS$narUR3{vz(mM^DwiW~Pa1d0Fn4lflst`i3f6e?aldQ`o8;1>p+*pjFLh|Su2Ue6=g<ZJ*f9;I;c1EJ?=gtbP))9r'
    '#dfni%?C@@FR$2MHxX}z_^X}a^WSoofY_;c_jM*HOLzs`>6AGJcuGVHJ@p4uJN{KeJoPah?sd8fi+@EI^j3FOUWoB`V$X(T7loax'
    '7ZP<fn`-9btA76lN68QV8*zoiiW0CJ&`8sS=i1q~u}y?ks2*H)PNHH=HTy;teVCcD4*6RhoyyPOip%UwJTwjfIgxZ<6;@@Fh*VO6'
    'jfsK)7Gj!RiR808ZAh`8-Peeo50HSP$B4;}fT_WBS%fRhUUO~GsKX)~<q3z9O!xado(}|EvWn=b72Ek&J^cN$L9-a^O6krJV+<VN'
    'DX8PU0g-gBD6WU@xb#gNUYxONtXY;!T41+Y<K^8dC|K7IS;@BQM!gs383t)t>=vv2Waa^u>>vFUD6}lX();cx#LaZ!9dHh2Ocn5Y'
    'yHn1SH4R6^Kf?&O>CBq8U55NZ=d`FLw_tOD_cjz|8o1AWicz(5s>Q1Ab`597+Cnk~x#xaP1<#z?kLZ;WT%vrVvZzYG3>B0EX2zx?'
    'BI@W*q7|S{@i>lFw0}y_ff`)YmHg*(`4$ftU=pM$iR+SY%ucnKO2(cKS)d3Go`d@eST!MIw7=N6Z5CBo@pmA~|KzR(BPsZ2NA?#%'
    'W|`m;;zBeOh-7ozfzE8t%{zhr3)r>!cNG`{pjD|$&_*%M7rT2O1OuQF8#@aD_K+rR8e(a^2S`UJt1FO%5G5R8GgEL|V}K=awK@uW'
    'XU`BCR-h||3ks`Ka!2O!ztMUW<UC1&y!g019%9$D+oJ;b^k?;%9^XrE=p0||qxveSPj$!^Tg!AU7;$`*Vn|gc1ry5~YTo4$R~pgD'
    '_8)tsRS6%bf$Y{H$Q>+4N7jAeq8q!^!tBAD+1U~9cKnJw3(|Z}50!<JuGwOEsy!C*N5c}==pS7+p4c*Y03n*VUdK#fR=#ey6gckZ'
    '7MuKi9$74;3G8ru%4#3ItP%!r!5duB3jD=n-C%WfiKyc@UTe1`1<kK~C#!Pv9W2}mxyN>1TM!I4Xc^==K{6hc>bzN}f8Gc#?gL9_'
    'CJ=<ZnNOLIbK9X9;Hfn->Av$7SU@3b&1E1#l|XUB3YCV1J^a)cah`j0)wHM*M3iX-L^aP9KIR@W%0tyF*qB;5MV3?IPgl#Vz;^=#'
    'NLFvv=f8tdg%3k!BwIW4kHMW?Z{Y8)D#vP7XP=)>GDC7^lW1O;GX2<of3Ua<CY!;0ZIkh^^xltNf6CgALq&lll!2ge7}+OiYaP$a'
    'Nt-pY<SGH$%V80dCxndJh8-gkq%DuiX*B3;>I3X*_CNy@j|c@=Z2&$^eFj`X3eHSO{l{=5vbM9)TK(tj5E0%Xe`#CEjWPG8haF_-'
    'c%K5i_F_MTU&h^_j;!O#5wuwf(zd1ek58!pkmzQZnts_pIebFyr~E$&dJuomX27*bueeaT$>H*~Gb)J-yr}5HOCN~n;iW_RT%9&g'
    'E%10q?$E$CSRWTl8W_5%Lb_l$fz1e?9T^2GRsl6W7mZxH<CcvK6AJ+cPGsHpb4bMSE`o?UrE|JyH^{5nHcRJ#g7e%PCw#?ySR;$l'
    '|C!-9^lP8!h&_Use=g>fplE0<fH}lmXq{HFq5CyM2_q6fl#XOnNhFCPP5^BL7&X&LR|B}HKz1F}`5EPF=rVpY%kNfD=zT`&#YrSI'
    'y|b0XM$cMEg7yyS-D*0R#thzH1=eLJ*|hU^p_I`%D!e9$@KE9^pxpJ!*QcDrLyOp29Ox>V{CDwxy*Z1RxX4BuU>##lrceKUku7gT'
    'Bp&)!#6LRZ4h~8h#=s8v2dBDEu4rgLz&UZo=d~8-rZ^?`3y#yGyPFH#dOTzdN)Z+NRF)7*Z32KSqBbiQNEcr7<2298^Jb#d_}xwy'
    '8qR62kxhYxwzvOOE6kWw2GmupwD&}{Sz>rso8IU&-uGER4Ny~L3O#5AmdrWY0d;SnrM|VSd^$VaPeJG0i_3V0dY7!ge#w>-cv=3m'
    'I6+&w7<>nV)g`;0P6t4VeFgqa6JM^MiU*Zus$Y@CNDDTG`wFsJ7FTTwenfy<K3*bb@9h{`32r&AX{X~>3)qZxAlbyFY(I(%tlp0g'
    'U+HE|VBCUgTGUg!{uaLj@PoEN;&DBx*}@G(a{@ngyw~v(Nq;nWLkit;h=7R^o-OWyhQ3z`@X{O2$A<)rL?mP>Cuzc>8urnB;O`VI'
    'G3|3<#RchO;1ZL0(u3vZRfp0}MMZMdDk5_YfV`Z|`^#mz+)=@aFs;(~%B`Nk4NG)%pEGWX_-vqZCn?CQ#JJ-z`+ox+z%tdpmdlSo'
    'W;OsO(gOEirm>(P&C7qbaqy;cp7~hX5@e@eeaa7DgH*V^iBLz<b&p}v3Z)C?*%|0vZIL#D6oOy}qG`#`4_D01?hVnX=^J9H*MhPv'
    '*Ypp-BEIwBA@>gYF@3|4pGWiik*L094{7)8nqhVQ6=ZB4saoVibR+q~`voN8rX0?y4Q<tHb3Ki6+7gYgY&@#1GjSb@%lS^VUCwQC'
    '0=>xQa=la=3m`4X%wnFdKDed+4O~ub-FL}c>c7>ez#Y6~U^2(SV-UkXXsivUdT)b7xO%<JXIL5xV7a*57WF2fKlpjRY|t(EXWygy'
    '$O3lmcwwDdrtlv8XX{J#r>}cgEaXH@^rig82Ws~qK`n_}o=OcVb)t*%<vKcjE7aOOrR@av`QLGbY(sByF}@ui5L=ElN8v}P)vt@Q'
    'b~Srs0J2DeuLssx`_A*5`OkB6D(ct;p3{?PC}z3*6kPfoq-8!v&tQch8hs)ta5Y>bkYU}fbwh0TKXNq*O<LqXWPU4txCxD&0uU%c'
    'DjQeF6#P8cU*qvr(m}*u>{JOj0`V&dKc-I|cT9XG#bGtT1_ern9`C~txidtk9=T>@P58YYwupcX6S<WQVXfc&-u%$;S}Cgi*n+WN'
    'b}Q2SVELcvN8@pY?1;v==Ut9U^`Q{j%ldW3rSneVbfwRL-uPO8-$-QuEN~Yy5TP<TQeV~{USXI0n%pcmTpD5vo;tSpbhA`zPIald'
    '?tOkz%rskiIWDa?oXxfcz=D%(<K1E8s7b{NCMp;ibt>EKKu4ho|CNnKa_n1@>-F{DLnNram)V#!-Ns1dVnnx+?tbmBMMQMyEqFZD'
    'V}akHrI{Sbm$bLQ(};=Pq6X==_NR)X9qRI-!d965@E0B~_dJ3yt7}WQSR|rPTf}Dc;}ZX%<8+D{qFVyA7L%9`7!%T^IwTFk%!)MA'
    'ft}jgk$=>Kd9H%MtF1R9#0)!|ytH5APuqad@fT2a?gK-ncqkRLCZ1tXh-ipUl}N#D-$nxS@#hI;wyY#c@@Vau7bvdyc|?89XmG3!'
    'cf{Ie@x6FLADnfB>L(yXD35!y7yVJ=B1CROh7%#lGI{NASSI>Fg2NEHV6Bs#T63QDD(@bp!m8rby%!uU{y2?1bt+5`m?V>)@hZ?<'
    '|0#WwrsTNMKd)#7F}0sdPH9%F-zv}GUH0@+7IE3N`Xxw^)!s-h_#0-I_en9>pk-_jmArYXK=vsM$tF^b)IQ8eF%SA%+zN-O!mX_R'
    '=04@>OmS~Eq_P6ee^+3F-H*Sn&*olx8;YpXqfBc$LoKV+L1~Cr8klVG1!w1e*=5e3bEL%o2wmV_i*n_><tnsjXIEvhzu7txwhRo9'
    'S_}Gu`crF!;X(ilRZfFU`|eDcoQZXTsct<LU~wK{s<@Q^uolNw)$;2D%-F5&5k(kUMgE8ER@ZKQs^muhxQ$l8BB*ijFcf?AdO-xh'
    'Bi#cN^;_!U$dzb=qA=h3LU_q`%?GD+rv3r|J7~iM@fFs89R@*DM258k;IIAE#3m(Xv~Mg<M+B-;0TNkAh;Es152>iAP7e{~uLpd|'
    'D!qnBB(`h!w3!E)J{YcDrgQdIuH=}88EC&?H$ibC_tdj5W3_Ee6b&;Oo3W~q4IDZch<I4S(TxJtHcyF&*)#^$akiu^*=Q3n<xfMZ'
    'Dg1tGYh?2!a=*L^Lyu5%2x2#gk0@fp8LX<<+{9~V$1M6Vnh6vU-o}Pi@VM3GX-C|DkXa|X9nV&$d3-*u51VT_qo59A=>Z`7&fYkU'
    '77^NbbCyxbZol)yWog{Pf~J=J3o7$G|L3=bxR%2`C%TgTp3Rj(2b1Wx^HIehNaDI=Y`W)y$iL3uN>^}w&4>pOhN&p8`WSg(eJr@n'
    '=AiH`BzeV0!6ZaiNH=ZIlrKWSxAy=`x)|k*B{;^)baQ?GV@-;k<7x}I<|i&xhr6PrvY9YdF^v2W34Lq*ka`%5zjOCHZKm&X-;Mr7'
    '4>h4ng!UVLZ4#{_6TFurqw@~?y3^yh-C$Psap%K^Bx#;lfAbkR_~fr~=LI5q8jr%YcN;UeNfL`mVp2dvU!5~owLYifA72iOs?M~M'
    'dJ|`;nY$#&@iA|*B;|3;&Uh?33VaiI=hDyg0T_+^q$KPw6g5XT@)h1)8lR4`2y)ywm((1IOp<8>jqVVPg->Z}820!`VaS&b2_5m~'
    'E5!N0<<*QFl7Z0<%Ais+c4%lUy0Kk|(oOFHnGvfSO4`WZTlCkS#dvMY`mmHeT>Hk=X9`#kR%<ysY3MH&KxTWGOtT!fzILcHbt}y('
    'd94MK5)LLzHkzEjJckyx?Lr&pN67z%-E9OVfz${EIdAGN0Sc3{-*#GVY=88u(0naue3ZZ*ZS!jb@gu?M0cP={`Y1@D-t^iPP_+kx'
    '?`pUN#b2%@EP>8B?e4)u@2vHT2^Xrxd1&+QxuVwnmN>o^%g}B*s<m9dBkMsT4|2;&=36w#6#!t39>f`jFWu%(P?mQVGxD!PVfY3P'
    '9W}%`yCca%*KM>1rz$^)%WqPT|5t4ZKm4N#zmGH5NSn=F+^{%hr<o--hTZ!6>}Mo>XCfXH__Sd7SkV&SD1a~=%ywY3+U_>2$lp(%'
    'mun@pJ<pCbL6N}}I~pu*r%{bKnq#etl&>@tPOx1fjv6Qx$w;UO9KOS*l<mbeKQU<iQ2KFc0A&1~nin%QPt>RMb0HjT|NZ=omZW0I'
    ')2(}aWcQ5*H@Fg$h0^EmBY$OK{-zDVb>K@MCXLUnmmqKpx<IV3Wu0j1JaDFCPwDkJ*%PZ?<l$<wV=D%&M*W{UNUKdB<U>yh3Bo%p'
    '%05`V-+lKKnQ@pYOHzd|>ZcGO<hZNbHLkP!AJ!4@oEL{lO<+ctZ7wAl|0U~!NhmIfpa)*SN^Bv6EjhN<kIMhODMgnGEa-(k+mqR*'
    'd6z4gvCkPmM=aX85J8P@9-x04ZiW!1HkF(gR!bnEfM0ue3Nc0`wM#UVZ>SByg?!M`1dw-_8tU~ds5hrq+ai+X^+ykYgy4jzN|xC;'
    'NIG~tVtOU8@Fk;W$-1eQ)%ZrjL70gnQ9qRVf1j=mO2=312Ly0JcGahlpPc?*Ela4B>9TG|oY83MVi>iwB*JQEw?hAEq(%exYgyPV'
    '9ob0gJZ7{HJ%MYIyA2d!98KVb`&YY`f!(lU8Tb@bLyI;>E)UJ;!?9s-sWA0|wT-Txhs8y})FB=Bkm*Losukfm8#XOV`{QKPB0~Dt'
    'd=}O{bihYbmVQve-ccz^mq70l#UtZBC9aW+k$Lp^jc71^Ha8C=ci%9aq&T_8htW@OQRm5!ykpv&p~l&NX^({cIMEa+4``nD+{P|Q'
    '$Um8TAe<z%UjE_akae-nf6VOlVHauRVPpE6>~P1M#&*gbmrG)U3ETHzaR(A5I-8XL3H-8Y%Kf?Q1yq%S*_vzdHNc_j=M-~Ef!|u>'
    '+2w<iUhN3`i&cG1mFVTjvNL_r81vE4w*Ch7gk5}b7vupLP1(~eH$);=R;rST5I2~-7cG`@w=dC;?jB=~X{MsXwsA@_xPd8iL@|H+'
    'Lr5y-ox8I(-Cm%|IH9E}5U5dW9JPMh<Ad@-7Aoowlll{yOn#5+JuOmq+x?392j>!d<gna)bvy+eelJi}&$FWLgsw`;R@7k4WS0x&'
    'h@}j%#6>Zb2@(L>(HTrDg?c|cPq*Sf?o;4jxmPA_duuUvK-qZ2n}m%WsY&bOj8b6Yc74#1GSx_H&5tR^G=a==litp!9Jkjh<~Ti9'
    '|LvF4`9=HH`}uyWmW%l36<YyZz9FC?7a_@0V7agmpdW?i>__wcCRu_Lfo&;rq@1RcIOd+(8D5v1t#N<8wRhmK6b`H33s&uHbfTFm'
    'a&1^W*<H8+PiREGzlg`04vTA6!YeHT2ZH;eAb7LSjKJbBp)qfMb~@@f>96!M^n|-vs1b`M(mFkPY_?jU<YVeyvhFftg!-5=bTeML'
    'PeOfgo$pd?&k{n3Jgf=2c9@952NULIx$l+?P>Qm(EEK3vF;l2p#_QQv0xKFU7e<R*lB7OWHUm8JSZ9Dqld>pEor%K<g<oaHalD?H'
    'fjikMf(_`7H#Z?I3{ah9T9_^+eI0P)U_VZ|_3%&J^%f<507&KMS+smnLoia)yu~*`VQKJP@Pp3tJ|^)CEvBZ^%rR|TF^i&ZsSJm|'
    'HEKI^AP><{=_ZIu8Z@5wxGy6Um?XllTyrGgb8MB?3eK|v=e6sY73|G9D|~&<G$-o=KVGRy3@WZIcUsq5<vUV9cS$Xi>F^9i;T6N_'
    '1C5~+Xzz>3SD9OzGk>@d6*(A4(97fQoIGKTkE5^e+d=4>59!cKTCT@UBVj*LWRcqway{D3C;uj9k40|tw{i^7RKy~Xyz>vnq}!dv'
    'IL}1xggUJzO8@4bZKc80P@D<`_|6bTXr#5YqXsc?hEE%4T40aqlCpL$v;&`>izDgI3R@s+mK-<!tR2@v$))G9eDEAGa)T@QLV@%^'
    'uSf$MdQI^p@WgD~4-!0f6G-%V*kXyOPM`!bfS@|R>R0}f$jG!QXa@}`wm9uabb{>U_)Y*LLlF0Ke%kYT74gwbK79=fXtID`%(jcn'
    '%cu&pB~>+Em9Yy_f1t!Xv_xoLil-#^)jB1rGyw)VG=MON#HORiKI(t>DT?GmJzi>wnk@Q)`p<1?qek;^t4J-quW6467{<9c@9fIZ'
    ')r39BhM+#)nMvCuN7YUJR+g*t(|~M`^ST$=ja>I6KC{w|ey7`UcjhTj!V6=y1d12VFIT7)O~dxH;<NEMH=4m9l72dGuDsr4Mxp};'
    '`6*OBI}41Ej-Fnf0qAU$WHB4@rW;dO9E5bU-RX64%M=QjB<yguV$%Q>n2V*7_G$<eMr~SwcHWTl?$TR|pqY0VIe6t}v=GX4w&^K2'
    '`%lKxP0D63+2#$$&;ja_GJ3PnwX8@$ivku2F|!<7MnFhnns~Uk<hdTv$`df=HSf9jR)*2i%35EOlJ~F@hqYpU>B^Q5>>L3_`T=Jf'
    'Y9Cj$O%1sh>$2kwn;uT5wIO@ZD|9xUDdo`!r8_|Ez@PqFtva{zr#3oNG-E^Ba1J6aCHs&4Nywe@99d}VbVBY0b?L}(FiG006_~sg'
    'W(LP_U+se3G&Ge#22t^B4p7F3_?j6hQ(#{|&Sk~(C{HMW+DV1~h^Jpf1zwd-Cf;^wSUdM6ZtitzIFF?hJ0q4=n3I|5zfL0y9UqXS'
    'jW}L=`b<Di^a2^3-U~VmVhC9wPlzwU)d%|1Bvu<skm7Vij0#3G!~+Azps=<^GaM<V{o*e^|D7Q(nmsgjC9Yvbdh28XDesglV*W7d'
    'BPDc#vPp^D()efPn*a)^C+`c3n<diOu50%jeWhb;r4>-dJ?7U^i#{+*v1h1N{4s%*0AdzvKy5mFu3tJ5xF<T=VIpZ<=VYcT%(yg-'
    'j!X^r;FSTDi564MbFWOZkyzQzBQu!<k%+Rd$c5W0m_ZA4Z5%sOw>^nqiB+dU&tga`E^K_4d5U-4pnR`=_aBO(SYCb%j<ohv!99@9'
    '&1~~=b0-(+uLu@!T27vO1cbof9U9M`Y#D$Oc_!jg2M}jk!}XX@POhXOcL3>P$9sRKNTyt8%2MNZMm&jB)91fUx^a#ugs+LS_UDtv'
    'sA5aiFP2sbRmh37nl;1jyqKP4Omm=7cY{xmQVuaKo7&+hD#%Q_1;DJVVDv!BQTIqV1Ef*LQtq1C#!m2t!gC@^lYmbk$+avggYPt9'
    'tY5!$6oei=8lhFahD+7=WRg=~{QpXx4!PdG9VtxOygaQ$;Y)7$K0_Tf+*>0l+>&3e3~xzZPeJdguPQ(8q<bKR&^Js3GJAtNC^($V'
    '?Fewkv}z@ycq$!>XbEGIMub1>MB6smP01$SVGo(>3xart`39QEvPAiX$lkdfs5+=&@6WC--!?c&5Z0ngBYM?l`uuNy3RoNbi@R-2'
    'Di~X~PL{m<Y`wq^2SDGAKd+I4s3gX#>2mX+U@&%p{@xvFObSBR_>~aH`fWxpYY0Bc`wJ4zSI63xyvYN;Lr|cTCw+u4pJf&uG);(C'
    '`BgC%wM;|tZ1&|eK4)(a^+_WolnNj9$&&b7+hpdbO+~x81^smA@L!<sC`2ol9rS!+&DHDSB|TIP2_~ono-D$@n25z4cR`D=exHK%'
    'fvrheTiJg!+)n%4W>VjoJgefGAc~sqbHYAGqS;M5aXcxQam>&^GQn8x$RTiPNa6hyjHRWTD)rr@<s{Jz0E)}5GL0>jmX7pCkqpx{'
    'kh^nvWad<918k-r;-_YaY-I)gF8)R`sCyvVj`1>18L#~d61@|?kgn1^yK?KO?ZRHHZ!VMat;CRJ`pSkvWmzihUI*aDq$8XqtfLL-'
    'k;VJApk+qqDQt9suD|i};an(|<g0?!y};NXDH0QJU_)apa0r`69`h~$DQg~5lsr!w0|TY_w|#_gCQ92bfe?o#g7s$Po5^EBsC}LY'
    'H#nTE1)gj55xEk@n|AW`=-{dLTNlki)>dK#`v_xBF6hT0tc@{iOS=@XEHOrK#<CqO`glYS4b-p;i0ALjlYFJ2;Bmka+7<`~41K9Y'
    'BMhahvl#@M!H&NZf;?M&MF9yp;2t-y95u1|#33<nE(-d{bbMD>w29YVtGavaStwzlBHzL&Dkdk5km$MTm)h{*vC)O=jf+RwG(#1&'
    'ot=s@fcXCL+GX2gGx)%seD2x8KeDr#4L+m_N0RvKn0tx<l>gue%#&uuhpc)c%daA#X;no@GqI>}-yPB`(XomY5_3}&+5Ph>AeO}?'
    'tSEj_&pSh%%2v2t_i7=(srvzMx3r9XSr?~pVO>xDKmFy7x&>f2?04fY!D{3HIFr`brJMV2W<3UPM&A@Cwn5=AnbQo3R2jv2R(4%E'
    'Mmh-k2z@D`o<FOkKqwW8vwY#U2y5Ysd_~YH$uD?$-Q5h1B7ATD@x8>FkMEmqPx4RXJ%TD*8W<T$5sx0onu~I!1Tyyd7RAO@+Qx*R'
    'bW5gLF;gASvsUY?t9gsmEQ3e&F_VR~_FzMDvDt+b6y0eGB7^5Tc@25(zFEGQlBx_#CvMxy2jSWo0IcKXWiliTv-H|^wOC@99*_9z'
    'I@cFb$?Lyflc~FoyPtGKX<)f#<Roop&hbu->5m3&*)~-1omJj(YpK)kOs{}qKeJ7i95_p8ZeLG=G=1zqO%i$-c;>Jwl^J%Cp(~33'
    'IDA)KAvm-wHmFth4s-z`R>SvFc8srBwoK-cw5M#V@h=@YH@Sd=88VCFe=SNt_Jj0|hk>ql5m?!Q@;;irhsf*ho>KG*w|r=<OQ^8S'
    'SO!?Pv3BtH8DhM3Qq&NG_<t5gt@%7vPe=T8G9Fvydg?P*DGInsNHt35hMjWAwQT!vWuCq6ycK9op7LzFAx`fYIK$1GJNx-7x}V}Q'
    '^T_1t$3%j2Z$!)6WSz58v3x4shpOGFt}|Hy)){TeD`NC6umU`jypWPI>j8Z$fU)LT+$vT|#2)T%IoB}vj%nP)C#FmIW<w-45CZxW'
    'Wr!5Zd=4;$Paml#vSms!0TwASbcRQO^3M@f54-}tRow%4E)luA-+9efr>jL8e+rIbx0L$K_5&dh`3R@RBK>V3F+CL8JYS#oBV%n5'
    'F|?Xg8ZS}y*#r9c)<J)JHGmu&w&k2~+nJ`X4>e78$ytDxe^Ioh)WX^^5~FA&W^y$49NU`bz{|hXD(fs5p}`?FvgzDHJA`#!_G1I$'
    'PQWzU{}sf94S;?eB>WEoHTRiHlO^J_Jm7Voa;rS1Dm14g*WzWV*4G}??A1Rgok7VZ7*@L*{m#aeOzP3N%9P@qU_^Rn5Q8Iu9m`_)'
    'zS<a%bB=$u0uu_lw67(TN(NQ1g1`TO(f{ZyuOWrx>wKDZM~3+@d1TzLK2`P6s~JxIJTU$Y;eM57nuVcRFCDh+i}w+w<9v6s3o^XJ'
    'V62!mfrOIErk0oPN4C8KWK%OzK54kz&?acSD$G+*v;Nhmt1eO~0oPwgE4|gQ6C5@-=C(8WQSIL_3vq(f`;XA!(=`2S0B(u6ulSG1'
    '!fgzs=`}4s*Ue0>FJwH`U7#HRdAsqQq-_TdI{|0I4cztVACJ@0BEf~Hd<bu~J0g2_ct%?O`i708E{D{hBWP{GM{57tiNQ$P*4O!6'
    'SA?2t$)c;K56r6D(Y#?)XmdC&M8BewKQ;;@OS?-L(V>ThK}Y0@)wYN&PB+LSpYV2|!Ry19n_I|u8=U8l8ORTk-9PwtyBUx0-c*j<'
    '4Tp61m-F9x+}@=!8z1L8GIsQ`MX$+9#L3V&VT|67B<2X#;j@z_!z4B1D3XR)@!PF5)OW{(aw-c9_|gzG#q<7AgCIXGKcd`^=|s4Q'
    'uLY*Zv88j~cJxH@2K43EWm3!>{oCepO@O4ghDz{I#mSqT3s?GdxXG1B?@oMGk1kE$*pW-7BqtCB60%<gjDkj}Cn$j!wvpzrdQpl+'
    '5&49&gm^w*ZXp_M4FI{?c)o|J$i&Dl4-IH<cjd)+Xm4V?C5jN|s8za0H(q)h!A1t>Zbp_z>z3{AB#YZ#b<9!jmr?XAhrTRTSs2R('
    'H6;^LQNBqC$fvu(K>2Ob8trVcd5p{Lm~EJ>O3I13*2gcs@suvD@r1(#==Q8K_`O@INX4X;uh7y&uP|2Z!O2(l=*c>05M+;aamVAC'
    'c^lh|+L)4#2P2c}7qOnKlb$aGVY|UtNf)iycZ_qg8%NQP`2w?*U2dD1^iz{xhxtZay@qTnA<?o95!jUC;tA@mvbIfXp~PqQMwCq1'
    'ilV<)d>72VH(Ux(4ua;SE7P2z&Z|y^s=(-YYz~plUxcdkF@9NT-~|kLaenNIAkRmniMU>YHgnx*7yJ>cTKb+o51Np%#Mem{xgkqH'
    'BZA}?h4yBt@sCH&Ac6WBLXylRiYOZ%#qyAZ+?A^D4^fdQ3LwraBqDLAp_Bjm?uz*EKOaXgQ_|TMZbS9ZY+E{23M`nR8m|vN=AY$M'
    'zte)#+<;F~VJExZQV-3??<F=>PCza6MWeWn+7^?XVFrVfrPAB2B~#PYqy)^$90i1>W&zZ^((BCy9t<Lj1ezT)d_^C{bv2pru;#9@'
    'k3HFsEJwzGy(?%pI;fByj4<A7ugwh`@C(`R5C4fDZZV!R^8p3D_)d|CZa=#Xtaqg?pf2Lbi%Je6bd9W#jw?<v6gFHiWXS|8s=NcI'
    's+&=|;AQyV%n`;NwuhZ6nR+9MyH2a8c8@2|j~ZW`S(Vm$r{=hw;)FjjW^E8A{E$=78XKC4IQe6=SLM^GhP=Z(fMQHvx^?mruBB@m'
    's9(TQ(=B00mElVoT770Ml8r=%53JnVH>VM9?hR@r<u@SGyxlv?t7=liZ?g5scBDhvJC!k=BZsHw(dRKT4}d2&r}G8d0Bq1#gDR2S'
    '{2C7$1=Gug?TY_Wrz@NosKy@AQ2B947aNnGtxi13`UWVh1YDQ!8#Ze83&4d!el`)och>4wZ-xdkhFi;6;n@G~qU8w|ykwaulQ9IN'
    'kKB997_g~St+K?)f{0cj|DQzN#ww^cz-3cd`;{z?XB-XiZJ41&AYZvtk_`z!dP&wgbSH0}yACk{S6`f5pp?`yUjhq7O4<Hf)^*3w'
    'XL2&j;0xoTY<y<?U?XO)a~jss0o&0N!Oo~<afId|5rQ@-lTjNTnPca`ZJ~JkVlmV(Nn+Zy46q1>GH#rGX#Q}OY9yuEjIe8axH>j;'
    '7Np<8?I37Pru2BGE6R&b{gjGxyPU()yKEMkBxn81m}IZYDe>QzE_o<7tJ5HWq%t$81okSzqbuA9qk6GRts~`-Y2J(hdQo%fx#Jul'
    'H&nl?`_!TT*SNx}sqC+XlMk8d!(mg#Roe=R06ng-dUrhsV)02Bkc^Zif@qBbjB{n`@OW7Mt$chzYas~?_i?T3R*}#GR}WUccClm&'
    'Pv>^TvJEq$1>{jL6YyvuYhoXdr$z-^_Ftuv%IHKC1U9OEs+t2eqU9uoI5IK*`NJF#`{{n(qL8}5@cKIy#8L-1wt9VGM7ccAg%Wa7'
    '&8#oKw_{f$uJiC}VtZC6ZSwiy;E6W<o}$=w%wZ`qu0|}T$$cK*L=$h&2*U6Fs|lgBe~6_^Ef37X&Bfs5j`}um5SzU0{suQBcUV$c'
    '_Z9ZnAgnr2RD~^`bCV=w+%=2^T9*iQ_U}kPmVe*)c43+1_@UJeIx-s2^}~Uhm^<{@szN|p>Xg;`c=6SH)j7k?Yqi1_&@|3Smj6%y'
    'qj@q#j0xL5h`m=f!h}sJW3>_Mb9X}MyT;*WxSka;_dz#Gse<<6Zcd~{uK57N;+wOVAc(tVrln>A`uA}0d=VM8Ehs0o_6rr43e{hr'
    'R#5ohD-P39J5Ut0aj}~U9y-YM2$c*u*vOT2+x!oOQdbMo>X<wW{^3ZrTojre8daNeW`a&4Uu5sZ%5`>N{H~I{hR{&;!X31uhTFv3'
    '<EcOzsV0dEv}bl<n>8$}2NF@<QmiR=qsrKByUwGylQbkEXlhknuhvb27>o+w)c^4TK4Lw1C(oYbm&1TAz$iyvb7|gN)ByBuoYAyo'
    'dXy2@9Ldi!Y6jN2XVR|@su$^6Z`{A8p0@VC(CIa{1hb^U-I_N#evTN8_%!gDKs>O#n3%-20lPjed0Ht)eoTu|{$le{idgE!pD3<0'
    'ZGBIoFcP>;Q-ku2^&lGgj*{SE?R+SG<v;NlgS39Q8)hN=Z*PmEPMTo&Lc~RntwzyD*V=C4^GKg@SdxMH#e!}?f|Hin0&KErg&S)|'
    'NN8~s)X+m^jscNgqUS_7neFj}TZ#%AVL5U)%$U|qV2{whbTHL}H|6sHN3(m~%teGMkPCrKsTg$NBAc`7IXNQCiWy9fuTpCV?jvW$'
    '`O&Sv*#yZENn3`@Cy7x8LFi%G3nUs^NUs~?XWllwZ-u-UAl{SPIj)tQ{$4FJ<LBH)QxTL7!|#V?d`K-MlHC#H)Ywx%pa@kB$_Bf('
    'IL(=%FEt31)*W7i^3;6q@&2@P%+`?R*#VeL8_qG|kwDOjn|f58Nov%zZPJclr=liG!Usr&rMPDq1_e)?1$qhol*8C0mHNMAY)3?Z'
    '{NmN+<D!y`y#4i830aZ3&1h5dn~WtczlZZUnZwOrzQ5}=4jem)Xb}3V?8S(@AD8)GQ{<TBB+A=;#LcRG3Jc*vqAl+A_Bxp7DQWeb'
    '&VHD$wBx&9^@Mkrh~oM+1B~xS7<oNC{ll9fIDV$t!1aX-Pq%q?L12f3n7NZW3eoU%>ttB#S1(|1y=TI>283I60GreKxt?MaT)lhY'
    '0s$6bN@>DQ4;%lgjk6Y9lbY4y+6#fgEZ2Im4~lev=gMKFK&E`ztOGq;Fj<w9$uw$TqyX>V{D;<n`E|5A^?~M5(*l4(Br6l=m2nf`'
    'T|)i0?cc3zPyUlph;;W{UaL;OskyF^D5V$okJ(48Mor_1rDwn<?60T*tex?lfYIw5tVQX?tu9Z0*A?H!mm0dgl+viMa{7QN(7Vu-'
    '?T$J|C8DeQK5M<^b-FHU{RgNIz19C{OTxr*+OKRU)eIeY&%aN8)zNSo_K%kfxNa*@e=ADW(VDT%v-g%*)X`UD<b<59G435JFp=x;'
    'G{i=N-h1oJx^S;lQ~`bWQAJn1Ep3zS1*BF2!Wz*6n82WP^rk|fkqqWC{eDRT#a+l=w`G-zn?;hduuxOv01S&9b}1=cDzLgKwN>Yi'
    'V`oCD9y;+8Z4hTOH@x5l$1N1;tYPUMIG+et6(~SWsX_?1&|S?{l#cDzVEg815_5aBI-za$vdE|Q_`@4Uq8@sSI8R1wd<M^IGgb1k'
    '3P*x6l$zIcrm;S`%Y|=jd<%KyqB&D)46fKOtTqlpqxIC30U-JNaW{!q*4nXRrO9oLIv}0s$u-xrumx{lpt7%w<*v`)1SkV{%f`iE'
    'r-xTE^h!Tv1y#Btar{@!K-IxZuD>{4Z@HCma&3Jr7l#-FxWhPwME17-&1lJM*ijLqv1Uh(k#$DS+SB|bzu@U);t-P<%nBpTNf7Z6'
    '9I07`5*6Bq({(7ijde#PGd$#)bpIN^7FsO5C7@3ee#zKz;$tTcEm*$PEVq68=E4CnX@xik&M@{(0$`YyfG_uU?)^R-@bv(7L)~c-'
    'ZA3dAmRE8`*RkFpXVP{HzCwp7xT-a+j-#uoo?G=;^A>0%!3Ixw&m7O4uW!Nwbeel8#<eE>g#ooHXlfz1-&w;pIUGzp*dD5)@)*($'
    'BW0rL*#<iS?=M68A`#!b>C6~La*a-Fr$yzc2iyaSlfeSvJu}eF8@@m&U<H3I693zoo4n|F|Lib#WJNb8VRVQX_ri*x7%!?lsVrnM'
    'YZriEj%ue0;CZsFp}I$4jJkfzz%CCwEpB_Y|Ji`W^QP)<GGbeTqt&XJDzqiNEXwceV@zq%&6bMCo7)8&0iREOCNLbUhlpq)mbROM'
    'iq1I_N0HYR=_~m9Qc}&bE%#tGHB2l!@2DGwYKkFw;SZwwwRle0Y9~t##KiZYhk@v{%T(O*wadN)=tFOcU4t((<bnuLH@jp_%lx>z'
    '7y?(6_s4HC1gB;!2{$ijqY&>Yxk!#Z)=n}#ZK$NVy86lS!j^#xioq4QZ-K=57FyfT1!SwTX=c!6)WvY&n+6(cx5$z7RsXzAIX(!E'
    '^<cBl+!7L^P=NF!E`}uMC3ePjWSFOuo?O;^IjvnX3zRmA7ad?%lV=$M6`3STc@sonkD&0~?iz`cOu|3{!WIjY2KPA2poC#KZLKMt'
    'gB+sJ3Kq@4Yb9S-9>c~{-#jc6Jw6*18V_l9Kq8>^sh3&Av(;|RE?{VNAyom5js8phBHZBwb!Zp6729{Y6+L;^7LEE9ES+&LK2gPx'
    'HbmYcODKfDg{9)7H~&ED*E^aUp@0|EJp2-ZF<V|+tHcsh6J=3)HJxWLAnp8U8S6##VCDUj1LMWdIbvB}5cVVW#o4SbnVP4^@{jGQ'
    '@6s=3gI{}kXlHkFMV%hh7cHKg=r2reN?}BS04DwJz5uwfG$e&MwoB)Zo&1_$2C)6{u^<SMe+22IQEGHvj-qh_1SB8@xe>^$$?aLt'
    '?C|}jiUG3xa#UE`XLY`>9fy#|-QU{aU0pd5AFZ2>pwQ`UO-R%{4tzQ+lSG4WC0ea9u{6#7_GG-4s@O;Jk6;$9<J6dz7Dp>BN*quc'
    'dGo`K)Sk$E23rFjhAyKG>Rtw<ygH&2Dj<gfgdRA)kN@tI#F-mQNBnA~Xs8sVAPn;50bhFqzVG)`nyFEH!!PWY^^iMecY8X(abLm5'
    't{yWpcfpb)>@gLN97_awKa@Y7;85#h3jiLM0>22JF?#_f0fbowlxI;_n!(_=U>mN0_tmfcibpoF1u)>~3s7ykdn?Z>#eV_9r(g7G'
    'P<iV=uQkram_U-<myVJBR=DzQ%UMz{spFJh;Sl7he)xJ{?nTj4r=9_t4d32mdG4ecRicx&0)6Gzrv2kPOAS%C`P=!r!DpkP>!|y?'
    '=X*v0b13X+An@T^7-)Q4l2bJ_9W^gqOIKm{pSJ9D)G~D`?74Sp1`=~QAvb3t&BK?-!4(*vZN^JMTke2>oT$}z8FC&c$n%2ca8~vi'
    '!>eBSTkSso7HKZcE;k<)6Ev5p`HA9mE0GA;Y~VmEwb5Tg_(8TlIfk&9$f*z<QJAN7l6sOELkH@|vX!ephruQ7L|E8ZTcZRVV|X-p'
    '$P^*kaLl+n3I+1zA^tvD^6sH%0D^oTgEM>(W#ycXC$sOV#ww&Oi1?#ZNy(^rZ&<lTlW^ht#{Y<{zu#OSX7>%Z0iOn)1(=H&Yu1N;'
    'e>aXQmi#g1gO5(K1%gi31N@n_|8kxg|A#OG>*@xn>bx26R3#OwUHI_Rs;AD0b`1*mkYikLPf1zyUK>&=*DmUwO`|~CSd^JL3BUC2'
    'EBycHu8&htBwQsMn!V1YFOixik3iD@N386A$lb*Z42plj*}m<vdn(pOF;Y%zip9b&t4Rhb%4&b5B-6cfGGvKY5b=w>`NbXclA1jr'
    'S!kX#%CmQ^Q)YQ&ekI<(@o-xWmbOPHutxGYLY%QD@iRYne#3DmyOFMAg&%0kDOB=uMltZz0-mRB8l35hpo^MvE|hl6Wet6&U$^lU'
    '_<VT6*r60oLXcvEtapQwg$eH!jZk6swVKq7`6N$+28-dS-J4(M>Z3@2q3*?(+}woQG@{c$HmglN*m!fdsjJodjY0`j<w5HnBkFYw'
    '1N7M;C%g9@(l@41heP7RY>FYz<n#R6*6m_hYr<U__`Xq_11mpnY9Rv`=dBj9vkEEXl%=ff8*7kQz9QtRQAuGhQ4#+oICGGwlqbvS'
    'YTy~?{zAIw!OJ`iWQ-!TE%kNO8gVyb_v??AAhc(}A(c}G4Xssjd7U-5%~$^uY_F#oGGz|<wn7GH*(7Ykkm=({)MFp5(}K#-p4X@b'
    ')k6yj-DRgo^kk)2f5@VzN;9KXbNkd)E*U3>Q)mgOH}U^q&AgNBOhwWXZRP^IReq;(Rubkrjf$L2;0|WS+cZz!ZpS;+_0fKrhmAeM'
    'GM${gncqtaB0U~SMH*WkCBv!))^ijc%D`)WtU-+iJ{E7|u2RSR_on;(-L>~#zZULpRA{aNbv450Un&V5TRB4m4x*l#Lv=8;ZPve~'
    '36|DYL~rS3GdKG9oh$MyhTcvppRQ&+?<CaQtC4TqR!VU`hM1i_mH!E@knpGLE54&t&fQ!!Yt2}#Ig(UH<wdfi*Ckf&<gErpQcU#_'
    'f!oMLU6H@j;bcIkUY=$g-%l_BeDPUJsw5BXDbP1W@gf%d$OyPVDmo6hd0Tv^t2Yed9DZkp6S;IrVc@SGUz$G3S{<(mcElleNN~4('
    'FXi`*MWYH{!rhVPIvXsL5Ik0n_Q5PL8J!aa4EJ@&gaw<aESc1~;tO#;M;;IW=(@u`>V3wY6-_k?yv<qTYjqxK=39EdA^GW1yFwgN'
    '?p`s;Em!3|5Oi%Sgbu;ICssXT`A}yhrpH??RA9m%W+MJAfHoQ;*K8oQz#QCV^Wm|L7fC_t(GH0rUp`h<?V7<uH&R3DD%yS#CuChR'
    'F_Etmg;eIAerf`!5lC9Fr?@81SQs#f)4QGvIJ*cUBDQW-6z3HmO#IL61*KVi`5YTwXS;HADGbxN)*2WUUX0v0--E3^9B@+7vNtRF'
    'g0Lgdn$tFbK<3>ae1Z(=>@=GXrg#Q+hDz9&@%;;60sQ!I>K19!=eaiVn*}m?<%zkjNh;O%tv-K|;KeRZMHuYRLC4wXnfp{Lmd{Z&'
    'bD7S#*y~4L<8`Aa@UNB3Q-*l>M$-D+5-2ovDJD_^d=_U1H4-8$bYA%iUhBw)-CA@tPSew!swKT%w0N3k1HD>h@BwHnO>OP18spXl'
    '6Y5Cq=C;8&Uw*aMk=Kw)VJ>&?1FNV+zy#b)JQANKjsj!lxd&p~#$y7Jk!9|qF^f;_%1!=Av^J^u*Oy#AkN0mO3v10k(1tp`ccE)V'
    '6HN9juM>g?z)<RkS2>>vcLw;%9=K0K!=+sC1;yB;x3hbKaahuD#PVi$(t2FbIt<=CCxg+P#SJKpCq%)0mtUy~{(~OyN_0@Hm>PpV'
    'axkNSXsH<s*&xN(d{#D_*J3M4u)pdc8f>0KhBZoP1h;j`M_UbZ5vo@YOCc!w8##TgezqbAku)-TsO1MKe}Sr|;HNdK+xM?dMblu}'
    'CES9PJ!nP}d%0Y7A&ALY2M3-MZl@DoJ;nP#vmEW$A0IgN0HCsN5~LW9T#8e<HYdr^#Ozd3oNk&qVl|%Ud=YA$lmu229!;!1I9HE0'
    '2&n+mWzhV_vvs-eX|q!6Bys$CG32vQ$#tjdkdMg0x%pAxrLOm)=z5QU;JY_O*I3&`XY-(QH!5v&NW0ZCVfNLz^8f7R$1|@LFl}im'
    'CsfcO2f*wZzId$ys^B2&qt}(_{UU3}M<k0}O~;ajv>c=yABT85jKfh-!WhGOMGJ6Z`c6jHg=6?vn#9L_McR^n!}{3_2w9}{>zh|y'
    'U5Ru|&s7W&<4Slz{urlf)R^&>jGn>|IXTOnL`PF=b`n$AqGZ^SfD#A+NvSz}BoX%;VJ~)2>nxR!5ql9sy+C|p6WjK^(xWKX&mJU&'
    'Nqs8osGi{X3uxFRCTmQ|F$CGY3J2;eCJFTQJv7=6@I!Z3vRHn0hHN;aJBo?+vm4^O<+5evx{k&g51YfMueAG?E=@!qV{^9&mH@A3'
    'Y1eU`ek$eOvQ9`Vwjoy+ySEKB@wTWowdvfNgVgD*69sFN>25=|f*s$yxgfRQeDNa3G^jk~)~zIQOMQfop;=f`*h~kJkzQsr>k-<z'
    '6M#hisu%**wp2&(peE%9eRkO7aOIkvPo<ima1wJ_fWN7R$ZG>!m2lTTCmY|H(6&mSdgwG~3@ODgP;c4jx|%9W^)GlF43VrU53ra;'
    '1);zYwI;Z@!}CcQZqaO|(E1TR+e_=P>ZTK~Vbw@pj1iu*BgU5_?ADGf<wT93j!2~g39)G>*rj4QuX~S0;I{@}*tUqK9wFmix_csH'
    'GBy62+bu%rMi9tkD$IygU2!IeHU*_>mw%G$NXJ#`7hZ{w0CeyQ^HQjJ6lNj1#BOxRt|s-%Q(9ar)|!_c@Q?sD>v0(kCUN9Tja{E9'
    'DkEMN?(ZeVw2wa`U1o~!Bsafp5-uQk`1or?9k(e2h;kd%8jp9i_=M!Z4}DaFpu5gn(JN}6tXZW4xQ$!8j3D1A6vhH}>Bg3|VLS7>'
    '8M+mx^3zJcgY5`tr+Sx%efSfw2qxDqRn&ZYBQ+5xs;2r_UtJRD4@CyDjbs-hS2(^#_0OTpFM;{gvYVK-iIxTEB=M`?1Km1nt_Yem'
    'aOQ@pN;&xy48xC}`e8yr@z)Hp{x&DI3d00w*g~}zJ1_v1>cm-{{%!5%$Im%PD{}7!Qs>uEqhAIb%7kq2aP$;QqD139hW9HQEooQT'
    '@sX{a6>6H19WhNlnlpNm05@=|vY9??o`9o-=`HFItgZ#6L<rtjieUR6Jy?WF?iztlsACP=&K6M`qvi1C*MOXBkPjxdnH|cSlCOuK'
    'LBpuV$1l#@25d(!nONYhZ24(~)tATp&apD0h6hezxl+8qc~4KOe`ll>&J@>Fv)fr{VaV#egmsnA|CO$D`ZuZcNL8S<)jNZ=KEM-W'
    'V8l$P;u?l{fy@^-tB%fF)_L%tGLP5^6JP#5VFVT8{jX$koRbGCWErja=#RTy^c_4xmqzr~dt9NIh;wJO0{T4&)j7??W*aXbsA!(A'
    'oSxudI88YK{&QM+dqk%4R$sv6t;Wxa*z4<esiq?v_g2!qUN$S(O(T_D^+=9NjzYUWdmL`V%={}@rxDC>y?sFI=W+xMv3Gl{>(m+k'
    '0}gwP5veXvU5yGlM*~xHcYYXC(5sp8!%oO95~>NwHj`V8a`PJ>4pd+s>fn1dOeyQb;j|;1%3!$?lldcpb*^?BoC=QWde)+4G9>yI'
    'x(M+JyKpJV%x$IkQin%IU<f>`GWW>L!M)Ci?~(Zp7f~fxp#1(4hpESDPY2M^f(h|XK~elk1-_ZPfdc+@SJ$6C$Xx0nC&|)5omGN~'
    '48no5AhTNv<X?4~G1NgW=)^cIgf!9Ht0y|U2(Nh1Hl@2AL($;aHd^GL2AdcndC1Q7X*&v#b{pLTEF2ow7I5XPt@7d6)Jq0WD9k2`'
    '#?8c}RQ=EGPT;|&wHtuoFna_nehSs&OZg<Or#2^COwLpW3?$3wuAV(;^erw-y%vjyA*R8dp&7O9jM^0@eh;fM%u+vs_;3s&EJ`o!'
    'Tgkd`SrojW?(F4;7&~_r!1m>T0*SJ;XW8&SMiO)K%*1aVJu&XR=I42<a7UK?V&l9LlIZHXNE3<=dZ#XEtd2hLN(mJf#QJuU7VR8A'
    'U2a)zrk4g8)z<ar;(3{|%YJN<kdy+gU4nw8>EVqK@s@?HE@*F2zL_28j?%j@@2Q0=?fw(NsHoSEWF3tirbi728@ms>W(@C=C0Qzu'
    'sXw&a3?_JO|Njhuq6oKt+KUVbCiDe-4p$`rQ-?LdSz}ce>gdAsCl((F-*g7|t54DhNN8$iV@J+b=C9m)P+jCa-}`TP$^mw*CM+Ih'
    'DGmMBF!M>sI<o56XaiQPAQ7ix(s!&zq&o0U>5o08SysJ#jSb@#Y}1EfImgf8+~fRD@%-4XaDZgxs>Tb6@e?!m@6k?yMlQ;RTDL!_'
    'anXVWf=5Fq%y<{j34tNiZ;Y4;hRoHyF-~1ZC?ul-V!YULNOZxho;l?)Ouism)*N=Ch;_O`Nf;Y3;{Y&Fo}0r}n&Z1lv`w5nL5!|b'
    '9F117<31jV9s?%+Z8-Uhfb$6st?wq}MnmU02bz~wfeR7qtM(kOe?RVSs!I)ylVD7@V_%HT6-GbA?ktPQp*C+{Dq+8)4>c+dUeGb6'
    '(}88!lw+JjIZGR^OO=ZeCPR%^`kqTbT@Daj7wO4>wtDqIXS)3V-N>iUl;$t1M^e+L7ZUl8&KQ~WYIjJ7z<?4C_1g3ZB-ee3;wmg-'
    '-I(knWS-33`nQvFE4pB&rtt+7lzF#6p%-k8Efu((Tq~XEpA8weh;N`5hyOGOPktpmD0<e6@o_;?y4#HfIdtXFsEtzRPI_}+peO<L'
    'Dt6<f-wLXkJ44E8Ql()USwG^y^{=V?@)T5|g+p`%mB|`ZG9Fve-i&u9L^dPvv!Vbwv;3ZLZ|3V&<n1KsV>)oQV$_A9L5rypz_MCq'
    'S{l?!K3a7q1<J5+XiaMV40cP|+<F4DGSSi+l$yh@vI$@f_Q3xbjU#JCrW=sn_Y|lr<6NC};_qlE%zMAjYp9!npqcz;-tm@(4Hq8_'
    'kl!}gP;0zZQ|K0nReUPJ%*Wrt$OzaSW(EUrTeYx^AkVJ}*)<NoPbPN2;Ap1?nSNPUP_n#dlHX%))8R7EJ(fd5cTG1nds)lbx5#*7'
    '6S4nV77ID^euPQ}#K_4d0{4b6r*wpec2VH`hOaUaM285}hP}G~T!-J9+VV8o<#c4PL750ubsk~`_NP4O0IqbY>Kz>vGBNhLX;KCT'
    '^b_^`{cTsn1A)E?7>Zi&2`jA)x(^)u81>#~T>Kt3gOUT{Z|vv<YI-Q2=7G?V^IF9CfzJej(?76Ta<;}T6Pb)_)!TEb0fr0O7bA=T'
    'xF+@|+K+-6)x^1=RQlDnJCKemLnLQ3TUaiin4gAOo+vk)>AA<iMzq0d+}moC<OJCYh<`yCcj`uQ$#Q+haxX;3QZl~mVO-n>ILvlr'
    '4KS|$rC1D^xCKfy&GM9+JnqhHbhx>^vl%m6FfcnD%|wvVvM&_Z*QhO-`~TT#@7o|=p^6B0!`Y7@1G?C4$OF`zP4`Mx!6N{*%~Pu3'
    'T0H}uaO>CaI3#Xlc#6ssR_^!vS)eIpu#$c_@$6EyZ^4$5z3yJMc_z#sUhPRg_#uuKnEum5|6UsbLh)=S`&DqXc5EHcgYl(hkgF8U'
    '{ca9T7w$AZ(T0B;NZpv?%_>vkDiJ1*vDP4J;=-17n`<osambPob6Q>$46K^ywX*@hAkYy6@P%a-Xe^4IL5>Zb-b@Zwf;)=BJ*5u3'
    '2!(sk4wNP9w}}%By26K7;)d@bxLRDh!dv9=nzL41oIRzp^wWz$oFA>^uH<a;9_%Dw0{XVByz<JFWZ}h@2qSF4f58Pa-2AMocU9^K'
    'f&QJvP?SVc9_v=A2i`tISB+HBqbNwmuh$skggXCXMKV}8ZwZA2m+ZqywZ9}<i?bd_gC7V{)+lND-61Ef1NUmoKp16`FUjY`uF|7~'
    'a20fDdVI#z60qwG%BqpLN@vBPlpt$37Sg~!{HF>%y-aWBVlV$r4}H$EYO+PcPMZ=-aZmb}?!vzz#K0{J^;7&(W*T}$(qP&H*>0_u'
    'Q9fT*CYqzp5*Z%3hbyX<&=GiNq&#Ae8L#!axU*o#Ht_@~*J+2H$SqzB+uVi^!dzc-n<)r7+hvrwu~PexrMPpP=jjBB)vz>xr7bJR'
    's-v>Te~aW*j&U`Tjo>yYo}|j^N8Nx6XN;RIc~a*|r+!x!du|=iM+gA&q8zokn-@m@eC|7CW&%TSBFl3D896v)TjrtT2KykMCcmHg'
    '<8Qlzh=l1oB%qT5^36wNcZsi~j<5vNK7&EE0JTL6YF&duqFNMHM#)8l7KguxF%#^RcT$khZLM8$fh`#5<^ZaMMaa;Zv5zqd&3;a;'
    'c0$9MyePp2ziA|Q;Pt6faU3}DMNV^Kg3XXx+y-{AkXXT8ews#i%Rk#`-|N(&IbIn54cKgDQYi8oQ?5k0^+{+kY`<OExC}3{Q$I`*'
    '7uat=!bs}A#rc#&*K3mqd_G4I*b6v33W7|0v5sqr6d8T_Utuy9hFUYqd5SG)o_ronMGPBvUE|IzD7Z}^fD0D#%|=Qt0rHRWVA*6R'
    '|07~6=v;=g-7sO2-Y_{E079DpMO|i4HpQCjc#E=kEn9caw~DlG%6Ra^C2pb%AmW}2rnXzpNi<Z6*s(cZgf8f<!-hhpDC&ECcP!2f'
    'FNcIbZ5Dpg3SuxXHf7+LjZ38+oW!H(@R0zreQ>y3U)cxj?m{&Sk;@V@XUhBfyVI;NDts3Xxn{Zf`xI-Kk+W2oAo}y=$vOcdcb8>f'
    'eZDHOCFbt2IylhwTXc*J7+P(Qj#g{Y>~aaIrQ71f9Tt)eaV+##s|}bK{wjfZy2N6(i=`~oRize5%#L_V*Y6MW<q3l!FX^98@+Sw4'
    'DDnBPeGCyo(Fk6O{#~(GeL<UQ<2Rn3pHo3A=RDH49~RDADw(y8tyQWA{^)(?*AQEJ^4Y3c?8!IUM$}Cns~Y>HOLsb!Viu2<%7MJ='
    'n|x7+yS_<$Zvmur;cg+~2HIHV#F8NSIRlaV{+a(|2EFl?+9CQH>y#Z^qo=hm3+hp#qGX3HqY%nk#}Jk1H)ZdkmY63PB~CNL=T8MT'
    'B`5}m8b0Fd3WV=-X|~`=SJ5dMO*H(S&WKq%kE}kD%%Vb!6Djl%N*1@qg>I=ApTbw=q_SF_1h8r&*;qGug8b`uo?06!+?m<DUUY~K'
    'TIRh~^KwF@uE08T-lR2e*P0+0_4bi<=M)n5%*;3{W65T#1L^TSQ+)GIKO3ZUo>no8{rT+JAr=KQ$O+XctkIorS+!f~sE9P)uFHFy'
    'P0bPOzb+)IJwE$pQMoG6&@gxGYPxPn4&2V35`z%Y8_)s#bCzK&h_AFyq;F>87P`$og8eP<1rDLjDtwEt$JbX{R1J?sG^b#>_#5~O'
    '@*mDQEwY#5*=I#^5UWU7LAAtP(tynLZ81xE3)_RRM&rWXt<83E#P@@T3`?02_lhQ4$L!TpD1N#R)Uc5|Fcp|voPkXr_jtLjos#Ya'
    'Oq9(8^atF)BqV*;n|T?2qT)Kv(knegBO_w#2B;6wY;G<!og7A;(yRF%Oo!buAWWpNv^&3@PmOFVVD6jiKg;@Q48F!#(uF?O+C=-H'
    '^O@Kb&6X%=d!1sdjpohMO){^kG)7_m-u2s8u=ZkIi>_!5fYB#5V87q*PB<Cn8a0Hxg6o4i1rM#3II1j1SNl`Ix#n<%Q<|2O{_zpY'
    '4d9{Ymrx4A7KlhsLy;5hYI^$yH&k<ZkB`rJ;cCeIe2KBj1;E8FkXaRmS13qweVj}l3fv!KYLSrROiaru?70%5a@8UgRzB4*gq-|@'
    '<CN<9ZZFgeBw+wZ(^Z>vn}rs&kxs*A4m5zf$)2h)E|BH|E^k!sa1ueDfXDEw?-w)1s#G7*h+^_O2Z9;|5O92b57~u`6g+^37oXly'
    'LB{XN;qfv?sT1YT-<4-GIn*pmQY5W&+A4Q>$-#qw7PLD#U<CU@Ab=Vu0PcMvl2?+zQK}UEOR%W*^%naUu26YA<AU>woG{(p<sU=y'
    'R-+4Cz8c@J9!1skA7iPr)4I{W8R9OiTb5zGX;xj99+0vWU}5qDIe_1t9pmeiy4sI@yHembF-(peZHPg>9qwYAQ!d(?kM%F{9{xLH'
    '`Q)2SHT>Q$R1p6@f^svUdvE?RvLx)G!{LMFZIfTUq3kkYw1)%&5)4d5F23ME^(WB44{J;olh4B9r=w*%-6p9K8*G<4Pawy0HQwFt'
    'oWf4*QDTcJt&p2{2E!LU+)IM>v(u_CiCBNAOQsNDT?6DZ<D!g?*ev+9;?`UznCld~i+S2_hq8TDvaO0irbkJPo!brF^re#3`2F0^'
    '2a)8g4n&EottE-6MHT1-&tx)aOmv8epC&*vJ|C)yq0CcMLi+XsVW%S$oGZ{92IJDYK}%dmx2r62uZJtIHTzD%`UR7MS6fSegI4M6'
    'zl5b^FtEyThGVeRjplXIxKabw-IM6p+>y?l7I*;LWOTV`yg<m#k0Mbu)*!IyQyxc(NfgijWgZW)GSbAWTE%pNH7OKx7>_Sej8kuD'
    '@vsE*^kYVY;PbjjzI=PMDsSVj%`D7vOnEJX839fDO3!~}sJm)dQ`un237o)=ZfYP1j<Di|kvG$z|JUWMvq1e~EBtvcmuJfiXq{`*'
    '9`WCNj4?0O!xfSV4LqrBA&9HXjQBqml+)7l%mxbKuUt71f^v#i+~I48nV8;NQ<T>xReSZqI~*_u$^5_loREfn6%Va(;_%zq_saQf'
    '4o2iP+%0!h!-{RQ<A6^{SjQ4YcFagViQX5n*4JAk9rF0D9O_V^G!#dX(0~Q!Qhi`poE^+Hh%oCh5Z-a=lp`qPz&ZK_nC<sCYtk5k'
    'xluUdR0L!UtgTR~DB{D_CCsvi%3t~l#p$B7LE8qi8XJM3;r*EJCgI!C&CP_^orwecko<`qLx32R%IG^2M-9gx8L?a_wW(R^!rQv|'
    '0D$d_RDE+PLmAU3;Utc^5=9gntnHqQ{O5#kw&+N6D)oK;>fNSoQlIB7m71Qeh$j1Ect;Khg%vk+=0yTaMFLjH-Po1p`jlv!aL7aw'
    '=9V$YITm$`OkmN5Elt2wg*;MeVI8B}>(?R2a=$Wl9R=5LaiOy`GF=efO(%{}smjih-4Wc3d&IjLZFk3*qeDSkH>%vKx2E(Oy?UcV'
    '*?%>h@HQI#bY?Ehak;U#vQTCVz_FKCdWY@dzU-*VLf-E8hbKNzP}7=HTKm|%SPA7HoXCy_8aLDC#&j_QlK~LL$w)Lb-<4I%E&%TJ'
    'jI*p7f_TVe&DKOPn2Cl6n-aED`P1moVM1?yJwjKDhPKAC4eBMkE@;*x#^ouA=ROJAS0^F(Lfg}HA>t938Xy|wOexT8<7afbRS1B>'
    'S8g!HQqQgeOtdg0%zXD_&?!I|cHFJK;%v%o-VJQX@dP&JQvMIY<Ks(v7bo#L2L4s(luX)xrR9EHI>Ku8O}Jqd3OjTkBs>HT_j7nV'
    '?b!NGfRe86z<+(!g?Nyts9enX655QSXXUniE2gD~^=mwom<;4_nn%bIt9k}%a0xu8f9&|MC+&X1lr?nug*$uL=w@88bUhs8lNRa9'
    '-t~W1H7DcsPix<8_$+0i@`U5hfSy!A#0E08=<ALZtwYYvt?cF>X3)J&H9~#`5z1k2HFRt0g<Gl?#?5NM?{IKeug?&lrr%w#LDnx1'
    'iF++BPdRgR|6d8F;}(GVD>R<0TZnfcSu`g2&SStDAO@7zB~2t%^Mi-OD*Z1tY~@dMyn%aDUjxq5hyj;Ob2o&pi*rb^tIFCJOEl-{'
    'U=~*bFClr2^;bw|m{jtRXLPb{w&h>^tE_I<TopGc@s*QJ$j2L=OrWgxxSGGPhoL09FW%Lah8!b48JfKr^mOHZd#;WK^y<Bw%f69J'
    'v=~eu0xqZU7uVp(NHI7m`pnjXRIGhKgsjNTF(-nL$GTg0lIqz|qQMul9Sz`HP`9pl{px)`maqeNS`OplY~WlLILq4(Rg*=xE?09h'
    '`Z;6Pj(>LEIEU@_Nv=6Rz_{|_Xfq&GYNn{8uVKrbh>~ai(|>((jO^jO))C(w*=}7hd8@_6{u5qr_fPPUrNct#<&ezB(kw!K3V*9K'
    '<aNUGk84IKEq+6-kSr%xq2UZ4+M4;1SA8!dR<_T(mXNz?v84g~#D#i0d|G!TqnauO1f?Z0J^vl^*^taBRff#nt6Tsmpr51-{l!#r'
    '(dM1Pb3?o&r-(l7tTvX0cKx&G7ElEXa<i5;B(Zv=%3rt?eXZDQ>1q!kE7MulNPo#g7dK6^5mDIomqBGE<FBzF!cH(;Q|uD{dcdX1'
    'pKag{Qak1<@&ZSL&PecjyF&^aOP0K2p~h8E!N5MF1b<;xs(n5Q#f{Pu>8u{F|7UhB+BBWt%&y$v<;GY>+pXc^75Moz>e2!mp~cx9'
    '*!Y^3Pty+=o~iFJ153atu4;bEs?zV9-J78X#I2KB`J0;fzX~b^<m#q9Ap~o{4vn7<Y&ATKWNp+Yfm7rcAZ(s}fB6l!a7!odXO?N^'
    'j_yt0bq`%D(W;IyUW#%xC08KDMTpk4yWVHkX_IMp->FD0w%q1>w<aRa`}Oy0QUm1vrWlgRCv4x+D26{2t+Uie2Ph3=;B<{F^aB^z'
    '0yc&f1ZLQI6<$KWXGXKTe|#TjH`}0i^IT*&uQ-JJbvfLedGv`X03D#?k=s2Ac*ZuzM3(ynhnB}jz}x$f!>sm5s^sc2)-U4YkR+Qu'
    'eJM=7#XJ->5}_4UnN(iUO7B33Q=*9ng#yuCYHM>?w$F^N0AEvz%{G@2i&pxcDzrmc)NOJ;$WSu0_8hYKb2GUm?uk@|$$;jVBvhKE'
    'd@IP95r7SG(2lt!f{vruhwwUyzi4n`D{*}6x)87smQu2`=c2k*Jl>EKXXsc5h?kvSUe;%ozqn$9;JbN`?L9pxghqI0r<3Hb4Ih_9'
    'ZXRZ%QOeC<u=-Z)xHP1ZSxG|-v+BRzV{qdvlBu(JB(%R!z$S@C49<MO3krzU&%ez@Ox;DC;Pd<H+#AVu*LJj>Kx<V%2ENsBy;*q>'
    'ZERK-*cKRf7`bCA3{7X-rYqK|+Ut3PB>L}8=R?|%wVw}8C)ygUFt_~<uFcT>mss|Uy=a2LV5KV9u?-RYR(hK?nqon&Dsm#Q;ZP-j'
    'y0~bLOJ>D8;>9LgNv>1UJw%1wW*UcauoPQVNrh)NCmw$d6NBz<=~JW9v({%A@*UeC!h;Lng@x-(y-@KeNf<ufJo8R^?AMQUGGU7#'
    'L5khNIEFDm;?fc2K^WW<N~GJ`d2N-AukYhI<Xx4(toQ-{zFmvp4vr$cJI%C>=mQlc&YcER3~~{IAnfyB@@{r3Rp|~aw&0kB>~B8H'
    'v1Jbl*pO_?H*-sKWu=kqrbxjiDq4MMvF4F5<yh`0dC6fcJKC7gwK54*hHJzDKj9mUd_)H|%XeRZ`|ZQxBm|jAS?S6#8AMS(<$;v?'
    'uc%C8h<yXMSzM~rL}t65M|1@Q<FNcZ($WuOjlMtn-DgHnGGaJAdJSI6DV@E*k^WO?%n^gbN%b?)clqjA>?57t3$we|K=a1O+m#WA'
    '`ydgC6zCT45|xPAQwJXZwublLXuOB)1K-kLPju3x@abZy$))M#?xCeqyuY~6fLr+<t>}WFE*dRk=1ZnKizOgf>VtfM$so`P0u#%$'
    'a~%=)JL1t_?Pr~OIQ@x&G@tS(?aE1p3$@;Mrl4-&QCAZQBb50lY*G1?1<>$}bn33-bwAZ?pke=H)|X1$(>)%4zXnGi9buq1T8@%U'
    's-Y3ZO|qd?hcU>7kn^^J^I-P)Y|pS0Dt_3Oznglu)~t#rD$;iQM?A;2+&b@LhyWiywI&Uy;6{5cDxETAaiOtuza)*((A(VHn|_R6'
    'Gi08%bjQ#A{l<A|PkX`gMl84Gq<8tN0d4*j*f0J7NJ6Y=YRV6YZu(WVuq3ohZA(;1`4gEqPQ?9hatHsUa%Elz1#mfajD<&4O~Y1m'
    'reSV_dBJ07mS*Zcc$T6LgS)KCR%EzjtZ{Mj!b`m!VDH^S7Gj8x0IFaED5_prmntNU7+&Pdl>I*5(Q^Yp><e;!yH0*$X<;NC?24f1'
    'I<g^8KA@cwomZ!5z>$FNYSz9RjoF&6HB|OCwt}LdBZGdRl@Q*7v!?Vaif>$(cNqqY^J#0gC_Jm!=|Rti6Pu>QOa`gM`jqxFCU@E8'
    '^^#htaaf59BLpeL-RDHzPVbH2f)@r#iG#UIdhOzY&dHDAXqAPo%4l#zY3gy7DmvkTK?RqGuXRBatCT4`m0lmaOBp;@L~r#Np)$%*'
    '$vkOI*jAqg$Jg|zgEqhlyaSB@nOJimEv;Q|5waU6<+3jWJ)Ob*2zbu2!j&;A;`axc2`KW98*M%LSVA^WQ(xI`l9h~*Cbea*c5-65'
    'N?*<QsZHILT@w0jxa{<ru1B-pJ1wnPI4?j})P`7CBt}T?Yan{NF`v%aZ<IoYXQI@V>tj=77AQ^w@#|FH;tv&8fE7J^jKqG5Z!kLi'
    '+P2PUZz(vc?esOc+9agGn?E<*MN<6eaK_oc0&<^w!eRPMeAc|EA(7-|QgS#rbmR<1`^KImHmBMI0VTHj^s6=^PH=ow@z*5@pNUa@'
    'o?_EueOFoYcv7#=ek486nI7$fEK#v%C0f9htcf-EGbhw0`DV}LY4r^?C^wLI>Za}KRSnBj#5yY;mE5`~TlgVg*3Ebpt~t91RF-VQ'
    'r*@>sr;ATd{51KkiGzhY4>24e3UW50Q)BatQ~!E}#6ka60Jt(ZF9E-yGl@s$FCxHPnT)k-zr%NagW|flqwr9k$;x0AC>-35*4bt5'
    '&-LkM-DTsK)^HRsHZajMF6l6p%-odU0-pIWL0pQ77MRROj=iC(Z>8=sh@jO>7*b8bB}!^W|BV~Ll9z|i_Er|Y4llS@_317^sm)S%'
    'P1=}`BCe4O)<HOQ1(aNTj&3VFIq#n~@Pa?GS{`#q27Yzh`4BclkZ$?GI8eZP_9N(HDhK9>es+D1_X_ES_yIRMGR2f?bhKnLiz<7r'
    'a_|)m0u>Yv)Qe(9@EDqY)3)InpT3F0#S$l<M=$;NWRcnQzsn8^sb19LX%wJ~)0KRE(s3OIlXe2n!=vz)ycD{MZ}dFJSP#SZIMvzn'
    '52FT6Y!Ow%AN<Q%7j>22h7?cCtIqdVXqy@t4Y6+%xdEGiVQVrz_}{BQ7g^cihqEvWnD;4C?0Lv~anAgDs&g!zkVfLoXOaGp@yWxt'
    '7;jP>c@q?H*<f~@Muh(HU^O4BtK}qNTAS2Bc6XQ(X`$FXR2OA?16BpZfg%&=aCVsR!ELWz;_@zVg~IXu-XMHa*pt|Uv`_aA5#_|K'
    'XqDn&0TVvSyeL$TnzY}}oTEP=I=CYdxwI~U(v9G0yC_$zQ_S&<HFsZZJWw6?#X<&z%=>EX4*&oFF`;e>65K}@00FnO0lbU?zK<O&'
    'vBYQl0ssI200dcD'
)


def intermediate():
    from base64 import b85decode
    from lzma import LZMADecompressor
    from struct import unpack
    count=71*121
    decoder=LZMADecompressor(memlimit=64<<20)
    data=decoder.decompress(b85decode(BANK),max_length=4*count+1)
    if len(data)!=4*count or not decoder.eof or decoder.unused_data:
        raise ValueError("invalid length or extra contents in the frozen bank")
    values=unpack("<"+"H"*(2*count),data)
    m=[v<<10 for v in values[:count]]
    p=[(1<<26)-(v<<10) for v in values[count:]]
    if any(x+y>(1<<26) for x,y in zip(m,p)):
        raise ValueError("frozen categorical law")
    return m,p


def run(request):
    if request not in ({"part":"start"},{"part":"terminal"}):
        raise ValueError("the two parts are start and terminal")
    from pathlib import Path
    from fractions import Fraction as F
    import runpy
    evidence=Path(__file__).resolve().parent
    if (evidence/"files"/"moments.py").is_file():
        evidence=evidence/"files"
    cap=F(500,3947)
    coefficient=1+F(1,2)/cap
    if coefficient!=F(4947,1000):
        raise ValueError("declared coefficient")
    logs=26
    q=1<<logs
    half=Gaussian(1,logs)
    dm,pm=35,60
    T=-10400
    anchor=intermediate()
    if request["part"]=="start":
        m,p=run_layers(q,half,T=T,dx=800,dw=600,dep=15,dm=dm,pm=pm)
        if any(actual<stored for actual,stored in zip(m,anchor[0])) or any(actual>stored for actual,stored in zip(p,anchor[1])):
            raise ValueError("the frozen checkpoint must dominate every regenerated row")
        M=runpy.run_path(str(evidence/"moments.py"))
        G=runpy.run_path(str(evidence/"gain.py"))
        E=runpy.run_path(str(evidence/"energy.py"))
        rr,uu=M["compute"]()
        cost=E["energy"](rr,uu,M,G)
        initial=E["initial"](rr,M,G)
        if not F(cost["energy_upper_1e12"],10**12)<cap or not F(initial["initial_lower_1e12"],10**12)>F(1,2):
            raise ValueError("start comparison")
        return {"part":"start","checked_fields":["line and repair","initial majority","bank origin"], "checked_bank_rows":len(m), "start_layers":15,
                "energy":cost,"initial_majority":initial,"part_still_needed":"terminal",
                "graph_cap":str(cap),"coefficient":str(coefficient),"threshold":"-13/50"}
    depth=9
    standard=Gaussian(2,logs)
    law=run_layers(q,half,T=T,dx=800,dw=600,dep=depth,dm=dm,pm=pm,starting=anchor)
    truth_controls()
    literal_minimum_controls()
    terminal=root(*law,dm,pm,logs,standard,half,dx=800,T=T,dr=680,du=720,dv=760,kd=2,kp=2)
    if 2*terminal>=q:
        raise ValueError("terminal majority")
    return {"part":"terminal","checked_fields":["all-plus trajectory","strict terminal majority"],
            "terminal_upper_numerator":terminal,"scale":q,"continuation_layers":depth,"total_layers":24,
            "literal_site_controls":9,"literal_root_controls":162,"root_payoff_controls":72,
            "graph_cap":str(cap),"coefficient":str(coefficient),"threshold":"-13/50"}


if __name__=="__main__":
    import json,sys
    print(json.dumps(run(json.load(sys.stdin)),sort_keys=True))
