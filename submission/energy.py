"""Active -7/25 first-pass, gain join and starting-majority proof."""
from fractions import Fraction as F
from itertools import combinations
from math import comb, factorial, isqrt, prod
import json

THRESHOLD = -F(7, 25)
PREFIX_TABLE_SUFFIX = (
    'U)*38=M4Lr_9&icjkr1%IKqc%cFYKGh8y71P(9q<uLFxl+KB$Bi77P7gS?-TnfDCZsc?Fqf=&$-'
    'G`l7!&Jy$uU=C(Sh(1Cena|w&_vahw!H1cHz3IUttI|_>Z#U(=ZAuT;qX!4^-'
    'tNwssn5f1qAnA8;1cTMUGeXF|L4s35UwGN{$pnN5|qELD*aB|$R=+snfK!pwC6oS<m!vqGvp99eA<fHAuF)LZ#F)'
    'wa<O9GaAd#9K!kQOrl=C&I4>5jUp2+a_hC4?F9enK{gM2+8Xn|SM!#HFSerN@UDFOTv#rp$nK=T6n7~?Pgv;mk(e'
    ';rov_f<+>8uufxldyUX&`e3btzS0h@%R_rYQJkt)S9&X1y+g(d-'
    'DL9|AA7@%c_+Hp}4qrw6xU4tApl7x02P!KG+G!`NE!@4?4;(0I<wFdo*L+SpT@bRPO6&%J;!E+Uk;{ueXzSDVzD|'
    'MD;Y@-P4LFaPo{|MD;Y@-'
    'P4Le}&=y^513fzx;O@L2~}(zhvCM{CAo9FaKQ@{tsC6|8sWyFS7IhXB3ee72N!vaO?l$+$ML(J@Sy05ycNs{Ris*'
    '0*(Ixn*S5D{x?U*@1NfpHz5|pp16`Kq#Nm8L1G0-|3UJ9(Br?M=YK)3|AOBC8<9f#kbYz^$srTTEVB6T-'
    '%=vQ<TSbA_xJy!fpbM}lRFjMCHE@0Paah8u!2Yb!Q%?b{u`eB2T%WlXa6toihLyHf9J&?uq5szm_(AUWDuD~Hjq<'
    ')zqd}tlJWlk{H^(kWD=Q7rjV&*I?1nK2AN4_k=bNU1#>Hy_xF2670mw+7F4jXf<^zq;{S(O>Hp8~J@2hJ>wWTt7z'
    'Y5=NEAsYGe`;f81T>EtLRUPNnifwOfrZJC0S$y8AZmBabyCSNG6fVWGb0P@+&^?OfrkiCUeMKGLICH`4uc63oBSe'
    '7FV!@EUjQ!1<U_$u#W5^rQ{y@MI3m}dL)$;kQ>A*h`vA$63<}1m%l%!T?^Pw+J(@ELU~v5btU8wIZjTK^W+M-'
    'PHvO?<Pmv7o|Bj44S7dCkk8}``9^+_pX3+${r3k7`L!kJpPy+}5jD{unndgG?|cy*qD%CMJ~61kkQn_hU{8F4{`n'
    'oUZX}Nsk$0pn&pem>B7<uE{XsQ8Q&N{_|3mK8=DF+e{B`;Idfc-x?oR`rry-'
    'E`_hU#7nLwsje6GJbt|IHnR<etfkVE7cIYrKqi{vV~LGF-'
    'y<RK{|Pswxgio7B3$b0gUd?sJWH}akQAmtTj_)Y#0;aL=-'
    'B5I=XA6SxVB${NC6U2&V9z(Qw=DDN=&#ljMe<80c&R&~mPvhC85zpKh=on6&n)0>HII}H)z!>02EO4L=u$iB0KRH'
    'Y8kQd}DQSoyd5*y-5s*(UwhcqV5Nh=ahx{xH2LNZ7u8A7tiNHT^@AbDg8nNDVq*%fCfAPdQ2vXm?%D?<PI-KRCAk'
    'gOvI$rs|sGp{CgJoE0qKX=1(hx6PX_4)pI_8~m`JwEq7Ji9f|{xX6(^6b@k_FK_>J<r~qXU}fW=g@)AurpuZof`G'
    'z-Q9<0?#K5vh`zwjy^WkAkH|M-'
    'z^~<8_n*Im*O;^+NhFhuBvVKMSw%LJz2pcvL#~iJ<S}_kK9FzZ7g6#3YY}~7LM(_4aUjmb_3!zJ7YQQ?WEl~jc{B'
    ';(xs7=4Zw>zWJF-'
    '1^c1@l=muDXnN$+jWeU75eE$OGNfwk=bGoJlOC(da%?p+Vw$tj$TzVwZOoaNy_#AxQ0iF_~lyr1UrJuatS>wtzk_'
    '#XK6dXapxliVQX#GYTj7HLD$8~n3_OGptpP41JAM1#-'
    'A;_uHg5I+(|qDVYRA{k@|8BHdUnPdT3P6|mWc|<DltZ`&5@#dM&ku6PtlRUczXMPRO-'
    'j!z`%Ci@8?i;oOoZ50n;+QuQs7Du|dJ<pXiybG8*(H--b12U~k~)m%Ii>>F<^aA+c&<Y3<aXw<{hZH}e74u9VHxk'
    'AkDMP7Trd@kauMw3dl0^d+GHTvKpqiOzNZF%e{YS<BL~SdqW5<`Ni&i{#*w9D54l91l3&D_=W!$cq!H;t=8$WD|2'
    'v*@D$(Se&m^8avk~XrndeUB%s1h4cjMV}IQO2Mdp|z=oqYDOeD(#&{CGZl&kW8UpS>2JeM>(3I(+ui`0NMHr+=;B'
    '%x$8N?x9Z|<9RO9gYVM=-%{V-'
    'e65jSawTC75d<_9%;+RoH$ZTDvY_fJYI8_X@2<d(?`;^lM!fkR$B+xeneVkH*+AY9AHMHovV>eFT0BD_=|(1z9pw'
    '7|LDF|ebK(Er-'
    '{<pT?=5?Ul&FwI8d9l<Qb>_XA&DqSOGBlkLDG^UZ6cx}iK4U+4GnE2m3p@{e$UVE`^V!r9Ot}V_wl&*KJMe3*AOM'
    ';{eRzka8Q(~z_?kIhtT&uQ3V2TqTu>`fq99Td7I$;W-'
    '<380=tR8J~Tw|ZMOJb=)O5pc%R_@459mt%SCorD;|53@TeVv8+*m$92UGeEgs{7&`Yi8SF_;GOELfdfC@>0TbrSw'
    'H^bJU3~EytqUJMfUe92Dkm1xthVsV@+J6|{8gj(CaRi5PEZN5Ktc=6uH3!ooc!=lzP*jL`?(L#N#q)nG%3EN#B+5'
    '+SND=i~l$*e_Q&hVsd%yqFvx%7NO;HH~<JhUf2L#T+g5!gx3!MtQb;9Sy3(U8~%$Ez?2gKa_3f`Lt>=3j65i9T{2'
    '<-p+>?eW$wZOk$;NK<i>kIt3g@CTWzh2<)C-'
    '7$r{QU*~1yT%WdoUdBD`xD*kPyhQK9<2bjp6BWhCy`<y6+k8C~_DL;Lx1Pac>RB(;|+VXB-nW2q#AqyjBvj%LxB{'
    'BRCI~z=CBG$g7k<l&U0b1+H#UwF29!ssHyL&=;aI1ipR(-'
    'z`z=1x8PSafr}#FM)HAz!@)m{k*_BRA8+Y9v3h0P8D3Y6+WN5K=e?|+(>Bth~RyWz`kD0{(;c`VxfJ5R54qD|DoX'
    'jV!{8cQ)17Q3*Od<JtX*lN$~&DXR#wV!$eJnZB`6(#QgRB7&b*S4A{)DLE!%*@IUy#@K2Sa-'
    '(Zf)nH*!ba<pCIQ0?MaWljj0Nhsb!aCk;AGLgX32nqO{mB2zNNkj)pBE3WsgEXXI6EB6wtx~A=lE$s`(&#-'
    '%1|<UXBZ2v|z^p6q`Ut%Hgr|KKm^TP+%SDJ;3e1Q9_nig;^K*e&Phh?&FnbEjZbI+FR|-'
    'F0Eqp6kaAPCjE3n_$CiJ*d_{lz@z1;u%Y!UfV!O2SDJN2SB4}_<`75U<)*lBVM@<t4*LjV4r3^PI*(vld0vKc;KV'
    '(@*=@I{&<+Log;kYnTqj@K7C9(?3b>_u>yLs)m15Y|pO+)o0RS4g1ijs$YdB(Z0$B+_0;;^BBHG@q3cI+I4z5oz2'
    'rmBEDLGLUwXg+{$B%%bFQQbiuTQXa3DD8SHO5l#ZXp}-'
    '#~IKEE!dzaw(GQo2jfn8f*?<25>2<&>o?=}kTkA&wHY!SXB=Ke<b{%65`mxCf}=L-'
    ')K^Z!ry|B;J=Z`Xu33Jx?s75nCc;LabhuZ0J;3jbd+NNC=hVORu%=>~=@;eTbd3?+gC?V21p!#VCoacs@uC~xLSR'
    'wDHGB>1Eep0*G?jU~VqN?=2s1O^HGit8lN)-'
    'H)jzEVi4l@fbi8V+^RaPyKu<zpHATPBMnJvo%r%i+i}d03B7zyV1`oPMZ?ta2rMJEn}+St{6lTorpytD*g^I(9)5'
    '(?uWBMIT(mzB3k@HxTnbBzoX1bnm`Z<g^TtCwB`y<p3H2|6?)#n*zVI@PPMM#NPx5UfmNJp-'
    's$R^g#R1|6U8uU}enECU|gsoZ!P822a6*fE@8t?0`RSgf9vmOn2c(5Ixv;kYmn6j&5zj<N$(c4#8Y-'
    'z~4y%TeeHUszU;OCrYCAoFpDwNkJ`73YEsvI9Di*GlKhj8)fh+P8Lq4a`1Q{C!Vi7s%9u)$S6fjHB`cJQ)T?pR6$'
    'cORfG&x!_0Z=kUpV-#a}d0<)aOy`W_-X>Z0O}9@5qtz;=QWcIum;W0mkc(Thsai*BI<-'
    '(A3>Y;lJQAKW5%V0Bj9iXsO%)``1A@Sxvgv0H=|j(rumf*7u<Gh|vYBsvQpm?StL_CZq$!}((j6?X(LzA;QU<*1z'
    '?^stu0?}E@lH^(Md!oOrfcLU*j9|=rYBZ1Jj643LMMCT<*Tyv0uPl*)D1qb4rq)`wdgUNcb*jXcs_N{U#@{~uLxd'
    'JSu6>;dRA`bsjLZGAy@-0;Haf%wu)77E&LjyV9S~z@H8x*92<UhJNcU~V!Nrw0~-'
    'x%?crntSZCthB*!0?Z~5D?S{{Bd7wnBNZp4t59-+;9^8m?ARKs-uD<MWPQv6MkYZph<YF@Wms74=vq5kOD)T0mFa'
    '&MGr<Z&`gHt6%4EPGjt0+OnS$lqRsJdG{=9-'
    'IO@tdWPWgwjR{jW6KsSZG8YNl+b4k@GLnL)l34biBn($c3EfCzPM);LcQSbVO9t*2WU+sx9JU6^qrbZXW)D`x)=^'
    '6En5vAw%T*A5UKPE$I<(_7p!!S`{lm2}LR$yjwYu0|qz}zBL%hl|Mr)lZoV$7=ZJs4A)>>h0v^8pUY*G2GKTLku!'
    '~UNmrrJB>Q@$&H{Wlz`m)-'
    'Hg(gRoJJz*xavQ1><8^R~Og+_Aji2W=);zOs<$Ul(*MMgYgDg42ep>M$dy;pK{CPQzL5stoQ=p+1L^H`1rt2uOQI'
    '3nbP7N!#3oFIhCOW=5n=*3eBlmtq`;)x_mBBf9&BQ3I-'
    'H0}s4<XXw%=XY5Qd@KjqYI#`QRzP8^BHl<SW1qVUx>Hs0UP2uXOEvIKMhiYywec`s2TK$5FlV{}=Ik`We^n-'
    'MH1CNIhb^GyWrb+6#)=;{`0sUpsJJ=cg8v|_dg+WKtA?W7+YQ@Cj)wW<ad@0M0lLp8;a@L*46B$1HJK1tZV1KDYv'
    'K5&5`}4kC+CD_ylx54eIj<xC(#$7l{vZ$GIk6Lycld22u^HgICoC;<1NErJ&q<n4!fNk%SBdD=}#ECi7>m9kQO9?'
    '?j{M0n<|M<ZzWNaB83t=Y4~+YqwBm3y7tJTe!U!4Z<WWXg9?amR>VgcWt>>5g3yPmXpa^?sI7@}&$Zxkvj^g?>q7'
    'OlJ_0<9khjqU!$<XmRlNnO&st&GQEM39wH5cE9ljqPh^JErqmQdAZuT32sx_mqt#&Nl#7uy{zArv}^2ZOOAb4ubf'
    '=5Ug9xs^>&-'
    '{h3d>)S}_m^R~^lHe|CBwX66CQhQMaYA6P_@YCLOY4C#U2s587O#?YRxchtk8_`4%c*sEtkdq_$qu<?9ihN#V)N7'
    '8B%1Y#fu0VpA!ZIOW@fH31lr5yfBi2Ym*dYa-'
    '^Y=AcNoYWU((!4x`iL(Q#J+O$JIZT(6ATda9^8qlQ)K8W=iT3wslLV19uvvLy^qvdu`~GsO!{b7A(q5I3|BGK*}G'
    'Y&-'
    'xPo;aYh#R)fm4Z*hX5xD<&6kJY?LtM2tWbRGIs%?P?5c7Y1MTC%&NZj4N5bJI&!T9|v5Yo05OE09r=Ezpix}ErEx'
    ')&QC<iN@*ALG(aATIwjjI&CS9(Nx2cL{Sv)@*D4-'
    '|P4VM|2qW4HQ{IWUSsRg|CQwb?zC%EPamqB3G54;V@JsTwX@V`$4$7S_0vwl1R8MiSfInz~@WDG*kw&Vuen2%c1W'
    'xdGryxLgk?nY!U@8#;W0$hX$lWwV-{m2L_t!;bpx64jeIt^+7Yt_+<f$iM?^=S6_T=?~l-'
    'K2Rv9Z7(ZKGp?=y8_K(J3kg^xtUrfS?Z2_qI5{#Z^bFl4EBzEkNMeNe0_!OUn#)TVD?3aeUaXTSzw+}_`xyYzJ2J'
    'f`f*!}S=Lfk4aF{&Cb1FqxE&^rhyzlZB~kMP9r84{<y!r(sdpnjnf?48(8qA$z!gr7Pw+?>W>pUe<a%FzEK!{vS)'
    'y;g8IiyiAZim<kl5HV2#%fE?iRVs<EX;O$vmWJ*=8B{dLVpv~!Of6SHfu9m``>UYbNezD%YT(Xe!Hp$4aPZZ~v`{'
    '1H=9ogk(E<m5Sz)!64RodtfQphM!ah1<hxu@9y*LVv%f~}LZ6byl1YpONV2rn#i{-'
    'XasC~N#?G?*$V#hl0s#M%4%0MReqfciZUUZ&>!{-'
    't#%DRXJku}Jys>l1mO>i6i2o<i+ktW}c35!0$<>h~n7YoMt4r4h1((K<g1*XH)SnyjdHu<6+8*6XOB8A^}3%<M+d'
    '2FNT%NWs@<!2eTe`YWl$ni6o<M$h(Bau5=#h&Rc?v2`klCV^kLW1zr-'
    'o0foDM}VzFXd2`ssNV|!I2dzD88fytzb>8wC(|=!Fu>|*bp1tO`+vxf!xi#5vXd5L$>xf5#od=Rzq>x#vMNU$6@K'
    'tiCFz00DCuupfn~Np?71DzilZ-j#`65D>uU~Ap`$4?T3r&5mYWY4ZDVO_}g^_ex^5(M9m00@dO{d-'
    'ymcAM=ZSd9UU3}5D_TNE*?>06?ZjRRgpeBG1`pr5-'
    'XNaY0Czlb6|^yJF~pQ!&pn5JA0ZxmQ9^8fo0y9#B#M5#yc?F_Gj>3&G5a5p`%0iO@9toahoN);E?nu3~D3XPnN*`'
    'F_I{=mclG2X^dMa1DRG?)UJ|8uV6(auT_S8hbm-_X~64(HuTK&@Ug=XZL(%iK5B_g`>gRF20;0@BXmZ(;<uSQ9Jh'
    '=|$QWNV+?$3EvqLemF$(>>mq5{EH9Vs?BmBlr*i6a6!uQ8vlT!-4gv&_ybpzARH^cnPGiZjshq>E-'
    'IAQe<tIA|p^AA<_S67dnxnRb0##^%@KYO;`$%Pr;bYt)Tjb$TMCNksw{%o#MFl#gkWj50D*n^B{R{m%ao6HkgS=C'
    'Aw6tkA4R%~F)`!TqQ%x18R;hFH4e=itbT8Mtcb8K(n2=XAz75OY-'
    'kpuz;OTw|I6b?E_Lw2bQ#($H=r!(@nQ7QVOqJj_S)lgWkiJY-IsCGAiU#<zZrdyzop)~^f4FJ162u}@$;fUH8G!='
    'Q{!i@mDD4qrNSM!m4b_sM_)*yXB8eTi>L5zMrvW6F9T-+r*%D#bD-'
    '&#<n_!5iQC#<;KjXM*hS<P`ZW|wWiKFsXJG8zUjuh%YY%AQedmb5n;`p2IYUkqVA6eF0b-'
    'a@wNQUaT*n8XxyHn7<9TUfe92J7LphsCJnum?L1Ggg0`y>BjLt|=wVO7$EI^sQj4#ck&(GE4p;hK7d>bz(p72<Kq'
    'eIGztA=wBolioLwaMBHt%Quu8sjp1`;5dKjXE|=tyBC^PV$tuXQQ-'
    '{=AEu4zh#o<eaaH=(fah%w*A8m2Q+!3oV6u%CShO~zdEMfwonmY&Ax?*5tzXC1yQZQ6%7dG|EgI{HlxGgF%YV|Fo'
    'jcvtn^><iR^&JWkl1yo<8Z$35VpY?u*`zmvSl6?WY-'
    '@%$`_~%CzLd^kvEyQx$;3o9%w!$=?YV``7@x`hRAe*vrhG<u1+34<Que7w1v^(z&3c*KV9|Z<vK23y+0Xf{?BRws'
    'mY({Wd6~UuOR_$)vDIIgp2#l8_cG)*FxaVy9X*}nLm7v=rO0q1zkMGkfxQZnC=mH|lA<&WgJodaDGT*_c`&Ag{s&'
    'YbTcnOzqqMOwNe@m_jN$mj9G&B=F=v+@icdMiVc$qR-{J-9&;aPUgyQG=7_6^dfy28u!{GUD+^af@5d~$?P^-'
    'nP^(~M!d4uWuzTsM%BzvT;$#k`QvWUh3tlDZgJJfdq`|t8Jw&YznYu&e)ovmNZ%F4E|BhtIsf$m%u=3Bt7PA+Ftu'
    '2-=)FK)1vo156-idHso&r7!C-+MMg_X~@!{K-'
    'fg{J*|ZJn)%3e>X{$=grpQXGiMswarGn(5xpPKc*MIA~J5Dn+$UmIIfF4_o;{@%}8Y0TtdFP1TINSLSFRcv5GVjB'
    '4w~fUJh%afHwica}89nI!XgtQ+gockv>YzOwrz_7lKS|v7>7ssw;=WJktY}<9)HzZU)>RM`Cx@Qta&7fYH{wknk-'
    'Z6GoN6=4~y^8y?_C+dI^K>BdogMOF_(7In>r&Ga71JcKtY+Xu6P`B7}ctmW)})+W|3J&Ub<c$lReEoQ#5l`Jj)2K'
    '%Vo!hCYu*t#$8S;&F^*tz?E+1ZoQJnxM%U+`L+`;9Z=4T~)JF2lZjWUL(@xyFe<)*s48Y#qseuJ_<iW_a=ZhZA|e'
    'h98eM3FON4VkgOSnEG?*i~Q2Zj8K1!@X%A_IbBKkD~KE?d{^e649-'
    'rK!?Y*`q`y;w=W|s&o2?1GMLKAlXb9OGW;o;88|?}G(R9}d*SH({xOn1TnLi@*LUCSXu!$d6!P##st{pyru=vw>_'
    '^}F0YMU^~=nbSYf1>X=1*ZPQklpa@&m>FSSl+S8ta-'
    '*<wk2!{8{)o^)g)vwd*h=lV{jRZJAaj_6gRP>U)$K8)=oCB^A{VuQ;K(NRps&0`rHc^JVAujZ2^OLpVPzmu`^@%S'
    'qE>f=Hth2s0Z=20kiqRc@bPoI+||@Uc`F`Eaip@EBV)TYxu$J^}OcLCf-k(qkj-'
    'b`$Z0&0ffPq2r~1<EjUaPw?|2#dxbPkcFQ2{rW`(eQ^4-'
    '6$~b&L4f?KHxR#^~S&>z=`<UZhK_3|QvV+DUXAIDA$AR1lsB;a(sNP{1_;V4mey)XjLprpp@?aNu7QScdP`9iV8{'
    '0oZV>dCM6m6Cp-iL(_8_IquPhx%>=CEPM5}5bg&FpvbepV4w$c|3D%m&`S%eH2=v9&iov6gB7n2f0+AG%bR@BL)K'
    '(`4<qg}W<P;U0Wize#*+(KIggDwM}giQ*kg7IUt!g3H^j<F#xv?=x~cpE`LLk9xL`rwuv8=UW`%5=|%g<iV$Tkyi'
    '<i^)Kfa=W#r^BW|LxggCLM6Eh^RHeM3L#ce#fP8xcXWTD$z9{a)+Vaimn)l3~vuWRA4vL0OYjG=Ya9JyiE2z+b@u'
    'Rs@^dg+enIB#fsOv886aM%oq$MC7in4G;6$wm2C*HDg8we?ss_8GipeMMHQ40}^%!0yxxV4I3Zv-'
    '20Hv6qWt*b{|hHd}o+d-mug>(gAx91<E?<&;;fV#7Dq)<c@ldZWev^t9m57TI%Kn-'
    'N^oaRR@)A&@7y%;nPz7xJNb%lPD($-HlQ8n+GD#WOb@;Md%a@VN^M_#&TDKIgy%9(t;ZU-zixT`4zt;pzr{-'
    'LZwgx%h~iH$LSns$X#ailbt-'
    '*hhawwmT;Qjs21cJ1vDJKcx|{P8Kg`$>V;yA}%|sV2HOm^t!b0p}!t#w2U$KwmD))S>tVi9U^+Tpm>El3RJw2b2b'
    'pZsbToMbTR4=uEX(n8Bp@cht#}sTs~Tl_|?x4VEP}bG~`(R4kKnT!+}k69>;1-'
    'Ls<Nu#mqi?GrRXThi$l7%Je4PU<Y-cvY2&UY>uoX*JRpUdtfhq<M1H<C}j+<p5Vuy?+E2<<zji-ot1o7<0ih^G?P'
    'z#mc#u!kMo7$CERj&1rJWY%3GpuakH{}d|YuW?~;DO_gJ>`ch5U{^Q5m_WBN}nrSO+eSix!b8Y%J$lB40V1QS&Wo'
    'V_9O|3BkubA=RQrDU*ui@29V<S}cjBJ>SZ@X=fy4#is7*eG;#$Ou1;%`tI(A8fxg0B1i9#!=&ucoXS`ksSdzv}_L'
    '4JmcW(utwZhJMemUE{?a9;_vG^JpA<-'
    '`PE&xtR=%*UKy}A+V*V2ma%NUbqMQ~AJ06Tx3DQex$MWRbL{Wk+iddBmrQ%wPqu%uBF7qIzM9+dQB5QHs`|-'
    '%rO_OIf8!$leau=OG-^BV^vUM-FOKsEi_Y>FnN_^?VLh+OZsxHSPx)%sxBSY8Ph9WxcW%`7SMWrVJZj`AJ42NkzG'
    '%@SC4DkJVNBk7=Co7OiiW(frt!)WP`fUH`oog&6MV7xB#oZ)WufCB4;c?d#MUU|T!9+Yb+pi|uM3|uhS;mzQ)HIj'
    'NV4vaRhmxF{5u?-'
    'zsI3K#}8gnv!L}f8oSo7go{@i0>>P{&PhdB6H<c+mj@yjcOq*#u}LvHOt;>aIpvIEefLafV~@nKSEn|yu3ZP2i&r'
    '`8cj`9#mHvuZ_UdNwzm>V$)Si6uOGp0q*;qc>a5_)j7{%{~t>kC=rSbHL{amBq1a}y7j^FEhjT_oF@`&9}xa_EQK'
    '6=a-'
    'zVh^MZnII6Y+fkRy=+a|_)DKOPMcA!bZ@GDXiJYLIMDSugK3=7P)bR3qkf5_$z#+wI$b1zZJQ)<VWSjEFG%B>l`P'
    'C!<q+~j0s4cLu}wk^D?BtYMMDQ?lMFyROc4{?3&#4kIC*X$`i~n5sm?KC=T8C>g3*>9iHiEAST{ZeTGhMpYS{^d4'
    '!Q)lt_DOdeuIMZzj5`V8vE*Q#oF|Tu{FMaZ0(sSR%f2fw$9qeR$My6VsdXVi+RtPWZVzdGDL~%`<w9xpB(wWPviJ'
    'r?O?9IJ%&%8vW7R!P3J=kbNNhmhJU$H#c!3~;g$Re@8k83FP8ttWqLBYxkH}PBQ;4n*^p*yTT+y?Ek&*vNMGg-'
    'q44W&lv(FN`koVLSAs7|+6GX0-gH`07((q!LP>g>B<A^uyc;KtklQjaiIl?vZv|Z0s)QW_RncR-'
    '2Hf6g<EWfI!aIy%yVL^HnyoQI)*h_`UEs3V9rMa3Ak{1o%TLb5(S$`fvuG_mHl*WhWgc#;m7!+o4XmwsimLNp5w}'
    'B*jp;IGWp5o>ox}vDyK@dJa$dnA--'
    '!%A_asYitY)n95!+kU#a_$Ha>Z|kJgA>NH+eaRuk4z}mpDc9fG4Z@wR`D2MfEVBx~P=*GPugWDc<APTwn083!nJ<'
    'Wq){Kwk#?C)gZG1Bg(&GMKQW|<p08%f)z%RUD-GqTs4u*{Q}5)WiT}xhf+{T1SL8ypk<e1sq0BR`4%VA0!1l|H<i'
    'Zya2c%ZkcD!YJUSZ{(e9%HMOSsS?$bhdnl21<ji3_R6Ypm928Q&<il2iZePuZQWsbwUgeg!<o{7NH1+W>k0%6su2'
    'wA%yx2F_h=kP0N*KfgvmiN%{VJtLFn=Kk|!}Jr~+0nLXEb>?!bBWx{W{k>ZKd+x>lfoKV{l<69brIu}i?q4Q=Dz$'
    '`>j?gJlOLaWd>#+~x`I0#+{R}dIK&?`p5aQNHC##brPTE$@42Omr;Gw64N#;Rvvp|VpPm%nYD-Id52li0H(FLap6'
    '<y`CWGhGXkFjgl=dcqR=LKI2I8qEVi~OpNuotR*3tLDDWq(&g|^g6!S}Z`zK6@gTvHzZR28vwt}<?$sA1$JO=SP+'
    'fgN53=pJYamp7J(546FIVh4Et9fHAQMx*tN59ZlS$F}F;*jf~i%(KaO|9&SD`X9l+>~jeEdlS>fJja`5-'
    '_hZz$ReJbvHAB0v(K&**`w>>%&anr*=A?4RcWVLm2Mq7l=qCCtNYI8XDIPmuY2;cV}p6@vkCmss@YsWX$g0~wTW}'
    '@{k(WoA&*;gnWt(raNk>PeCM-I+|>0SfBR8^zI@W5Bo_-xebS!{TU_X~^Jx0_-'
    'J4RCr;^OJne=9R1X)Tiq(2{)(7)L$sczyrN`IX~H?ftJFQik3P8R*p+eZ#Nq+#}226hQ@=r>UTw^u2lSY8z~m<Fb'
    '-)JE?Wdf3;~7~{gsv1WB2<V6j@$1#IZq2dPrSL5;JsUP0{nuRZpF_^t`CI0qH!-'
    'MzPh^je_xx1=Svfw`QEk0mE0cUG(^<aH=*s@`bquASsAhu3*F?+4Lg<X1=$JTUQWY6Z`W8GgmSm|%#Z%*rS<yZaq'
    'kBy`FkAs2ym3lPa{e2BL@0rPywjbp!^Um{{rS*LE*H+%XzJveV-OU}^WT|?kHZ?8mNmi%&QFEvZW&4k&1=oGZa$O'
    '+ZIyRd;Jm%Afphe{OVHr(%zlK~lrjS$iR`T)MNlokbP>TNnQhT0D`+pvx<@05*yiFF1_RAyfkRqOP75LYwA^)EyB'
    '<plAQpFGo>Sj<q-'
    'wVYfY*Bw?AdEU(Q8jc7tg<HJRj(lER7W6u;}Serw;t&SGO^&z5d;{XNA|p1P(S`0k1l`5<sFL5y4{TR{yCU&pNZ_'
    'EOazN8PGY;SX0hl~r&(849V=04V+k{Vu$$ABxzamx-uBy>o5uKX_|M^c3YKy`|1I1nIENqJP|Tf&)^LSh&HPF9Yp'
    '&e>AMbxbg8Fo*(ZDccvLDx%mM?N59Z7eJ`Z9rrjh{-oinB@5az1%9FCzQ)<#c(~I$~Qk(+-'
    'Cnv@#@1?Cxy39h*yE!;jK!eu~Jsh+_R@@#(W1ezz##vz#)#_p3tTpawE@dmyNfKK9lcqsYtx`2(!sYHEi+p9kYnn'
    'H#=s^2E+WfAmk7jrWIQAoFD<>SEHcLhArdmz+k_f@&-@zmF?dI&h*VXNN3$F!>rAX7bgY<wZ|thNBj-'
    'v!6DzHmN)&*?5696*aMqY8`CmNa8X*bouwi{di#DD6W|k$n7UZ^X{;<eEXA3o?UQ^7oWVqH=Ew#S~s3>&o`a?zSU'
    'pu@IZn5@99#KSuYy#%$}w{9!6`YkEhK>ev~*XguW}zBex!LRA0J`ib~hg!*QD_KX3<0vE7vJaDW=S@+dO!I0a2Aq'
    '*5xOiZSJ+F<B0q)fAv@Cho_LDp2=T$DAxJ$WPG4@B%})XPaUCR4XXg+hUx)BTk16g+~4uEc`POhJiuo@jL?mGM3;'
    '<=z8eQ%0yDs5%8_&Fe(2gy3e$sJo_8|j!|IFS*GlRuM^vQ!kb+m5XL<9R<RKeGnvHg0(LXw8gue|%4D;@v3*ArdH'
    'v{~yk+5F{?EXhHw+Ht4W|<Ls?DkV=II0c^Y~)EU!jKkcsFx)={4_4`^HOiB<Z!12CaKzLJ@y#=+PQyx;1??ojp8}'
    'GNw$YirI5X_G>itZcU)}pw*N=IfZ&QZlmyqUG!~cHkt3pBRTKm6i`%1nuR4~8gP!}3oB^qOnJN&8P7Ob8Ck|^SYW'
    'D&(!)J)G*2I$-X=)hXMvL^tl^Ms2Llgh?4ps#Jm-'
    'n&%luIkJ{xJVF(}GkiQg}`;LrGMEPPmqr8BOe@Ln@??B3&47}&!pTFmv1HA^mXV-'
    'E}inZo6T%%C`hy$U$Uth~;#ly7(0pUiePwhVlBx;8gau;HFkBY9_H09W&h;-'
    'RZob90lOTu<#N_nv>AziY4O?_WIT=9-;c^W-'
    '0{ze<5J&gjxOpI&s{)`5B?52ug4JZXZPKfSm(lX@SGqycJ+=~&AOdNXi6Rmi51V&qP8ow$#RZXcq4ACFQ+RsnT?K'
    '0_Oym(!FD6;#zyMR98saB7DVyrGJ#chr$1rwy3~U0CWH;c}0jC@8i<npr=54Ryq$okOwu@)%repNJHjAS_Ogz#`o'
    'x@GMV8UBphLPR>W|h;kJ6t%p>{Q*>s2#mpXZETX?LJMeZOJKF5Y;=@AOKiTCh@Jl)?J#>QYSy;_3tF$u5WnWm!YB'
    '_E*z=W5FJMsfIp8U(SnLJJ+o*(N<;R`13=jA^O`687oyoYxaXR@!j?T@c~%qR)^^iqvpG#k?ZCmY)Rcrab~Hi~XX'
    'PNX3o(<$}fT=MmeA?Y1U$*<oUS~h1B<r?mwH@>^c-'
    '#mv}R_BxP#*?(zxR}&umQnSHi{#pLnF@@q($P{y%>1T|A6wLLX{RPU6?IVA!vNl=Ou)J<aJJJLl@;Q~4tK^p!I@V'
    'Xo|rkr9|yk7f_zIfa_U#$&9hWoliQDEp8_--ybS5*P539>j-#LdV9g<Q_Mpp(WmFGiKCk>(-'
    '24S>t@e5*DUr<%A1q}W^KY?xw_mX&gFnpbi#ngMq&J_xaTpJgoWiTVNAP`qE4hW-'
    '4t~J)FxO8l<0Zk@`S8F;d~(Bk?!UI1CvKLbz0Z4)>>dmHcEygSI}D>!isLEK-'
    'H%>%g;0@3Bvl+(MEb{9&|kgv)M=GQN6L1RX~RDHm5@tY%a4)%!9u!aP)aLCou>$uN@`52Cgqx18Utm#oTrN44jLl'
    'QYJ*wp;ZC9v_67FD{qI&N_UMNZ369v7JQN<OW1zWrB4moDqsRSljO`VVSKHR%X3q@#I-'
    'G}4zq9Cr>$tq{G4f=(px#TGz11>cCTs23>y2ZXT0sb-'
    '2?^}q_pPk%))Dr8b|t&?uZ6h;ePlYzrTC6s2K;`k9dCE_;JSN)xZaRho;+_oH@vf(2bP`W_a0p00n_htMZM?z!tg'
    'Ghe}U1TH7YdofFYgHwx)<5gQ#qPJH_4hrjd^VY0|hkw8uG`Zfr}Sp9@yg-'
    'm(<3$k<LU*iCOf9H7L6e44WJBn3??ra>3VXhqFM@>y0zGHq8$R^kR-'
    '4p2e*H#O}2PjDti2lAl?(D`nHMoUX<)ar}3*X$rU*%?o-'
    'x}nU~6WwKgsB)f#r}Yc)VApaa#cjsCh`s1B`6LzvT*8`J4R|%?H4?x5f{(un>lI_c{wNJ$4f?*!J8&Lb{&O`8IJl'
    'ebj4fhHAM04<k~Zd?`-7dzQs!zZ7F?&l3(v0i;qH&-'
    '@L{79xm({fUKDVUxBfiC|1@6Z#k%*oT>4usx%>w|S13(S?6qi?QcoJ~*`IvB4x#q19`wt1GEKLfLGx!t&^3cN>Ke'
    'P8ntrV#8PzQmP?bSTCHB#@JBMhM{V|H|TS&tjOX#KMd1|0aip{Pj&-'
    'PlX`%+IU>s3*nrGdw{w9zj~54X=6VPkqv4Akt6zXAPV9^?q+@k24lY7Ab<OvIMXX*lva3}ap_g0|jT9QNA*t@DTA'
    'B3XiaVb@Sx@es}LI-zx+m|~O;>vP7IO$r{x;?t+I#@0Bdn6jCjH_T-Ux6iX)BO2L)oOf(^1>^OJJ$V1iHhj*>k-'
    'QsIx%>Mlt{}CB`<LwG_WO?V@ZaaTRQ64Namy2a`erAWn*W#QMJUp~5<TkaVMQ0l3?%Q9Bj{J67fsC#poc?d)6|hs'
    'H0Vw|Ex)mfB1dhcF#B!vu51?>-'
    '_NFPiHGUa`4e<^TM>1rmr<h2MVcX7MQNL_QrCs+bf(}I4fIh*(?l&qUepmDV*pJ{Q&^0%gpz(=+&*oGP#tI74R=F'
    'u-gq>Bn1XYIX5wK^6ntftq401DWW07mQ~5Z8Dlgzz*llc(eSu$#zvF4H0+Wm}Wu*;+SjV9WY-'
    ';NqR%N$>CEUwk>2W8Sf^rR$3u<K>4t`;}sd9W_KNId|=g4=rkLTOuX7b%di}{|M6rTEbAHRFJko#V&;u@w+{OFUH'
    'd~d^7e#B0K4Bn`b)^}qX?QKKxTFx|L*l23HJ&_hvO{bjRVYK>R3_18Dl7!(}dJwspzS^hLwX{7HyXqiW{5?Xwdlb'
    ';D#xs=daE{`gFHzI`E0pI`OHq;aw9D-dC4ba_qDc?9RqNs1S0faZ^u&JS-'
    'Z<vk4?}Go@#3p1c4d#oSyvw{><C0n{#@Koj6?kB)ll`=hSKpl_`IkHS&`M46x)KxfcMCj0-'
    'HNWlb!YK!_Ge+&I~U3vkKK{<`=e}9iwdaBdwJ6QM$zrd%j`?;lEjAKXvY#ZN*(;hw`=`zI;bhIJX<Uf*)z$#y_a%'
    '@v4<)xy<A`?i2rzCqIA3ALsqz{i<Y1)1?Q!Q?Q^7z3u2j$x!;RZyZhfGleqBLTE(iJUYE)5mgqfAYJeEBo&oLaY~'
    'tF?XsUjW%B6CtmD)Yc$)frEu|y27sy`!GJP(np~&BLWbo%E4LN_8o*dA^U>RLZd24_kex@i1w1mE5U&O7p!;I^L5'
    '%hWlX8ah3IbSB@nsNw^S<Q!I$Wn|*-vFy0nRqbk2(~wsV{NZ`^j-'
    'WEDto^`Zn+FwreMe{E$kU@8OsLLg|I>D3Cv}|R;I*{u&K$H*s+hzY`SD8%aV}fzjo^Jfo=V{a>8itb9EYDvml21W'
    'Uu2cnps@^<8khxTfr9`yUnY&J>wBKK5;)+P`Ifwz49}l%&y)vslkzgOh%GX^8_;aJC!DF2qo|21r+&p32l9qL?fr'
    'A&<yYG)UP^=jBXsDqL_RNIChdwCKi+Te`QprT|rygs%U=THCmo>gBHczCPV#3@(b^QHYt7FG%<$vt)BR*(;JTZ{o'
    'wq3APSbdBH_PL*gV@CG1sPI<<L3ss9cEJpjBYsx1gkCKW@(}0C&9%w}*|$iGPDBFS_yOjS5?pZ^70-'
    'a$(z#Phx`$BUs3<B=*%Wi>+@hWRtyXnaRPY%z67amJ+JSlkb}G-'
    '<?kUTi^tKEo?TIxVeN|)^6r`zq7exUJ+MMspb{;@A1(|uesHtZ~SqIB<c8S&>=5VT5!siEGD^7;leSb-ad)69tF`'
    'ln{bNwy^w51E~D3c9W{BRlDBdO)d%gRk?x16yy__VG#60Zx)K`Ec8;DrzeFKBs%hfKTI!a(N%i;dka>6$*^JZ0aX'
    'mvki7~~~5KCnKv__S?9qa-IW9j%2Ff<+q`yakAEDlC+>O6d$w**V0l3@{&fuYOukeODBoMW|^l=%p0b3Wp|j0Dpe'
    'tHUh&+A{fd?yPy*G$y$?md(19!p3SHWY4ygvqP`$u(qCWS@5U7%r8%qM?dPr^{R*S$nGip{eyYDxBM!uTe*WTEI-'
    'VX^~!m}jqCjL!$-W=unxZeM>lW%E=S$7btqWhk|qzgr++7g(X{mOr2W{Beq_uf-'
    '*b_)z;`jZN3JC4f9t7IXDeO4kV#T)`)PMh9%=nJP7)ob>2Km$N;`Xj4(`27bv9S2Ebcm8jJ`#_)(sSvdykZc>f_%'
    'iV<firgtTgJtm?2uuatqf+&KjITt=ZLasq<W1Ms?ZHvWB##&`P_@ZGu@4G;IgYRU;5ExL%khIbI-'
    '{sP<lzhkPC0_%Cegk3)5$h4HanCaiyY?;wAc0X(f8+P^>8?9f(rp$c67OwonVl1Wkt>XqfvCNKVF7n{Z6NC6-'
    'omgHezn&Xy*v)6TpX95fFY!ICclgW`ZG7Fm&)j$-'
    'qnsWpBsa{EoGSZ})agO=_w`8nmE}$Qa|5Z$b`BjijTZSWfus{xlgEn`GHBXP8Dn=-'
    'vPTY;zRago{Zo|mp_r~tE2o)DD@b?56?%H{8Wq33K?T=u)8++@RNdGt{@)OvQcR&e(-'
    'ICX*5H2z;Kgevq+}0=)X1^$s+ojT-yr;J2*=(<i&0K%(RzM6q(TouYCtjS-dE$<#ugm!eut&*|B&IM!NNMMSp4sy'
    '?AfL%Y|*Ou?7tIhS@QQiY)$WC*1GOGYtL_E{fd6DDT9>x^Hg(QJja>0eemYrYeIQp-'
    '=)0bRVrU$oWmXRi}{spHQcVcnTPql;XhrzbKf;mbX!J~#xpZ&7~YRGzPpg?uQAjb>r09L!F2jaIJrc{QqN<{C^&o'
    '_@#Cpf7nMPmEBDg!%tQ2_`Z1F3RY<n?OXzE_^Az8!k{ljXlf<Apy5Vt?G&Szhn%pMp`r1NU<{87VUEEDStl)CS7R'
    'O8nqAX+x{zkY%I@k-'
    '3T>YVKJ`0=q0vLWzMC#j(DEqw&GkzX{wAMK^DA(iD?<eSUqzkLfrP<~XeHLXffH_Yb!`j{jv4r+T%uOSeDF@~<(~'
    'R>h_C^DncCVdTmw*TN(B@X-`*NqQ5q#okf1djylBZ}Uar-^#yd^WAe|cEWE7srOa%)<7>eCLMv*9-'
    'n*d$NWnsw+}h9$M++mnX*aEksgo|ZfKljEzI^iY02RqkI*FHWr#x6uZg?Z1`&{>UU3lWekW&Z7$J6C`O|M2ct5(('
    'pGIME1E%y@IaNltb5P&;DDK=-'
    'WU6x9?G>>;p17Xo^_@mPj~hjXkFZU~0G%_ID1$*RdWLAnpT`+tbinJq(?0amdL^LgTA6TpgQ@k(mXk{ZWa}9*sy_'
    '{2HlKe&Gu#Gv$z;EZt);lTGzzj?QyguOTa#`@9Thee49Qy>x|zetgJy(`UA6uMA(*WXQu??D?mOV|o7kV7~8hEDv'
    'f}&;J_i;rdTc@=@<D@!2!(@{}IW`QL$Ed~z<MpjZ_;m~2R@%GMP2W)Rt0xl?tGH%+=9NRokb$kQvD<V(ekuzNMFW'
    'SdC+^LDyAcQ*}QkVCuw=97!VDe{*-L&KJrQ+9d<<xRapx9?pebzV<ZZ*SA9<VK2l-'
    '%L&F56SjJPxN|j1*0Xlcu?Vhw(~A9$sCE#qdhU|(-gGFhagLFK6dU&z-'
    'o{6==_v{Z#j90b}mK4r)wx&{t#c<Kj74N#`K!CSXZz$D+wCG9##3Xz^Dam<>F+vXWu^dsNoD7{ptoA{_#1xTlSOv'
    'H(7<(hFS0nZZ3TC1|NPPb`IZkd?`O2xrJX`k;51ME9SRu*KiN{7Jg{c8@?&#JMVo!iWb;w(iIgmstxK#UJ^s-'
    'xrPV5Tjonof`du#R5<xAiKQ_G%jj0bI&wamN<1cmN^AB~WbPqqR69mNmW8DEq=Zx*&y$v2CCRr|)5jroG}`wjjWo'
    'MUKTbB0*Pj+@)p$giUY6LNWR1rG1MuhhAavRcL#*Ez^b7Do$K*il?mGut%@(3Yb|s7yQjyR1!uj8EG<00RmdCdcT'
    'ik}qfd8;8PnMbG8M4+9_H5XkvCOz}23xo<o?Xq_!qm$0nEA^K?8m)E_U^$umb{R0`lQW=J?zW3t##v>`To4rYd%k'
    '#l*E<prSrlw`P^Lg9GBa9gBzr`@~k@@ykXsM-my~X>46RvBv{h@?e;WZdN}=kHJ(=J`jgVtne_U7B%R!_n4+^*l7'
    ';04y5O{xem~BnEppkkus)B}v`^4&g(A9t=q&Aie1Xnhx=iQ2uTtMV*QssSEh_VAAicVK)Ft(RROMUAv&ITz9BdKg'
    '=>Rn|7u-'
    'JM1|Qq;2wOE7YIlNRB|Q(n#>Rsut;Lrn;lU;cVX~wM{s*t%#lw4OdifUTuKdO_Cso#GrUg^AbYWhB6WOfRF!pW#N'
    '_H-1CtF;0lD)ZG%|?_wVyCuuF|H-cYv&m81Lh9g=-^oXH#L~YtHtqWE$exZ<sM#i?IhQDc8Oa%-'
    'sN6D+IWaW7k5oyr0uLi7XuAxNL?Q)KRJlf-jAf;ncif+H;`-'
    '%=g=$FX!@}#flka?O;)#3NcF;Y`tG=!YKP^}%_sSkBY%p1KP#r)Q_E>cLIw30euXBTyGB>u-'
    'Jnekw`tCjM#_HFOmb2WDd>AE4f3`|`QQH7J#7%=6Nf_4e>AS@d*j&o0N97lMsnu@7_VD~V&xPJD&2*l6OSONqYSk'
    '(u0y|}6~0EFFvm`UIo;{O{utRXw=Zr?#%(GaRTs_5cCTmJyZ5ux`$|~nv3hoW{tNbW&@c9Ks|sIv&Vpa|b>Y<qe7'
    'N7tIsE0mrM%E_3y=5A;Wx^Q`T0FHJot1oj~e!d=gWWRZsVn>y-|ZsR+^Gaw=KCJb)hpQV@Sr<mkvq<lapUK)f&W-'
    '{H$fP(|8^E&rYS~Ju~RSvc2?j)*+hv;watyTtF^)C3KCSCts;bQZ20}X{|b{vAs!!beHz+XrjUwEu_=^kgh*^OwF'
    'HcAo0u|w^uqNP}~QmGUM>qz!(1NLAWayj_@~e*wL7Tq9bW|a9}@Ll21V@`Vv<5x{Fu$USPBLcZ}FB&xQsXv)*?d*'
    'q#gH7}uQ1oK7rZUlymasQAO|VAw@=XIc}xrShKjE@9l@r#7$a?8}e&yK%R8f8NICb4kM_zGG23A2d0iJFh9{)z&x'
    'oa;;V_HM4`W@7;XaZ#jB2S%<i?B{eA6)BCl<C?jM%ITiZR*npX2kq}AlYK!TK?MgDKUQZ*QrqScXO!`~6pOWJ9D7'
    'xu5$zC~4AE%xr<$V`Oas6dFs&|zp&Am=a5x3}(Wdj);x<~Jx-lw&1A5mWU6LQGwkJ0xWv9rt-'
    'nsY{B?R_tp4fe;>_?e)p`H<~QzzUo7aP-fB&CXmTUOR)kA8XLv*@FDL@8G}sFGBuPV{0B;vaHx4>`2rkHveHb+qQ'
    '5O>u}9v>dvQFwM7kEjaIfLwTqQ%%5vBn@v({yylnJX?q?g!<x^w%g~;`M+V0)_lhaAw)9w-'
    '<clZwfIkAoRTkx4nUIkssR3_;X1DfL0hla`zqUk;(#a%doKC1`Pth7)%zF`6F{jr2Ly-'
    '%VgzA2<Jb~`<+$)e7d1N1o|pPEZflE=nkdM8y*xz-'
    'i*sJn_5CR`))6E{d_&uuDjYb1^1&BR+CkavA61*bkGqd0p=wGGBJi4iFO?SX{z6R~0LG?;y#gK=pKVWhGWKlg6NX'
    'XQO;&N_yHD(B&}wjReypP=%=XLuG%vV-Gw*@8e@rX=OgJ{kwIb*VAT!gd4O`E@@FeqO@t&)2iDOJA_CpFdfT=_>q'
    '*p9Nq1+nFzw_2C0<gmTTU1RilHmFK4&;IDZxpIKeauXo?$>tbJX8J}<5Zn-2it80+KUlV#b-'
    'j?2XIMW|ChN72DqCTNPbf_eZ3f3;9_N$4sam!kI*SMJ?^U{g^+e7Q$AEfU2M`?6w0eu->LPPS-QQ65$lowh}K{sk'
    'C_-'
    '#FvRo<c9zD*=`yoJK5ACmF0$24{7Gjg5ih<$&D;Ec$6b53|7;`9`ZEtmoIwe#>kIv&zP)?$(NcJz620Q$!Y(HM3a'
    'PO6O<xaJiWJpO^r&lDKTGhrid4rH1c<Jr7RGnrlg1h(={8VlZjn3*M9WHG~<SpC~~%yAgwJ63D+B;UR~^}z_PbIy'
    '-fx<&FngIDoyr+4rPix2ZHN6UC+|Lgp){Ubhm@p~R_(ak#)<!EPg4|2D(ppB#LXw%uDG&y4&(bFlEn;k+8)${1<o'
    'JCZ<d<7Yqt|!S6X|(y{PU<$@Pbr<bv~$vN3Y>D9?sb;ZQ@aav&Qjb@RW+2MdY$65Z;{I5yEJCTJ<80vPi;pZ(Vi7'
    'gsNSTFj7AN{Vu#`AdOikIBYbe|+*BA-C_4H^V|vIk^iA9d>zqtfo;r;0Ri(Ig?i#Kie}J<-'
    '9oXFV54*H9*q!&iSe>IQyP)99lDr~V^UYOk@a#-'
    'BQ~wm}>ZoS1#~!iv!CkCup$yO8WXP+X?D@oW4{qiW#AOaF<ciaixzy1tzQ^GNKjvJ)wUTdhW3OlY%Ct}X>D_-kdW'
    '#ZSpVOynF1@MN(vh@hx>2j{1R7>HmHJ(sP2(;`QGmk|`qL+g{^f2Yzx~^2R?jSY=X8LMwjZWzHYe$$axvwkl+mTK'
    'i!^Ur6-jAdqgkFe=+e;JRR5`g4hJ_=#?}Yqw7!*;y`IvNwl*rCFa%CXBXJ;lJa!M8j9m|d&=DJsxi8{SG&2d-'
    'bz2ZOcprvrJpr@!3n($Zg=+U_h&TQ6zxNEYcMf`D_w-'
    '}`vfNqG+CVm|Hin(`*ub)$?q_GyN?53WJ^TCPIolKUlWkq2%#&7_b60Zadu_b=xsusD>ctZ761|xxCS~)c=F_|*_'
    'X>Y{xrysde#PJQ`j2<cm7v6rY83z8nCcyEXw9d=)c$uAB_>a#1##1<yLm3n&5faq|CUmB)f(#5-'
    'b_E*chGsCJ@m)@AcZv@p#?uqQTN3&G{yKF#SFMaW?ffkP;e~`SY1y&gYQt(??x)}Z=u>H4{7*}$23Lj84cR=oSNn'
    'j!`0`b;i&G7$<_hLSDl679rGbylz{l9$*Ap<j@oAj;l8*Cs(M$jWqlLM-'
    '@nGguRo#Ms>lpOOqs%BM|Sn|c((P`OqMh_fk{=RvD*_4v(;ZNFwFyvY)SujENdgU{2DD@=xfbePYve}3#afea`VJ'
    'CUCD!{Zs!jM=W(~-'
    'vwTZ;EuZ=L0atQ($A{ei$=AG+p$pF1RQtCljnElDO?j?lvwAEgJef?UyJpb7dlA$#F^)DJT2ALYlIcbK7E0=~lV('
    'TlqmRRKsioi;MV1%R^^j6}aQHlpI#5ZTqid)or;c8o5LtN9U3&7ji7dz7C%3>y6x#0zY2SZF{u5r%gne!(^&N*k^'
    '^=fqJsrn<=VHsMg;=q91wQ0%f=hZ9qP89ZTTq5AK6UWmkD!$D0T-&lMzm<K^~<c-`tz>LWQZ@zc8g#$Ggq-eMw#q'
    '&!AWNBT+RAFddR+deP#!?O7qNc1HM|`j(>I@!xyVf=a)yu@T*Ge_+y)0T;={To?Ux^&$Ye9W1c<XL+B&Fwcsy*;-'
    'N@&VS42JuNS4faUg@?BS`17C(YFlpeJW$(f7*vWabi2wa%+3^TGz|xVe=A0(OyhVmAHieVDeaI6+DCifG>Vvs7(*'
    'k@`{<MPyth=lj>G`{pgWm()PfzwXg1`v;V6-'
    '%4h`pHSP5He!EYkj;xxnEA{LU$gvBHaG;!3g@B8C?3v<Yw)#sD=sT%Lu2A8gw3nKkIlDnF|rN0c3%;=N17eErN<u'
    'c?8g?qbZ7UH0-4N~7#8(zJ+qG4&y>ELVRKjBV6(g1*!h4TOmeXj53}vbP0O74?lLdl5i*M(oEpz_{!8Kf!#-'
    'XZQpj8Nt9X(}BUiZof}c45g+K4&bb5y>6|OR(m+jUfZ#a?PO?Rpq<3oDmr%^}g9Mak!O}QTvXz%saG}dwxZC2Vru'
    'M&4tcR~)0RysoSN1URus%L2Y+H$&{Q9%oWuh4}1*T_$@o(jcn>6Fq)!#@4Liq88V2(Jyp=R9{qgOXH=il&y)UYb-'
    'I8iY!eGK)e|k<mg$L@J~R86~4+COaZ~?>$39eVecMUpSx7{an}e(<y`Z<SNvSs72k|dZ_PgLP&56d=!+~J*)j}(#'
    '69p?&L|PGvf@KTX~7Oo^xd4Q}40=yq_`!i5T|#TsnLBv6S8EZD-'
    'W*lZhG&ajOg|o*%T5=ZNm+8<UUo>S8PYc#$hNejdpGN+<H&4~lrz*iIg+{EtuAB}$$33n*ySM(UWTMwf3HQn$7Z)'
    'xW<(OPm8p=lv_P-'
    ';+fnR4VAp{SJDr`h)y6M?ybG6usAHp*T?%cH7qDy}}MW@leG?Lv1*R8see8IjY^Qv3=V$gxlPP>lRN;bMZmEX#gs'
    'kLXiC<0+YgD;oOJ>6#aXPu84F%<ULj`%tOqt0t7jg;C)RwyuMb!zO)vT>>9BBPZP%PYQ@8l-'
    'AqmTAZs$!W;31`GULY<%<-Wev%PwoZ9ePGHm(R|Grq*Jkymoq`)M_-'
    'N$(T;;>x+pMp5n&E6v^G6!=$pRi6FXkXH@2<I8mJaoM_XUc57nkBF({dAy&mm_33Thl!KMc{#cjsYJ#ahiO}&DRH'
    'GMbl&bBZ4wS4HL(OrJ@bJ?*VK~vzFsP;_)U$$W1!tGhDh^y*l=+f{*`V-'
    '<C9$w%~FTw+2g2pJcTQgmQbH#3%95n5YBMLZUZk|_k4uv1A+Jy5{54Uk+`Z5i@CN#?ACsVXM-'
    '8Iyfp_GcI3lqL=o0FltMkN0>-)3*zH?~kE<II8Pbd?J+0`K-'
    'p3lcG+DW&E_3NU%}Ql1u+`_TvGgbI%rpNni!Y5}#UJ0Yz{dqF%dm+Z>H7NrKiGKCByoP{nJl-'
    '7SLC%Shk5*@vpjpv4IX{nm;W4$<_nJGaJ6R*yti+Vr}m5`-!Lhv{I-%B{_Udh+s7!}{~S3@a-'
    '#gHJ~Y?*1&uh8LUrc~Np4INjqe;F>1`Z$I>+JK7YPLISO~Kjt6)0h!X}?ppxJu}cMVUX;*cpE-'
    '&<i;h&>X=IV1D08?qc8K<xHoY?BQ}kKS{LD@G$P{S6LuC1FNPDyFZ@LW*fFMrnV9`p9A|FfYT7r<KU@t-'
    '(F*dYtZU#P_`|_;tMv4^mWFT&fnEc}Ji5Z8Kwy<u)wyuro``@nZGsgIKV`D|Rq7gBd56u@A`|?2_Rxrqm|H%UY#4'
    '+*k4mDZ9A1)^Xmj*ovRr<HGOwKIOHe68NM;g*?ZzorlK#;R8M*6tYs9P90I8B*g<H@25{I+%Az@t_!vLKcV=-'
    'X!>BCLBov7$U3u)wq6({g-63-'
    'C^!*U+a!@NUIzE#*W!HQHe8WCfE$0cV5Dk*w^PmV>ewaB6u*j#$8TYyj0ax3cw@^2e+2%02GMB|@T`p)>g{;+?s<'
    '!UJRK)Zv$5Ij1F}vPVEphBSm~5ws(lqC4Qg>!ya7!<P0;(&f@>1(_-3HNvO|us-gsk{oM_38+;Lz5YOZYba9=j<b'
    'r`#{E`eQ-$z{^rHOyG7kKrrlW}c$F?WHttwN~KeF$eg4Jwtway)Ez1_2k~6Vch87JAPfMf}7}n;_a7)k<i+yv?4-'
    '=9HO?;CL>Ln_vAFKP_rk+Qyye15=^n<<4Dsfhpg9D({$5rQtA3h&Gn=3ReTC4YYx7q$)jiO1}qZZi4%toBIutEN)'
    '{PIHpc=ZyKQjZ;X0IpU9fNWJ#4=13+2O4p`9LzGwCn!S^qT@V-'
    'j#6GzC6G>@Ixt9#g{dAa}kHV}6&QWordSt5qXsaUC|6G+^(xW{_Je5+mAiMCu6J9(a;n{AJ4OS6^aF_TON(+8!)R'
    '*`KjZk<9B)3Y!#N$Yx12GhxebEa1EV7giAGS8Zkagp-PVWcVR|_}CdPsc@ZFAM@ddo<(uNfo!fdqn?jj{EdrCj3K'
    'SenRI@_3R3J=ri(V(^uxiN_K4k}OWz)lxj{H>mrJIlmp;<dQ4OTY`zd|LKPnn6ggFx8u=0^c=dl$y8M+y@s(UcsQ'
    'WF`Ub&))98ovyN*sW&=$!bR=SKYx5o%<MX?}xjxK{$UY97nXHFsvsQf2Swm;GcI`?U0G<896xkE+4QeLT*<n9!OR'
    '~e|8P#_0?h3^+sIiZbr3q8_p<nAn&s_yB=-'
    'GO1D`so^Q(x=iXv;?g2}A63iO%UbAD_ne6e)ayI@zC)2I?#j5!@p4KJBhi_QPl{f6-U)_&!u6dq|t-'
    '8$})B^Zik2vm7_>uqD*~+(;|KcxS3e(!P^T^U*9XW{ZqezRBbXe^I1;*YYQxiYh9Ue)Ol+tMT(Grr5ZlTNeU+I`9'
    'aDL$gydFCf>Y9u3?%Qe{6js7R*L_gFd<517dbs=I48;Cjz^F0@1W4V)M7oRjMi0?o`~+^lf?+cM1tPyhqg^i!b#}'
    '>+P)b8eToz(Ga-sI|BQ}^7qp+w9#=k0I*IR?jA@x|cv<c&Gx4`CY8$M-q;CPuXi!3?K%)+eLc9W~DW0V_Hd-'
    '#a0r{~PhB$2s3`oIj{)-oyQJ|?O_JkME_yShp9OCuGy&Ex~T=&}Ki+IX2;&h+4O)`jvM(-'
    'dChQ^wtVdwALfB1PFL<l?oM-aOb$%O%w*NynIE#@o^fC3pIl_LMX$UelA!+4Oa8B~4Q4qyrs4NV9Y#oR&?(#MW8Z'
    'Gh7bCtcKWaxC2w;RN-'
    'c;4TEq)@MGp!?Pv}ErPttg>^5>{cw(fX5B6;jfOmWd>hdG7_3|q`tcb^<>bJOlD;;avvti2faHXUG1}96fHntoQ`'
    'BnHHSqlN}2E^1fp<KEZW_#N)W=|(}W$LraQD$tFm<^jU$B7C2^<=wp1DN<wBg_y<WsWvQEG(y)?aukeV!{P@h^RQ'
    '%UMI_MzT3+8bZGL#Bc@z#);0bZ-rQC+lJ_fQ@tYd8JW25j|G`F)rT27FX<ACBGPaY!EG<%7e3lHOu985U7ww!FN-'
    'bX!Xu|wF60fPHDW$!1Lh%oc{XPcL^QXcnVjj$%E`wgzMzn|Q!mJ{7%rid@p9`m8C1457k+yi_egm_;T`_ls7iOzJ'
    'g0@&7Ow7U%Xc&pMzpwFOQzE*TyhB802I9u%z<)wMMwb*~hgvCo?^eLYyBb}lb=W4{h(L#CNT;@9Qf)h46?EcUjxj'
    '6WY01<d+cTAx+icg*hwPm|2-'
    'EF}WhR+ftWl?e1vYoG#8JN)O&G^jo=NelniV|Gb|*jeQim%qwH$JoTl{v(6F#gVmd6g}@oIq<exva}E@Lx}!un^E'
    '&crp;U$Tc(1a;|gk0m*6a;8)XU#eG(pqd{kRQ9uwIxjcTVdc+s_&Ue^apUoI{B+#ZUx+b_S7D3NR(Q0kptI=^CLB'
    '76h!v)&^R>e3%l0U0bAn@@8%`=dfUWXljPDDA=j`VQ{T+oX=ieaHKM8Zrr()aBOjJqd;)TRVXyp}Q_PR1mH?G73-'
    '5R7Xs7GN^BN`=JF#TW~zUg(~{N66;)S9wE|BFmX_BxAr?au7RJZ6fjFPOhs5?g4W$IkAlV~M}|nATq6qijX_;YHG'
    '1^2&PNEVrLC3w>Vk)0zwAxbu`p&$#`gWL_Rt!b_iY@u@nDR@IA={+fj(w0<Kw7pl^`Is?+XZ%ukBuGDMfPiZG&$h'
    '|v*T;7$D#Kd;G5i&?eUk``7$V9B@l0=G#4E#dY;+Xd~<O&@?RI3&u<qWW~^DO4cUc%z8EAWxMh2GD1aZtq@r<eJo'
    'F!>pFSH6JM-59VhaY$xw@h(0MUuS1Seai==&MLs!lw$mlC`ZZqq3)hv3z>#`6rXIuh~yTewYR~nw*!5-'
    'UAR5Mg6R*pWos&&S*e#7do(wYJqe9w5r0yd!OmjlpxDBg<RE+GCBQY?rt!`=8E!FSD<9vc!57Jz@S|_9@X2=`@Fk'
    'iNym@a1pP^LEFUt(@V`U>qZH+iNn#hr#k`f(@J4E}Qn9%IS4iq%clMX!%rjcfGbTvDNVoz05)6H%g@!wDCtr~@xQ'
    'B$D)d=9n;%Y%t+fM&fi)-5^+y;2=m{WHSUyB0VYWdrA-?$$o-'
    'g4H7TaCe_CG^RYodizlPwSNf_sn;mfOTbyJ6eO_>c&NR{Gp#&?&n<*rcnMm6mLp+oH3Hjek#(m5Q${zV;9x7Bo@>'
    'W+(@yN!+KtJptk{x^S6KLUSGMq$54*lMlsS)%W8)rVvt^={Y_3ulGuZT-'
    '<u(d&RW&JYzhMPmkfqGmxoGpjCJUZs?#$OKJ?0xGzUDu_<nqDFMt&{$JD;>lkVZR5QLD>JQc&DQl6Q3Iz&;BKh`K'
    '@57a!7+oN$UWN~XOYA4z6R0}1`^ql+v5(a|nJZ0Vea=O?96ymAE=U)YRd(LI>HLK7t^x~O@38gYxyLqpmQI=+th7'
    'H|hz3+}^P(GM$o0&!?vIDSb*p*koQ6-9}d6!i`Yn=)ZxJH&4Dd<ZWnLeZmA2v=3WrlJ~x0d>$?-'
    'w6BAX2^86!uoGJj`npzH=!G~_SP)(?+q59?7^C@K4Dr85zOBznO$6&&vZKL7;WumKTZ?(lM&_C3G?})jCI^yeIK`'
    '7p~r*VFLA@UZrt)~5Vwg)<lEeexa6S@t}^*ApXDq<7o+FXf5Gc%y6k>3Jg7(C+b+<*8Mo<QrXSTdN75|SG?Lj?Lh'
    'AQgsO0bvyA6S7-NFz`ngMy?#mI<V4d0>ep0sWs-mEzS6<<B{oj(J)w-'
    '@l<(*YlPopGzc9UdzlB2exLkQ@xz{^!^p7mdg{aWGn)jGL5(RTr}m=%0(as~?d#zZfHKm0`-'
    'OO2|E}!H8q^Sl`o#4azNWv2H`#wGJFO)rCV#dT_GBj@@;?#ZpBdux6JaCYl?=b`GSmn=Qpm<XsC>*g43yc?<Ax>u'
    'G%ULK%KhX$#-`S)C`8pW>D)9e7KQ7hmu3oUh-'
    'V#tTI&xnl7ru6=eG*<_26X2}wA_TNGRf*Q17#3^bpvn3;2cdGdEl*F1}Q{|j&68cs_TPJi<L(mTrxiJ#=YekWFbr'
    'x#EWU*O#J^oDC0cmYj%#+uK$!SB>&oIYsWouZmYgi+58(IS%Sh>;%U6TT^<Z1{y9z<a0j#pUpC?2X0-'
    '{ORFI{pjGMnlR6MBXXDPstKkn3ZGFtt!kqU5n+@8c^-k1dpy3m`rGgsbD9JJGzkms0YoPuCXsm-'
    'PrAKzU*3X7<)G*p4pkbXK{Bc*{i@VCZX}0y($voE%}ms{hj5!HERd=w>ioqtjzhiN=F{2|A==kis2%q@43YNdY-x'
    'B8()(&nxuZupw9MXlzx8)6%=TZOW;|O=(tMpxn86{Ig|{_66itS2MTbmrG3`D6xIBj#Jt7;DPnlJejawLUWO*;jX'
    '1Jr7i3M<At`<w5(`hEAoCpJvM%G3`VDM9=?YZ=FSy7*!nf|Hu-'
    'z1f{<V=PDte86p+w+I3I?1raH8}*6sz*E_Ff@cCYGW~wE~q#tI;7_2f_RX$Zc-'
    'M*Xyk?^l!&|_fF&<?}o8(FK)*;vAh;fwy55p#bv%^qs-'
    's3%fE;Aw7;GeJN2{X6T~%1lq=4e&j<dl<*&Bv<)Z~p@(rUe^6%!ZyyxdrzVmfFzo=2bFQAQ2eejFB{T)vh;&bUj+'
    'ZtM9szQl(bt!q}dAh#WnRLmQN+cub>f01r5LieHCpJ-'
    'a(EwS9ay;=Hhqs;*Slqe*eX%Q{(YggE?Ns1&{t!z3oIrf12?`EdVZ&N`+zxhvs)rl=raVB1<YNfG4TAjlaEMhzp<'
    'VS2KAua$xGkyJl9h?@UpW})&c{-ZBHTwA^5rX$xu^!8`s-k0*NEb_X3U(?26>qdd=Tk^S#>vb?RwF)$%X9-'
    'ddLjRf?4q3EB3P}oxQqO!ZrxEvTs&{%<zH$pS^k-H$At6<M?KNGWj4^aW&@O@7Zw&y2tlk595w&-'
    '|^6_GQRa>5C2-vX@ksUdLp-oP`;6pH>*+j3PaLKu%?IrS2`l&PqK4kDEV#%1??@PKR4URX!anbtr!kR4-'
    'woxCJ9Z)B`B6#i|?}AAQiM9l~!61uh2)n&soH`UxZG;6$lR8#FNCkFq-a-JOO{SIX{E$vlq}k9fS4xarjl93^Ug>'
    'eEE|FCGijV{rw~4?i9njrwqD6Re1Wn24lnP@pMTO?%!xZUvwLCVmr`(w+n-'
    'Ndl1&si@I`mHpR@3WmknWIhh2uSTTp$t*l}x^<B))=r=p;D#SORk>pmH%XrPl?Ob%37Pk{I<2NL3a2F>Z{%uqg4+'
    '_uX2bb0Im+_zZ)!72HZodTWmXW8e2}<<r&mod7H6gXx4z#t$g9LX4Q;0|$86V7{)}|^N`@4(ePy8f<oujbKdNLMC'
    '%)y&!@(}k`fS{%Fkc+56<&+M7KQ+S4jTXo=w1NNd>-Zw?0<lU@I8F71PG<m)?FohF-'
    'j`VL;}vdgPJqU;6#Q3`4)ewD@pgG0uCYRxUMoR(aXDfetI(EMiw(LB=qqo+^x3VLq1cYU3Y{1~z8eSAdf>122^+1'
    'w*l5+K?0I)2TYf%;@z{@SLTx=e7~Rhzb&2l`o5-'
    'g;o5!7}t>v;SRk(VvE|0lo#e;lqbJfuSe9w(HeC_mnUU9pbujBu5e-'
    '$ATJ2i{=!c~;)vWrTd>rj`51x>zvgUV+=Bz^U8s`-*c-'
    '^YF=Sk+V0;Xdm2_)A49g0MI|4c%q)arN19?EJY2P2sz7v_J!M)pb#J;xv@&ED_aWix*mskk`9|?LA%y5cb2;)If|'
    'G6%LV4kqEmSi{OYv%yW8&U&5Jqa4-'
    'k!cI9JWcoCGXN)Z%UfrU}kLrz_X_o9usZQTsTq*eqJw<9yP6V;C0aFXxEzNAl>zSxJ|Tou9|j*4aRX&Eeibt(Ji+'
    '{zkn4zfQ-1-MtnRKC1vF;CFj#Ld0c_?*c`eA!4_&L?|t`%57_>RU2*zf!`38#{T!+JD^BVj{_&lcti13S_LdpF)-'
    'OsJY|<ZELtiZmNDXLOqhEWv7zhn__bL+)UHLzR=ZVLhD{(e3qYqdYeVCoxd8d_ABC2@m{<wIgGredZ?5<gCO$@7;'
    '(@6Vj<2DbauzP-'
    'w%+&pWu{xFqXzY$9KnQg#CR3x$()kTb7Dw<t&W9oC{0sk4U9rR2?tFROd>}wyA;rntIfyH6nC+3(g&AL$Fo{G&Xi'
    '&%)f3JKJLXdsXn~>^q8s1zhG;ZCbA&OTo#a0#kx0ivjn~0toyhSAKN*Tdx<aOs}i>HXFm_~XV=c~JLj(P#v2d$y-'
    'yL`Rx5*Ft*YW;QvLk?so^y2=2SXrB14_sTj<Uh4JvOorp}GFv}LwCIT<~r9Xnr>lt&iHt5whzy$;&b_nqt-'
    'M<7a36oo&epe`W`fhX&*$A3G<eLH|p`9~2o#t`#M%!d5m8nY6wLc9MKydyl&$9=G*%^&KALZD|90q1G25P3cxdM0'
    'mSJv$xE_Srb#^8x0D1z7R37*@*VSgBJ5mo2pr`CN}6(<XQ)wV<lC4P`AIh)n9jIP)IVjsArHuJj=+A&?!Hi)K4dz'
    'hnQd6tFV22IlampJnYL-'
    'h6r@zbG(|FV$SbU5D@CJB}RZWnV0LP|Hm&JouP*tGwnNak*S_ej_hFKFB}2kD)xznY8i9a{8OFgLWrq(VvrN$;R('
    '0^&jyf!z&@Qb43D;-uZ!YTWUxyzlZde|E8jW(a>EehTJ!EF)nf`9L8>hU;R!P395rKj=`nM7*W>e;C*?B-'
    'QqWJafvJZKHP)iA77-'
    'TJ%zqR7>WfW5g7OyJKGa5>0=7=^fTb~;yp&h<iW$F5Mq5LNLg5cVRF^*7+wd%;0DA@ZibF}D^^*wV~=4ctX2)NyJ'
    'g7v%swHtxep;8p)7FI8#XyClf`^3WzxN^EH8PGEs_`D`4&_8iB*eviN{7Ruc6Arq7C>nH5)!=zdLWz4CXN}6M4j#'
    'B7W{{J1@Eao39TOrq(-iscGF>I%=;%ZrZx!(`-r3pPguevk$Gg{(_7}Qs~m(0<u2SNXl{p^za13zn?-'
    'w?;!zKlLhcyw-SBETd;Jx3X*0XLiqC&u=X}V#$+oL|Fgp>OD9}A;fA}l_o4mS4`+OXaH>2Unqg7US@Z_B>Pc8AnT'
    'ja^nb=;OgO8c{IA~e~r<zij3sj=@Pc@9v>)@!~2)X!Xn0#u5AL~GFS0|(+x)HUn7s(Z$@N`o@ep*Mc_ZmrTa4?tk'
    'pQ>hm+`5^Q^KS-8A%5-3Odg}Yl=~Ga@ym}7@oFJcKJU#HK0oq47mk0<C5NT)RhkvtF{GF8P9zdcm`u-'
    '|7LnDIO;og0jfV9ZkhQusg>G=A?GaDN$S#`HyV5B;yp%)+T1o24H!`|43{U<}fKA&>bV)2hTl^Z_&s4(Yz5AiJND'
    'E)?=wpn=Swu%&#Edgn5T0=pPOd}U-'
    'TM&2O;2EX`WdF&c>(inF_;`4hdGJK_<A}Gj+I$BI+%+or5{meFtn%iGT1j(Vt#fF{#>rdTQp(Pp%xT5w87A+1BVT'
    '|a9~ysdXsvQx3&-'
    '6zWuN~8^c=dQ(3=qA^Y^Tfr+jkV5XYH^@JyKAGf(Y`o?N5Fn%|G)~~}K>z?BWg`K&?ct5^GB!(Mm47sp(9sgDPg`'
    'chwpvqJUx^FH|V;Yr6zVQ%wUpAqZBzyXE$%7nHgUHF^4K0azKeVYT3O(LMJs<z0hzJ2ReVc?G5wmeWUJfQJ74Tq+'
    'GKw#%;l?~2yfrjJj-Um0Y_h?N>T9qNbV0})Piz|S!Sc)il&%Uz__~+aT=fcmQVF;vmV)f)bSxSB9-k-'
    'VA*QAP!;h9gC!`z?Uss{)W-UsVHoz&o2~~qFm@}~*H-tJ7)7&-GXgv^`_X*{`eJ~yD$3?R^rujUJRfU$Zd)962%-'
    '|q9BPYP!1g7%-r;B(^_y+E)dVp_OrO!XVzQmUq-'
    '{B<&f&9?vcz*xUN8VP_!aaxo<hMTyk$(0pI<K~h#`*0cnJYT<<DWS>O57M~?+2t^5k}=&NtC6ZPf?9^)KvJ1PEGt'
    'vaam*0(lQlG&(DX~&gDpXwh4pUyYcv*1_Wp8LVW#cG`_LK*jQWaT<eIEJ$K+y=7qGckMQqVAh!PsLs~^7MjeZVr)'
    'wh0P2VB%TLxlR=3xH(e56(uVzEjoY@CK}&9!Po?W#j*TLaAZHRH{#RvdZUjvTj6oITPFg%Q0-'
    'bNPh6Uwyc?djO5r$;_np12Yn-VFTlP*rxE`?BEz79zrv@P^>&RuTkWsE}C5B>S^9S(t%G(y2n!z!?=EQ3O~?Q%4^'
    '4T^Jj|~iT|5Or(&h)+HwWjtGS<MNb1ogs|$3^{1*M_d_=tQB~9L#N{hsbiJ3Q(-'
    '?}ffS(V^gHXeQ1)3IXZBAjekg?<r5v_9I4hf#;2CaQ<F5oaKB_yR^RbAX$dGgNN5!}se0)OS9{9Q$Atym}5ti)h?'
    'zd4uu4l2DhLiVZSZxOgZRk5_-hwu&OWT3ZHZ-AY_nsloJ7^_c6~i07Z1Aurj6umv5kn$U%Z#oah;){Cv5KfzO@AK'
    'U#0@W(riiOChQ+t(V|vXB8*WkP(<XAv$%bNPp&RebpHUA#F@n+N8Z^C{|%eDMTdeqc-'
    '#U%|4tf>aGxUNgY!ltz#SrqR6+8Jaa>D|L<6pvzIl<lJ?cWYXQJTkt7)HNB$p$ysDwP)_;3+v&iA@08^@0yd2kVf'
    '{i1%AaIlX}k{iEVe_m`T&-tABEf>1DIu-VfKHQ5EXqD-'
    'b3vE8s~w<L+nnf^hd`25bV>Bz~#}eFzR?b{v3S^W5M(xkI9D4r4M+nRDiDq#b}#Zjti@+P&=s>sqgABYGV_kTwBm'
    'T)Cq>s9hl_Wg_%k{I8)OL57j<oM)hOaw*mYQYWO?w'
)

def ck(q,m):
    if not q: raise ValueError(m)

def energy(rr,uu,M,G):
    def ck(x,m):
        if not x:raise ValueError(m)
    rows={''.join(map(str,k)):(F(l,M['S']),F(h,M['S'])) for k,(l,h) in rr.items()}
    U5=F(uu[1],M['S']);t,a=F(41,57),F(28,57);u=-THRESHOLD
    def r(k):
        k=abs(k)
        return F(1) if not k else a*t**(k-1)+t**k/2+(F(k-1,2)*a*a*t**(k-2) if k>=2 else 0)
    def gram(n):return [[r(i-j) for j in range(n)] for i in range(n)]
    def invdet(M):
        n=len(M);A=[list(row)+[F(i==j) for j in range(n)] for i,row in enumerate(M)];det=F(1)
        for k in range(n):
            v=A[k][k];ck(v>0,'SPD');det*=v;A[k]=[x/v for x in A[k]]
            for i in range(n):
                if i!=k:
                    v=A[i][k];A[i]=[x-v*y for x,y in zip(A[i],A[k])]
        K=[row[n:] for row in A]
        ck(all(sum(M[i][k]*K[k][j] for k in range(n))==int(i==j) for i in range(n) for j in range(n)),'inverse')
        return K,det
    
    def cone(n):
        K,_=invdet(gram(n));s=[(-1)**i for i in range(n)];D=sum(map(sum,K))
        A=[[s[i]*s[j]*K[i][j] for j in range(n)] for i in range(n)];b=[s[i]*sum(K[i]) for i in range(n)]
        best=F(0);count=0
        for sign in (-1,1):
            for k in range(1,n+1):
                for face in combinations(range(n),k):
                    C=[[A[i][j] for j in face] for i in face];J,_=invdet(C);v=[sign*b[i] for i in face];z=[sum(J[i][j]*v[j] for j in range(k)) for i in range(k)];count+=1
                    if min(z)<0:continue
                    q=sum(x*y for x,y in zip(z,v));ck(q==sum(z[i]*C[i][j]*z[j] for i in range(k) for j in range(k)),'face')
                    best=max(best,q)
        ck(0<best<=D,'cone bound');return D,best,count
    
    def alt(n):
        lo=hi=F(0)
        for k in range(0,n+1,2):
            for S in combinations(range(n),k):
                sign=(-1)**sum(S);g=''.join(str(j-i) for i,j in zip(S,S[1:]));g=min(g,g[::-1]);l,h=rows[g] if g else (F(1),F(1));lo+=l if sign>0 else -h;hi+=h if sign>0 else -l
        return lo/2**(n-1),hi/2**(n-1)
    def ep(z,k):return sum((-z)**j/factorial(j) for j in range(k+1))
    def arct(x):
        p=sum((-1)**k*x**(2*k+1)/(2*k+1) for k in range(24));return p,p+x**49/49
    p,q=arct(F(1,5));x,y=arct(F(1,239));pi=(16*p-4*y,16*q-4*x)
    # Outward quantize intermediate face computations to keep rational sizes small.
    Q=10**24
    def down(x):return F(x.numerator*Q//x.denominator,Q)
    def up(x):return -down(-x)
    def iv(x):return (down(x),up(x))
    def add(x,y):return down(x[0]+y[0]),up(x[1]+y[1])
    def neg(x):return -x[1],-x[0]
    def mul(x,y):
        a=[v*w for v in x for w in y];return down(min(a)),up(max(a))
    def scale(x,s):return mul(x,iv(s))
    def div(x,y):ck(y[0]>0,'positive denominator');return mul(x,(down(1/y[1]),up(1/y[0])))
    def sqrt(x):
        ck(x>=0,'sqrt domain');k=isqrt(x.numerator*Q**2//x.denominator);l,h=F(k,Q),F(k+1,Q);ck(l*l<=x<h*h,'sqrt bracket');return l,h
    def atan_small(x):
        ck(0<=x<=F(1,2),'atan domain');s=sum((-1)**j*x**(2*j+1)/(2*j+1) for j in range(20));return down(s),up(s+x**41/41)
    def atan_point(x):
        if x<0:return neg(atan_point(-x))
        if x>1:return add(scale(pi,F(1,2)),neg(atan_point(1/x)))
        y=div(iv(x),add(iv(F(1)),sqrt(1+x*x)));return down(2*atan_small(y[0])[0]),up(2*atan_small(y[1])[1])
    def face_moment(n,D,beta):
        C=gram(n);s=[(-1)**i for i in range(n)];A=alt(n);deriv=iv(F(0));faces=0
        for i in range(n):
            for j in range(n):
                if i==j:continue
                c=C[i][j];keep=[k for k in range(n) if k not in (i,j)]
                def cond(k,l):return C[k][l]-(C[k][i]*C[i][l]+C[k][j]*C[j][l]-c*(C[k][i]*C[j][l]+C[k][j]*C[i][l]))/(1-c*c)
                prob=iv(F(1))
                for k,l in combinations(keep,2):
                    v=cond(k,l);det=cond(k,k)*cond(l,l)-v*v;ck(det>0,'Schur')
                    arg=div(iv(v),sqrt(det));angle=(atan_point(arg[0])[0],atan_point(arg[1])[1]);prob=add(prob,scale(div(angle,pi),2*s[k]*s[l]))
                if len(keep)==4:
                    block = [[
                        M['ONE'] if k == l else M['div'](
                            M['rat'](cond(k,l)),
                            M['root'](M['rat'](cond(k,k)*cond(l,l)))
                        ) for l in keep
                    ] for k in keep]
                    fourth = M['moment'](block)
                    interval = tuple(F(x,M['S']) for x in fourth)
                    prob = add(prob,scale(interval,prod(s[k] for k in keep)))
                prob=scale(prob,F(1,2**len(keep)));ck(0<prob[0]<=prob[1]<1,'orthant')
                deriv=add(deriv,scale(mul(div(sqrt((1-c)/(1+c)),pi),prob),s[i]*s[j]));faces+=1
        out=add(deriv,scale(A,D));ck(0<out[0]<=out[1]<A[1]*n*beta,'face moment');return out,faces
    
    def jensen(scales):
        n=len(scales);s=[F(i,10000) for i in scales];K,det=invdet(gram(n));sg=[(-1)**i for i in range(n-1)]+[1];D=sum(map(sum,K))
        diag=sum(K[i][i]*s[i]**2 for i in range(n));off=sum(K[i][j]*sg[i]*sg[j]*s[i]*s[j] for i in range(n) for j in range(n) if i!=j);E=diag+2*off/(pi[0] if off>=0 else pi[1]);L=sum(sg[i]*s[i]*sum(K[i]) for i in range(n));z=(n-E-D*u*u)/2
        exp_lo=sum(z**j/factorial(j) for j in range(12));cosh_lo=1+u*u*L*L/pi[1];ck(exp_lo>0,'Jensen positivity')
        return down(prod(s)*exp_lo*cosh_lo/(2**(n-1)*sqrt(det)[1]))
    scales=[(4566,3083,3313,3341,3316,3025,4740,6506),(4578,3073,3330,3329,3267,3337,3326,3019,4746,6511),(4580,3071,3335,3319,3283,3324,3277,3330,3330,3018,4747,6512),(4580,3071,3336,3316,3288,3314,3293,3317,3280,3328,3331,3017,4747,6512)]
    tails=[jensen(x) for x in scales]
    k2=(1-r(2))/(1+r(2));klo=sqrt(k2)[0]
    sint=sum((-1)**j*sum(F(comb(j,l),2*l+1)*k2**l for l in range(j+1))*u**(2*j+2)/(2**j*factorial(j)*(2*j+2)) for j in range(6));ck(sint>0,'bivariate lower');sdrop=down(klo/pi[1]*sint)
    D4,b4,n4=cone(4);D5,b5,n5=cone(5);D6,b6,n6=cone(6)
    M4,f4=face_moment(4,D4,b4)
    M5,f5=face_moment(5,D5,b5)
    M6,f6=face_moment(6,D6,b6)

    def upper_delta(n, precision_sum, beta, second_moment):
        total = term = u*u/2
        for k in range(1,8):
            term *= u*u*beta*(n+2*k)/((2*k+1)*(2*k+2))
            total += term
        ratio = u*u*beta*(n+16)/(17*18)
        ck(0<ratio<1,'remainder')
        total += term*ratio/(1-ratio)
        exp_bound=ep(precision_sum*u*u/2,6)
        ck(0<exp_bound<1,'exp upper')
        return up((exp_bound-1)*alt(n)[0]+exp_bound*second_moment[1]*total)

    A3_delta=upper_delta(4,D4,b4,M4)
    e5=ep(D5*u*u/2,5)
    ck(0<e5<1,'exp lower')
    P5=alt(5)[1]
    quartic = u**4*M5[0]**2/(24*P5)
    A4_drop=up(P5*(1-e5)-e5*(u*u*M5[0]/2+quartic))
    A5_delta=upper_delta(6,D6,b6,M6)
    c1=U5-sdrop+A3_delta+A4_drop+A5_delta-sum(tails)
    G0=G['compute']();gain=F(G0['gain_lower'],G0['scale']);cost=c1-gain
    ck(U5<=F(73670242301,2**39) and gain>0 and cost<=F(125977975227,1000000000000),'active energy bound')
    ck(f4+f5+f6==62 and n4+n5+n6==218,'face accounting')
    return {'energy_upper_1e12': -(-cost.numerator*10**12//cost.denominator), 'first_pass_upper_1e12': -(-c1.numerator*10**12//c1.denominator), 'gain_lower_1e12':gain.numerator*10**12//gain.denominator, 'U5_upper_1e12':-(-U5.numerator*10**12//U5.denominator)}

S=1<<50;H=1024
def up(p,q):return -((-p)//q)
def ppi():
    def tan(x):
        a=sum((-1)**i*x**(1+2*i)/(1+2*i) for i in range(20));return a,a+x**41/41
    a,b=tan(F(1,5));c,d=tan(F(1,239));return 16*a-4*d,16*b-4*c
pi0,pi1=ppi()
def sqrt(z):
    v=isqrt(z.numerator*S*S//z.denominator);return v,v+1
phi0=sqrt(1/(2*pi1))[0];phi1=sqrt(1/(2*pi0))[1]
def ex(x):
    # Arg x^2/2 = point K/4096 after 12 squarings. K <= 25/8192.
    u,v=x
    lo=hi=S;n=m=S;prev=(S,S)
    for j in range(1,10):
        n=n*u//(S*j);m=up(m*v,S*j);prev=(lo,hi)
        if j&1:lo-=m;hi-=n
        else:lo+=n;hi+=m
    hi=prev[1]
    for k in range(12):
      lo=lo*lo//S;hi=up(hi*hi,S)
    return lo,hi
P=[(S//2,S//2)]
for k in range(5*H):
    r=(2*k+1)**2*S; den=8*H**2*4096
    l,u=ex((r//den,up(r,den)));a,b=P[-1]
    P.append((a+phi0*l//(S*H),b+up(phi1*u,S*H)))
ERR=up(5*S,24*H**2)
tail=up(S,2500000) # Mills <= phi(5)/5 < 4e-7; checked by endpoint bound below too.
r=(5*H*2)**2*S;den=8*H**2*4096
_,erp=ex((r//den,up(r,den)))
if phi1*erp>5*tail*S:raise ValueError('tail bound')
C=[(max(0,l-ERR),min(S,h+ERR)) for l,h in P]
def cdf(i):
    if i<0:
       l,u=cdf(-i);return S-u,S-l
    if i>=5*H:return S-tail,S
    return C[i]
def Fcdf(x):
    n=x.numerator*H;d=x.denominator
    return cdf(n//d)[0],cdf(up(n,d))[1]
def pm(F0,F1):
    return max(0,Fcdf(F1)[1]-Fcdf(F0)[0])
def bound(T,r0):
    a=F(28,57);r=F(97,114);cv=1-r*r;rho=(r0-r*r)/cv
    z=sqrt(cv*(1+rho)/2)[1];bb=sqrt((1-rho)/(1+rho))[0]
    STEP=64
    weights=[2*pm(F(i,STEP),F(i+1,STEP)) for i in range(5*STEP)]
    # G_l=Φ((T-r*l)/s_A - slope*y) monotone decreasing in l,y.
    ans=0;i=0
    while T+F(i,STEP)<5:
        l=T+F(i,STEP);h=min(F(5),l+F(1,STEP));q=pm(l,h);i+=1
        b=(T-r*l)*H*S/z
        ints=[];top=b.numerator*(S*STEP);den=b.denominator*(S*STEP);dist=H*bb*b.denominator
        for j,w in enumerate(weights):ints.append(w*cdf(up(top-dist*j,den))[1])
        ans+=q*up(sum(ints),S)
    return F(up(ans,S)+3*tail,S)
def initial(rows,M,G):
    # Three-sign conditional face expansion for p'(0), p=Pr(word 0110).
    T=THRESHOLD;a=F(28,57);r=F(97,114);s=F(529,722);cross=2*a*r-a*a
    y=bound(T,s);x=bound(T,cross)
    if not (0<y<F(29,1000) and 0<x<F(13,1000)):raise ValueError('pin triple enclosure')
    M0=M['S'];Q=M['rat'];ad=M['add'];dv=M['div'];rt=M['root'];mu=M['mul'];sc=M['scale'];asin=M['asin'];K=M['K']
    signs=(-1,1,1,-1);R=[[M['cov'](i-j) for j in range(4)] for i in range(4)]
    a0=Q(1);sl=Q(0)
    for i in range(4):
        for j in range(i+1,4):a0=ad(a0,sc(rows[j-i,],signs[i]*signs[j]))
        h=Q(1);ids=[j for j in range(4) if j!=i]
        for j in ids:
            for k in ids:
                if j>=k:continue
                c=Q(R[j][k]-R[j][i]*R[k][i]);d=rt(Q((1-R[j][i]**2)*(1-R[k][i]**2)))
                h=ad(h,sc(mu(K,asin(dv(c,d))),signs[j]*signs[k]))
        sl=ad(sl,sc(h,-signs[i]))
    a0=ad(a0,rows[(1,1,1)]);penalty=F(a0[1],16*M0)
    if not (sl[0]>0 and penalty<F(121,10000)):raise ValueError('0110 log-concave penalty')
    H=G['S'];half=G['lower'];extra=F(half(('11','11',''),[2500,4800,6200,4800,6200],500,-1),H)
    tri=F(half(('0','0','0'),[1236,3201,3201,2560],854),H)
    if not (extra>F(138,10000) and tri>F(14,10000)):raise ValueError('forced five-pin and overlap')
    out=F(Fcdf(-T)[0],S)-y-2*penalty-4*x+4*tri+extra
    if out<=F(1,2):raise ValueError('initial strict majority')
    def ceil(x):return up(x.numerator*10**12,x.denominator)
    return {'initial_lower_1e12':-ceil(-out),'Y_upper_1e12':ceil(y),'cross_pair_upper_1e12':ceil(x),
            '0110_upper_1e12':ceil(penalty),'overlap_lower_1e12':-ceil(-tri),'11011_putback_lower_1e12':-ceil(-extra)}
