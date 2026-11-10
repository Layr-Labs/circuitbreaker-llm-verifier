from fractions import Fraction as F
from math import lcm
from pathlib import Path
from hashlib import sha256
import json,importlib.util,base64,zlib,lzma,struct,time,itertools
path=Path('/snapshot/files/family.py')
if not path.is_file():path=Path(__file__).with_name('family.py')
assert sha256(path.read_bytes()).hexdigest()=='f3cadee24242351389e7748a8b8346388c9416dcd43b5b0f1a70fc12a3f4fde0'
spec=importlib.util.spec_from_file_location('candidate_family',path)
f=importlib.util.module_from_spec(spec)
spec.loader.exec_module(f)
PACKED='{Wp48S^xk9=GL@E0stWa761SMbT8$j;TJBE^IZU8zC0LGDvkomu2;)R`WlZ3Y}ZJ1FRgv7o5O`R_(tE^-1p%^xR=fla}noF7~FOOK{+kD@Wz7W?=C+}71fR57qI}IBx^bTf`2iydyLTR$6FUKpg9>Bj3$tTSW5^rEWUm7d-KCG7q>L!zij*Ua$z@cgorjMk>PgaT&I^IEi8mauU^#`vUyi|&ilvYU=7)*2k|*$N+zPzuAwDLh>ZQigVH!iQPf6av5AdT+$~6M8OWZ~T-(C<oDy%`JZvcl<NpYq+Y4`$aT4%yeVJ4^9m|W=;e1lr^#a(7+%8x_4QYHfV9PQg*C|7Vx%;(zk-9oAUZztnu1+`<qX`F#o12%eGQ$*rT@wx|02|Wjq2$#!+R685k@Nh+(zTrEV%>Zb+04*6dX1-ZGmEL?b{Mg*b?%U&anQqep!0s7PW|7M7oy}ItV^CFdvnr0uUqmFkppY2Jy?!8(E=~LjgwA8D>>_zW2GXXHt?;;%>Ru$V4le}z<OS9b=rxRPXYi=3bNQbFT=>r#`g$YwPFO7u2unX)k1nQESPqB7}l<z(?jaYhmLI7qx{0ileqZs9*kq(NYH#gaEWdQ>Iar@tOb4d?3wPTVOb-IcsB$p-`KKdL@^zCpi1ErQV1E+7JWY%g`qrtkd_gvE<6Y-_W&Q307V)RjuJ8K*~ywn;c0Z(>M-SxBpkeywKp2UdF7pNInK1%Kcpheen@=4Co(ByG4B~+Ed3UnrY0#fTg9RsXOZ(1zu!%hcudu>WB;D+oVPj<W2Mua_8wabyw@SFtjZvycEg>}PFuwEo06PXgv*y6^fhFrRf!hhl+Qn6@&0a3=iU3I^`fJAFf#=FF)l~JSEa%ue{7e57E<D^B*D3T1<->C(gRRi3%u~W8`%*K^v@QwlzWqUJjxALVnF2M%m3p_=pYiD6+K{fao|mQW~VLq_Nk!!f94(ZU-j;n7@72r{iVK+M>+<g=eF$0<LH5ta|rbwqG3<+yS`AiUD2I?zQ~@UbA7V_=NrK(oOpi7SOiYb7H50m505DUa;>AoyN%w8ejYpX<GG?!f-<1IqVkK}0TRf5tz#U=LXSJcVi~?(usx(sXW+s~&PTh5Gkmb|N<~>OE`|Sxty0~c%Z0l%rE)jvrC5=q@2n_7zf0|OdmFcI_r#r?oBaFf#;85z$9MM}82P7QKYnD7rVhzN>H2!WIXh?U<iDuRP!Olb*t0+3)`jkJ0VU|&AJA+X^m7o{go?2%Mk(71Haae47Pz8P$?pTx)B60Rol?#Ii4l)1TF7W~xqmj6!h1EwBbI9oJlawG`2ksv472u81EAq0KoqK}evEO=!Rlt38eaOArbSSYZy-Iw#`zgMhrZ0yn$}&7B4Euc#hTgor*)=?E>yf3MG-_Zdt$DUzMCS(!D6Et=CDA}_Lr$NFNx|GGm?~R!E-!GBn;UVmxm<8CeZbqm(O{~EfVxW*@9Dzls~phohz@|{nBirjo{H{`ylnp2M-UR2Gv!+!}_{gpp~t0tNo@i|0RQO+ja~eVR%z+<TDqC2$_L^oMCk7&!++V`#n3jjUw$bjwzfWv{hq7-6=w5%liApS6Y}jWJl-_z6n<^n4w7p=@GlCLWgUZj~=exl`NnhInI|CTsG1uWCT`@#D-*CSC*JHa0Ohr9zH*5%i5>^#3?jdHQL_zBXz(K%oSA?eCKpWDOm~eAta*XE!W?dBH!9epuBlj?O2_cdagG6sM2R)N2sjG_zahR7${xX&RE%PF+zDleMrh>OQHm3VUM5!#(k@H2Zys6;Y)k1zd5$c^+8GB#=Rr{PZ7vFQ0S9Y&aRFLL`qxlNU}CLgpB$k0GIH!F_!i<gT}&^vEqs4@|PQNq<Px>g;9IEbP#xq)*O;8CeAwn1Voi$@UN=>OHK>{)EB@K;af0qU@x<#G_rh>E$GaPo7wn{fG_b;4h?)cm1LMm`=k2Oi+_&Yp*!0jK&oEym;6t>wjW@_2(}4XWki0~u)zs_HfuAU!Ns)7%@zb#8JZ&i>FCg`vQBwhAO%GA>D#GFPFx+`k0?xq6k*w(Sm3Ld$vW@m#B>b$2W)q!1!=g(ODzE<ctHY*$zsN!dgJE&cXj585a&g-vzBJ0ie2v(y2%;41tsmn)&#t5vt?ozc$?kWQc2l0DRY4Dm?e9XX=3#&ewL5^|C+6+^`KNVx?3bCbE^G5hqzMz2FtOd78&VJto~=3ePv=}ZMe&#oe%UPE_Ni{G|L7qmR0T<x4t6ky{Fuc$U!pcpgbG%t_Xv+-5}VsmZA?Nrm^SBLoH2D=gsOH2E$A7*NfJ+tln!=eRhI=)spv@=9xazS|dU9f{3X@_dd(T3N-g4SK)L81T~AxT>0hnf<B1+;CGteSL0E&NAs2ppzM39Sofjkj7eUy@hKT6pI+cs7CR(F2v4VMDEbSQYv39YyQHV97Bdz?pvx|leSpw?migu=1b(l{oEYvqVdoWrZIgXQCIIFnVU?{r0xs<2fKkcPBO3pS$$CvUkFG<gYDMZiv|M2CI|?H~`!6hnTS}n{rgoE~UpzlWm56W0CXpT+-}s*jmrR}2u@j%81LhFLTgf6|<WlGp1lOHp(whyA2asiH17Y-F9pY+aseyqC>U=Cn7E%YhK|#zV=JVDiF<H@C#m#Zf;*0|d3r<cf%#+7$okqaG%-$<Q_O^K6|7(D;l&Ov+mL;fRTBmaiA1-gYDvdva?e2fYj{Unhi?ilnf!6I)A;Am^hR_<kqcuxN;SxW2zIB?=#V)OR>Wv_da+en7@sy1E^wX*T1GUf9c#C*jm1kxF`egg<9Fq}XIfI2SO$DdE(Nntn19(ZB^C=!G&3=qs?=ND7j1+b~R`Sz=d&Ac<)_{HlJPTgs+Fngu#!|Jm%|)hv%1j!7!g~;9My%F<VhV#&{99AN7JA?ep5BJK?7)FUdp!GUsw{<<)L6_q(b}afyuM@NIA#ed;7dz9Y}!u;s#H^g{fidlQ{Ddjyn=VNGBeoz$Do$oPM+eOPcZ9a8GC>h6rHdZn0scuC7$y7wQ9V&bzUdj2N9I@cC~t+o{Nm1W%X5z7&LdJh7gc2`Japaul3d)g}-gzwmJAcPCK|PsBO8%0p6mu6S$#5ND-c?eovGdWcuJC@t*F@3@-G;3kZmCmxEaOkA3+5XROIC7q-TYCE_JL&b1=fo^#$rwoh0he&nPNC0m2UQL!SE!=7OUs+p#)EWTF2gw&!Vebi#R1ZI;iaqRp$b#dHLX1XY<OTtCal5EVtg(H4;=&5vXEdMVgU4*Ectqt8P<>J;~m8BO7mVJF9TqOgbxEwfD%gIZHtmw~*|EetBOgkxZ9$@0j28(xip@Id9s!pjW^wb;6iVdxkLBX6Ip0uj9eEmBf(N72eLxEx72Ujwp7=(;v>!_fHYCb2Ld&=0)b1($Wy><StC`X~FbgbF;HtK3F@bxB6T|C8+QbA}bKOA!;l61(Xun!I><Mowvx&K@oZtTx#7&fGrlB%~#CEKWr6|Wu1O-=W@m-n&?zkctNvfDPuMGHGImu%I)%to>?=T>!ZU?zn|7uf<-uY2Pxfp0mRJ1JMSNzy9`_hd4TeHiG#-x02w^%^TQBHw|%sEuYHF1Rype5!vAX?QOseg$|--48`Da#+|KjF(k9E|qdw<*QGmgS#Y;cq;}&q-=`>Xsp7|(3L2n2z9xW6h`B_49CkIG#J1Mcf?QWQbFpx3<`eDc%rr<0yo|R()SR_Rl9IMMZ2gxwm|P9jsHBBw;2nZS@Kx-<3+M4gyUrDM915U_1Ubgy4F3Xu3CqUcCRbl=#m91(-CYgNW9KWh^x@EEs!_+gcYCv>oQm<@EAJCZ_ghz=ZcXjH;AxR{JrFB6>qD2a$3d;J@wq}j@uGa8Ef-ay*a*CY9hc9&oOs9bG~3WfqlJ@oMILQy{c|g%8(@h4JKO?J2a2~LY}AdRO`D&*_Cn1OyzcR>u9T`4*oA0+1sasijOz5ukTNu1;NpveS4bs=jVI+gWE-XCvauL>G8?vKZ)2U^KkMqT~($Y`ja)%DX~K?k`fclPk*Cg4~g<pt^Z7f71^jAVpmt%JV<|CuDIh$5+wjTG}EMT_Sfod2+KAiOK}Jh$N_nXF+tD0KS)FcLe;nSde*gJDr4LNHyk4TCv&I1W$3ViKPRybK`Ogb{pI>(+YuH%e>JoWcM~^Gw{B4N`hHVqRKSFamnb5lmk1!9g+}~>Wb@ct@%OS7B2)b`dZ;F<Lg32GcWOv9(dd}CjtLvCXX6S-Fy@SillR{>Eu*w1@-Kv0@@x9nKk+U@^+Cm%s2G8Oiux!uULVS#x}HVHmWIFQ`dH@(fOLB?vy&Fe#ZK98mV#dfpeVA~OZVwa_ujL6u!KQH=rwK&dHOb|ndn#4caH*6VNT&mGfU+QK?n|oM7!1g=*2#qx}F5swwcmqYZ~dv+71s)EA<C+7wAXx_%G{K2|F?4Iy52RB#>jpa^cR>FMm`a0}hiSYn8^dZQv|uzSV_nA+t#{Brlk1;QP&y@lA1?I3dI`R=~xU)zP&K5;x{qF4(z3P-QXr!hm<nV&5TJYt1n4+2NYyTD?<51MfV@mi@(vKvvr9_(0yT@mVU`4CwMT-P2eeaH`{()w~&04odg`o7#Yl$^pr3a5sgR<3V4znRhtOiC>C_oY)Av0AbD|iNg2x_XFralUOcp9hEu{lbbc;NvRNEs2|zZIA&6gawBHd_VF5_(9`ty$M#bdUh7*6g)c+NeCEeqTH(oLU|##x07rFMN$`vBr7B!Yf%JZ`<L5xCFM{cVw+w+5%dpmjCiw1KmN<<^kX8FDmIA~9CVY*P`fJ|#UAC@tX8T1RWJ}=kAN*0$392d1&D}1H(b&#?l_?Abueu+EX7bLui#=o&^QQd+%}bzG1URdA)Z8FNp(EdvPh(@m^Q<S%l5pKqbI@NJ0KdORAV%x;LF$;Hi8I6&fxtM@NLFV#f=em^Qw-*-B*;^4ibL_?;s=C#%EJ?m=}EZady1VyZ2g44A1K~>VA;Qc66w?7knmZp5M4YnQv#=QMitepDrP#)+?KG5GxIsJrKu?`lfX*)8bEppR{d&n;UwWrs-Y5;uUzfYf(OT9Xql4CH0y5ng-H*!hCu0!7nP4aRPR15hv50?DN+n(x>;V%0|6OtTaYi-U>NPe{;kTiskovKI+(w`>MOD*yQ+2^$g!t~45mEY(bXC&v7d@R!Crsk$piqWYuY_uZ<HYrMy@*1e|Jw1>J^xipfn(8#`_vzyA?^)&&e!Hg{~@Wrfr6%S<lEN(z<`qpg8f@`WAI>3b3M%F5vwWfV?C|s2VILFc@t0t_W!18<6wG@CEd7Do43fV;FBg&pRdRjgsk@x-BMp@giO@S`=d=vL|aB+_+t*?y!LWiNbzAL!Wbto1zhy)C^{#*Wftg7W&UvOcQ9D@JVWDmiR@Re#qQ0%;G96l{wEYZ5Bowd}_>QG(!q58wQIh@PH&{F@y}Gre5-BwK$h0VD^8A-sm5K_{;g&VXewkC;hlDUbN9n;^dv^rAB?^I<!cKrY02Cy|T9-2%w2zWrHc4kJz;6Yof&fBK?Kp)mSvrynD6GVZK3CLO=qkyyvxT%vqNLCYR2t{J0MXeYz*Qk6~@+3fgW`HxmL)8fmS5Mj-BLz*$0UssGb9s#v%~SZKCh`d%|toKkL!_tdp#aK0#8qYzQ|$4FK~sPYX-A-BnxE!ur8%X}^Eth47SqT{mwmw3pK--l83=ToI&lxjDBo+}G87vT^)J^S`|2bs!j3%+}NXij$?%ktKI!wgl86%=K$>sW`LuCuhwnp9gK**Z7pgOujjT&dhiBWfbA-Zl|kCJM<gyIR;YJ<VuILuHlEo^tg{{;~u-Z;NfdoFcMS;~?mEjkSo&^X&lS8nT=Z)wfPMgD}Fv?*KW3Vj9AXIfo2aqGy{jbIX9`cc6U@{6+5@Jx`WaWd3o&Ap34k(ujqY6X?_f#kkR@AdqTK=A>GR$X9@dA(wz2=WzO+Q6VdF7$Pxa3?{Vsb6SU`Q?npy`w1?FuQaoK6<uszJAC&95M~*^E4O=G;+IdBhmVW3sRECoP35D1^3P+^6fBk%6VmjcfVQvpt1pEMP9J?a;3i3A{k>C)>kd<!o~bc={oWxL=NszPk3@jiMS{!Q(yt1fX?(+he#jOS)A}e?CErIrf{**=RAlfJm79BP2j#{HC!sg5y47#v0xjqem=_88kWV}4iT*Y*&C1NUHn)6B5*r{ncSi0T_GeTOo~OGoeJqoyXRkuB?6U+!h5YV2ZvAG;ArXt=?}w6H=7p>;vbq|1t0BCJLja(`ll%egd(FAf`LA8G7<uxnVb~?$MJm@WgzsLJlu53hB=uRap_qL$%N_jWz8Dowz}4AnT&v|`JB{Uj0sv=ks-Y$9C@+)WdOWRMc6yRA2w7i34Ulf%5H^ITfojxlw)@9V-QD_o8+p3AF1%&wP8ew)9VF`SU-&%8sR|bc<s@>syXAx{D|JDhvC)H~#Yaw7ll8yu?nDnomFK{XD{Z(VgI}E4px=T)r7^$IPbmlJUQNSmBAiPw;2Jg|>qL>2qsdUC`xC7cPBuDA-MV4p&ZvTPR%?xEf=g}ECdDsbLo+K$xcYvekP7L$+3vCM_OrD<84wx*Xjn#-8ZD)=Pvr;S)}iMder5m>^+sIW{jvMEl%l&I4Hosr>VAe5+0`x^@med;I1IB;<E5(tTG7v}SRoOS)N1;JFgxo1^7XM5&I3f0k2o1V;TP+;{j_(gJiJROukcpZV$~snpw!WJk>eUla=z4DHx}Y}{-kI=X|fcDqU=AKfDW9(2$XU%=orGV*tw4p_w`!)b)Zx>N85mF_%vz<2N1*GCrR<5Tf7N3LA))HnPMGKuiLd`VF0kj@~Zk7+uKI{?S~pl36#yVh#ZmiOG52G8D$-dcu}7BT};xC#7L^!izsc&C>XOdn?In5UNlmNmwjJm&plByd!Q4w#4abpJ$Sz$l&>^sLxgtpVkqOY$gK8jd4}1Ywxow06jY7)IF?_;NRwqDI2%RrK5%P08S2Wd5}Hff%`rbyT6Ohh>1Wzh>LPz=IT;pca)DNWEJTo+VF<3HG_U@Kk+;ImBDs+4)8f2yc`0xgA5IErVoGQ1TawWHC8M>Xqpqst)yUiyPpdiA!|qOCQFkqJ#Fgo3rY)XQ>@==wEjCKln_pMgKh7)JDM`R7ik8P`3$84O(x)SjelW{Ur#u{ySSO(U?<O(W`TSD!Pq4Hf%UW>z^oSfMff>^sd&&Qr^L-RD=UrA=NtI9w60Icf@~ofJ+8)J{?3#ecB&J~1qb~_pIm~=8%D&uU@w?Az6F-x8ePK3s*GK|bWzQm+*=iW(MU=Rh)rCZ!<)ayEZtS%l_rT}ZEvWm}KoZ7Gt$V@<#H+uI90^74vzOc<2Z87+hA@+R3X2Sa07=jVdwL2|Qw`jT8@oW(faQp#hiRQh7yr;)3K+e3JdyoyOUB0I8o|!*Gb9p=H6J5#pmYev0zCblE9{)z@|~(+Fg{Bep507<dy3kkUd6eRyvgB3lcesZ(J&&oJ%l*tyA4L9uCGtKhyj&XISev6x*-#_S7>&*kUG{4O({$r>f5A2s|}y3j90&hnQ!Kh9aBBEFZUmVyhQApD-G+95L+de_C1r)n!imEe3xlsRR(lcm&w!l5CMTMJS0tJP#zBu!TX>Cr=X_Z!nY{3zvqDuvqoNKh0iky%e~IG2U%Nym;}8ZvH*$Ymrq}h!;-v~A2rnCGJ_oDaCt(#9Pn3gz0~w&iz!Ms&(~gBe_AB=T75#m|ISovlIyqk){#^QdEa6n0g1pQCn^t~x3d<qujiWZcp9jDgeU<W@mG3M?x1dtoML8i2stcVVh8MI&+<Hr=|c;GwibVoU^hwLgV-6*&Rw;sIN5ZPKZ~f$n+rBMWK1O3R$rdOPd>i|Ek|;6eAI2dR?WdJBa1)dE}4f$$qcNHA%$M$Kn4WoUk(fi$`5APVyF<(&`w~29p~9D{*nnfvh?-n5Cl^FX8cq^`KgUG&5L2Dr$+99hyO4aB#YUs?*=rr;<5*5JsLy3nb-<nw1Ky*WgDHiOBvnqUetC;GkZ<SpDW6_T*@UTtFEa;%E5c4N25Ub`6QxPp7{&vQGGNE1Y(J4$FGc`T~8V#o=}aynV;dwUP?e;EPn5AL@6i9jFNZD3Y}_dfacwP0$25s=UGlBz^%iJRuimZa~%IAXdohZ9)P@CdEhv0nY*jfD}OM#;*kA<*hHV4g2?e85*DgZLsp++uf<T&*)JfMx8b{tw1ap|c9$R6s#|^L^T8(%yI?NgoPRsQtJ+n%Lh8>fBzH`lRlhhylkj~SS+m|-Rd~?*kwk3mtQ`@`+@YVbs%QV4T3Kg{sCSNKe=9+X+}lKYx27Hg4o+MlUV<!TT~33gO>%K`mAlDP`3+0Os`TMxiTyHV-p{odz+1)zWkf%5|I3-Iv0DC@!|S=&Z&G;yU>sk5&^_hCXLYt!r@s_59Az}^khkB4uv8BPe4BI&M#Q$w>)(eoTk?8B^FKNgX>9f!c;8s6+jcA)JCAp(Hf#oS+0krfPeJ>C@0@W09Qx+CYI#dV*`kGmrN)yJ^MK`Bou4B6R<OZ(u#jO;JWFn^J2m+Vyf9<Ennm^TO=X^P4@uyEsE$GmeNwe|Sgpd=j1OmF$N)eO=Wz$|+3~JbRE=u)Poa?jL*kSpIQ^++vILlIqQ@+bALuhFH95=}-UIab7j|emXbxG}IK1l|qvZ8XXlC=^tZQlFY$Wph_tu5icKcsQ+C7e0pHMhAlK#<Bgcc?@Ky<as8ak6~XB7<wjnK{%`#1z2yp>t9P{i-0hdurwm6ARdt=JdbFMeOPGAp<<WfJHbCyVbaWUwcIxmg*`R;khV?_>q|>zBfwfA{Jf?MzHofnTPMNjpndCu&c{=2fBx9logS(>*j_xnm-2JAIc=o27Tf-yF?+fGycKK7;D4CkFBXsaAiAWhSS=%eby#TpKU+XjZ=L7~9`xVe4}#3;O`0hU`38g(7BE!o_ft$Q|ly8zS@!gk_X7n;IYu)!`l{JHR%~^kNwjT1LvYJ0-_Rx73#K3)F}X`R_3o470TT@Dsqq3Qv#l9lNYpwA$~uB<?X_p*pxd_0{R%pQA1p#msi&EZ<NsnsSLT<_HkHskED(U=rEIWpw4Jak$gE8trG49O6}KuJ*M1cJymKbrR){DxNDq#JL=tY)A%xWur$ZPVk3e;hqp(whp@=Xk|3!;v*ibTO)CK@6JU!#>Biw>~2+AfA#wCralRwB_5N|FXyLZiM4n*7)RTg>&c+0USAy1O>MLjBsn!;sjn6YBUEoq+8>qo6VcU;uFYLQ>EGPhe?sP<6_G0q7vX0ii*orV{4&7YL-bl7gKC}P<{hr0A+`a70^39+zaQ?pyX-9(RF_LnD;X(m?;Qr2r=`gBT5-1U0kMG1TW5K{%^XrT%H>jj^SPI(a+}wDxMq&xPj1dtY{!|mEyg0O&K@bHQNZkyjnX4TM0k)~%xnrPFQe&Ve$mm)lbv;T(t7h5Ob8NrqENMfa7tZT@ndgc2Y%nu-=*V#c_q;1rs@@flZ^7qblQrCfnn(h9V=iJUSLmx4EeD(@?i&F%a*V=<+c+XXTmBJIXC};7U-EIQW-GH0|U=#TZ)~zjZ#i`=dCiBPQZ9hfs!zj#(&yit+u?{be%D01>eVpH>YKz?tvhxq*Lar(ciJuiP)U>TF-p`j(FIy3S7;YS3$i2BbZ@aOVnutHA6DRBp*w8D$ek@l>&DH6}lqy+fpaf^Zn6h31hEozV^5x97z!Ha7{Swr(OL0)?0Dd9j3NoQ?d`2CbMrykD--A&F?&YyQxm680zztrpFkQJ#kJ-Sp(Q9#}l2*8L&sQZJ*CyRf~CuZVr;|jEMPc(_G19#nc!N$#~PuCNkk6^!QF@_J66Jr%5O*&0QqVAAj0yp3%XM1Lb)w_$LZ?irw;Fa2j)5*c@~Q2R>G^Z+rP<MGN@o^CPEotq;e9P1h>QF2R>DtX&u%9>oQvgLJ6!yegP6jepuZn8&;n=Xj4+!VJr`D?@>V&V)878PEW1S>7vr!OGBVsS6<p5agd875{t~)c^4@2ISAJT19ee6$aL%XfPWwO7Xn;<_+YE4oyAZyomYTPL^<YZowX#%7>2(K&4toB>%sn2sbXR0Z{AD;3%Uh1fo`y=*IpfS_tF8z`;YI`P7%c53CzaRhaUAuljMcS4a^X5l|=|yp_aHcu2FI%0E(=WQRdPlcxM}F0+aTj1Wo&07HR0!<eT3a<{*nvm4uEUMNrMCUdfoPa!SIEBcci3hQMz)6BcmV7voZjNzO<ZQX+4FbC*t!nvhM^&4S<Ff0aNfoSMh5#XeTG|3i;uAZH@yP)-*AIA3pn#~KVzo5Xzu@L@SjHn5E1ks}qCpw+lUc^>OOQH(`DQ9Lg58IKZjbUgrgz9m^1?1#YH7`Qno#OUvd)>keg%g!K4%}=sA<rJHWHRO{2_8(?=^LTjncnHm=Cim+fkRwPmq`e<*#7vJmgn{3F}fUdxD>iGoR43xd91Y_*^t^R<0(CX&@rPuEVB&;LBZf8Ha3A$q~p<^L%pFEKaI*zqW7*$WV|U1FOU_O2>ETx%6ZJ;sRFu=)scaVlQWx-ouloOr|AMpvZKeLucChJ!sTiu(zc@OQx8$aWdEtrt66!g<Uc6W>m81b<&gr?*UyC2&EYvr_CzfCLhasQl<8BgW_&n8#YfxCwwSVQRACyDE~M^7i5bAww|XU!r2n!kVW33(N0y#Z7wo6{P*uJT(mWU&_{&`s>dWw|(Eym=hSU-nwR|!HWT0%?1(kTlJl!2V?Ml~Er>a{?sKxVB68G|2k;K_^SlMXxvBzMhZP4H6eLQx?JLgx#I2LMW;Wc@c4|K_>-_Uid?pY3r9F!Urj54-pEtE7!WEjb~+>>+T(?TAw9@;|@6}Gm3k;!0kpBFDaHca1A?}<S9$V=^9SAAzCxduUgIO8Ku<EFm*!xPaubZ-$XlsCfZC%s)Au+ptVp(k$C4kA50(<Hwa@OVHhZWvepjGY_WHixOpxYSQv`6-`hRVmlo7?_PBFWI1gs}nKDLjhbikd&f;ZE3kyALQ*j@^oTuJQt@Xc)shQUT%YUI_ku7T-fLV9NF+tFTLlo8s#5Hky&2fowAMqqspP?mX+I{UQRv5l9kO^<sdfA=t2!}ly*J_`b5=^mFY~ltK+v+|M<+RnQ$QGhJJGlnCc)a3&{6yBTl@J-vcmAYolo(xGezd-IEw@i<7pMCZ}?NOTcR)>QjZu(Gx)_BVk{?8LC`~T+gcp!mPB=EV&kRC%|N*r(^PX{6^tPw<j2~Tku+o<3sHql7g&E5Vf2<tx^qHNJGMOJl>(Kub84(F_Fe1UG<z(liGoevH?Q*PJ|kVH*Ex<@cW?|j0{$OKHqMjXl8{*n7Qb8OUedBz*YjPpc5@}m8^qRkBm<Jz2ce!Iqf+&uVZD<4emaCK1zBOyx05~Bk~T}<Wux^9WPEl9(YcZSqNf;@ASzgHU&Ww9=K)zNp<6Wz@zMxi&W0_D~Y9XH6~=PkLj^1J;e)SMZW-dAg!wF`nQ0W!l7yF5_65I$Dj7r>eOcDAAwL$-gqXJK&nM7l98|pd$w)RoQT<=31<k*Op{yp8iXN@A+7#8!^2(J9BHIhPRxBOn&+A2R9@?RcU}xY;b0AHV})EwlBK%YKq<EARZ3e<w#-K8joIquL!XGWXGkWEsI-AN>8-)t?chslDHv<hC7{N0D8p0Ufxe{o2AwPkKY~oJC87Iuc}|RVB1*-I8_S_v;8s<L1toGyr{c*x2203Ug)l6(cX_%T$3)e)Gf@8upo6H7HBv5_@n8%k+9()LlyiIJ+8)9GKq@Y7KK0|4VdlI1IVPj$iUD}jD)6r-)`IuQM;!8Xd-NvtH$d!MME*;ZUw&BHm_)acm6QjBL0w)1;6-vI7YpUmlSZa^$l_N^dGaM#i|v9XPq1UY89Rx&^_fv9n!I*N0*J$TBbgrco_cr&O_yBeWbm7e-?8L2>J;HUU`qmXsY$dFN46%57cajT8E|eJu7PFHXb|oKo)M%=m4C;Sfhr=c5>B7Z-_=yuKf4PD>C%v94+$FR&g~qAinb7LmP~uV-tL4^N@yH2e}ooA0hQ%!zI|EX;!|m;zK>|TjRiP-6EUdQZraOvL}Tp2BR7%t!l6ea|2j1FDvCoe+pU54Sbc!%T^8D*>i6tdn3HfULG2(T7z0=xv95e@$>bjLWnFg5ElIKtw<!+l+pY#NO4d)0E06(<6<|qV-{`@4MtuPBJM8e}a*C2-OOsw~y31V4!CyqvsYKSNO@b7S8}M0T^exf{OJEvHVxRYQnMP<a2`&|`1Y{MV@t+ElsVF5AswNPyZ{cTe$`Ypxk>TVY=o8KZyF0Pq4hM^rzV}jPE)}`cDIRHP?203Gmg@P|cHqjOI|$}Nh>1OM&5J`yK64pwg@5$w!o#?&JiUU<_S1Tm4!A_$^Mlis_qC)GVA|~#YX5lKFUx95l8EjZPj=g-<E#jUe^V+aDw_1Uv`*vsqjQl9HV^yq%4Ge1l{A^G=>Pdtfq8cZodIM>9uUN@(;o|Yf=yfrVTf6Z4tNG(=?(_ZW<Fj^o)(RltaFlhgyi5CIBJsk9Roz6Eu^KlV3<`lKB^?Wtis{}m+R50Ie+WmY-BoFEMLbwH9Lrfr1m!#6O?d8oD=coenAJnT+y*c1aGPt8Yo(DK-xWniMqqT-!iaLWbX<x$6$qtg^UaZCEk*7`Q&)w=i*xbFI>|pe#YTxJ-Uptv7x<nqT4n4wKFf9FVb23!vhATreEYn^1FqkMdJCnx;%ef=?>4yhD(&X%zyaxy=aNtl=jF6D??i5Pq&VK#sbX_FfpRn7E1>O{2KvGqz5b&buMwTiZdstU<s<1Chjscm<pMob8y9d0n6jaw7;t-q9$%fX%kACMZ}hMG#ba#<XgF)yd%L{WETJj?G8%vh4!mt5=^1q<4NNM^Q&~Ul(#n@F;Ot0i}rVMf+GDn**KfwVIkK!W9yqx*uju>50XZW_pI-<vjj(KJ$CW^<q4*Cgh06MRLM;wnL9_I@!VvF1)3(o;oQ@vV+_fSNWlwRvI@5C!a(VQ6yXI05NJ(*^m(FMy#*;Ub#aX`9dCCXo$X{7Uw4jW4-!hpoSDqnTFk!<=UOd74#)l5mKYUd4V=CJ{z-W-ly_vVUm2CX1=L4F&>_XHJzl>B(1#e4VZSutzokH^YrJBdU$C25y<Vv8Frhed$R#T(Pq{>BLeEK|)nK5`@Kf%@0M8x&m!ELE4a_mv6&82oL^7y3M)`l$NSRz=GY>~H2{+qiJvV7lU=)j4#stHf+Tx3^a52fXbbr#>o2pk65)k7~T%%{!P->TH#`bBHB}yWjzpkmhu)kF&Fb%hkH75>jgnAsUdQ>d38exer-XxK)uZ=It0AVZce#QmU+ds3Kws%3bOwY1x+c}Rvkly(Mz!B;z1@Qg#I$4}n(f;pailAwPvwCO&hkTNn#4eifaMyEl0A5irbk}Gkoz`d6P6dZHTyGz)cGy1$iMpen24NXj@JPVDuhA2o5G_EVNE){Debn9Frriw4?HNKUXUCHzsU#>h$@|n#t;W=N1qox0xLy?S*Ita7s#YqrG=KV6zx|R;n9|5h=QOKQ@2h%?Sk@UIpiP`r5SyNEz}jUh&L@^ECFP#wtE9&xkS`Uht;0rU>76%(aOeg?S>y%u)ETbNhD1W`{9lvYLNWu}b!gTdMT~==mhIwAGzOclZAP7}=q)%o3@JDiKPbiwtYc>~fa1DOw51)@&?DGryj}PUE+Ghkl#U-^Bn4gnOG9yi{6h6jhI4d>?QTOIp6>J%Ap2RXZB0rC%2D}R%E4on1;x-a%x_)bUrLa<l&i=~a$>_*wY+JGe_U)l1-S%sorqFH_~D>28b1)U<H}c4{coD!zGryXW2iiH$=h$}1jqA7CAw^VOEt{wD?Wno&=hyv%&6>>VK|Ks%K6Jo?Wy_TRF39_Ny{^MZp#Z<z%Gn_>6-MrBrR9pCpc)%q^k9PQ!y~_BFJTzpOP3TpU}bWwL_eTz$rlCGX^UOOiM)<HCPIG6>NRynn_<PkvZr?7V%K)MxthM!av`QI{IreKxgWoGOqE42Q{(+P`@ss0}FdYG$!l-fl6gIDYhNsZ^zceX4*_yk>pyMBvi*}pRm)86ze2Qd>oMF8nR~?emXMw_NjVPCUNreFF1OSy<X_*%_z(V`Ndb+8!Gb3f~jKO4nWTbQ2Xd+tDvNr1S)u17S<3Q=CYM77_DRYJ;s_rgpinTe6F->=297}P1H}r&>>>=@-LT?Va4?WF7}aZvHJe;F7*#_Hbp`f5ILC^NLtM9ZN4jlt&h_Mula7*Yb^Oif`E4MrmKm6N08Ub>|3~Q9CXC%p^}#nf^4CRiZ_Vc@+q$47!Da3Er@<>9>mP46?R;H#eqvZ!)8qE4VZ*o1)l~H30z-+Z8zfP_=cvWvuGkCQ^x?6MR908du}6=E&tLB7!Qv<$)(_0cp(^Ww^IT7Jz0CL3$UH1(^B~JXyDwhwa?t5sp8F5b}dU!4P)SqMdB$_Wvr})(W}Wp(zUDOc$1047%o|vn_Nt^ZDCl(_hK}nTfb>^w9m*WFr2M(>00ro8hy#zu@Y)A*%PEI6N#2%#+K{YDHe&<@(Z}{pU?z&E%keD%D}Nchs}13G#As9q`9O_Wbo>jj=-4@6{=-`;6=bu$lzlMSVep#{J*{D-y50Rl4I63@78QkOOYWpfD@^=tf`!t&^XM6Tzs{q7h4Pe$QE;;efL=+nK6!U(oa2iFO&k`QB<;XZ=1!uP1*mGZtzgO33;n0!gE5{c@)%qJ2ha-E2|>R`4go6>av*vPl^kMOVC_^vmZ^to0V>P3y8uLSh~&j?mp2x#yal)Byn-RorGh^ulJMv!}&XwqL=E=>&#P4=opK}S6@&>bd&%&No5PFYeZ><s{P@0=3Dg?%F<Z8FV#FWm8Tv@T{(f5h2&LL@5x3q5mjg8ATC*1S$2=X7vd1AMa8y(+c-|pcsK>xJGp8=&h#y<=Un(nNZoB;A{&mTO3~dl1KvV0Dk7|9+|RKr^va1FpVXtH_sN||y<f;S_xm{3-GPYNB>rHAV5-v5z>2Cxq9;PIY>sM!Y@L_9Sp7=8smV-JA!erOUQv&L2{~1NW2$POe^Ibk^|-=DE~>MV2HHmOhBnXNt`Q5=4*uZs6L`)nJz7m27sj78e!gENDwSOTT(jA@&kiSihV3{8xU=O=mGF@_I_Qw!{t@m1G|oR<vH@U#p2{K5B^@9?Qy_V~;*U7|<i0<JY)@oEw2vE}Mtv-HjpAI;^nTcx9$@Z2QJfOGC(qvOd=L7XkNz#Em-^qm$gf!vHSq!{kn-I181d384xpG-Gu5{@K&PfBn?ePPJ@medzgJ{&3sv_~^1EN3eSgxy0WV?hnS{2NaGwXmc9bn~i}zt3NN!}$nDhA2RW`g|NOX3>#SCH7ZriJnb3D`FmH8@TfBG-V|Gh$+APo|74s!Rkg?4$(7tI<#ADcNu7tVgU7u3Z7PGp1Mqn<J~-Frm7=3@?j){9?|?|hue2lH(LlHF^S&l$(P`lTP$>>jKQa2r?Tolh5XspJX?y}e|oOg`=k45F5>?`K|E!Ut3F8`7xmtFr-no2199Muh~s`jWqd{X_~f3SPd~@(O}%0(t023#eR)eAd!>d}*6&lFkLXx2omQgG1MIq|jCJPK=$r?!@n|7F%^Th+Ut!kMN-elX~oFSYMvH;P(M*IQ-(7jPMWz-V7WJ!y=2D1_7YN0$UQp*)+=nk@%fA|J5}7ANY0^!>7sVXMmBqwO&Et3`O7o)eIN<V8y+vA6xXk5v^8od=-NJ9=@25OdY4nl1MCgUmHyQ8lv(crsD_3av5T2)DHx=NXx5CwG=?Nb9*W7IdoQ7RBs$m!Px;^)nYJN;pe$P?U67waW8~zdHV-C$`z;fv=%O#voDy)P(|l+jIuASerx!Tf+{qc|03Nfik7c-v;W)Sd=+o_IhXAb18*NJWHOH7{}k3@^XU?HpAE7hJX46vpP%?jZDl^4Vp7j@4j%BU;fgjOVmast&z!tVf-^(#RLers@F}}0j1|&-Ui+vpEKnXPCIRwfp^)z1hL_})oOlzzMaUPZd3&Qy-#C2&=-@FcZePjXR!hz30GC8w0!lhQwV6mjkqiyv_m4$C!n|;3W0^I;L1KFEu5W&d`+KO|Z_+@t{M>U9A&Z7-^gt_68vd!`0Jk{1({=!=FS>$grMf#j{5?zS4pcx10pwMPQqbVqLIZ^O?Tz*MfpV9`R1%Bb54QlsIB}0*uXv5fH5=;t7~yq7;x2+%n=PW+<14dKdf96ei<~Gf#R6;Xz?O<Yi6x_W8vpCe80|$b93t&BZ2!^0Zp9f7Yj8OnZB|3ykrGdD#C9=;O~5H66mfiRvH^Eq?;wxqzzN2)ypzPAuX1BT`FUaq$MHt&(PkL1QV1*Wg8aoP>P}AK`d3N=AdwAn8H(NCtE)@c;7y(Ie+@W6&!L}DY$6<jnEIrgv<n?pM=D@2&F-MKM7Vl0O`^uIDpi_I<osba49IA+N3c5Vh1?{iZW}O++u4yU-jYuLKm>z0%Wi0nwR<0sGuvlrfP!$L57T)b8x+n(Xa3SGRBF|A4~qIZBe*jI4v=r~?h_k77b5_4Qp@egcd;18#mdm}3I<?C`W#2cTuoxlSMwY^6a}VsWiw2%v{q#Oce+Uqp$_vy8YxD(D80OovSft{BV%W&yLdQE%gX{q%M73gJI%qA0yizm3f7<XZ8&SB187PeFAEa^&mB5)kNgI_1_h`IpH-O(B#11#w-4*+6ky>TiH{tTVF=b9YnHs|6Q5ZVFt`ESne%9&{rqsJox#lFrfcrY`d#k6Zk<#NfG`P;9QVh1JYgTmX=Qeh|GR9@V9bJk{FK-x&e+^(nV-tsh<Y_x@(F1)%6N)UQIJUr;f(b~AO7RJg`{!$DIY$^q%X;C@WKA8vZsb6$v$L~%Jc1r5nT?P#EfP7^L&qWmM1mJ?89Rn@x$sXbZRM>;c^_9BT;%G3+15U9PoMh{XZ|nyu@AfHGX1G<SSHXay`U13E0V3Dn>Pbqqh60!ab4d8(D9KkC@aeNY%WKi80!p41`ykO$AywW^JlEk;p;n?OKJDRm1qBc4Y6~3?jhG5GZ9AR~SZO!5)sF81-hHwPkO(axOKXW39Nh1sz3wb(&%r;&kZktkiux?YL6DCjPE=KpWJx!y*8+q6Xlw@1cfb;^A%3^3s*foqs|ef*tP<&up|yY^QoTu$<6Za{N-FpUZKXkq~>!C7!o&5#q2WD%Xb30{|FANC7a79d&F8m#UOdpeTO|ba4lSo1O+0attYc;NanK>N<N0dUP^NNZP^9ecCu)AbN|-Ey$75F6IQeam}}<g#u7^)x|&MW|ltFx+zimRQIVv_UAONXB)Gr(}Q007zMb2ddIQFnunS%Wf$E}{-2z`dRf}v?184vb5r#SL#fCDu~aCz9C|Y*1bh}i6QJl4sWptV?<C9(D!+Wj5!2+bLh<N2i+2A{z4x8;Q2_%oK2&{l%g&`e3t;w-#^6MpVp<3#{^fFaulvI+6A^Od-H@m`VL57t?3bTR4|W~`2ES&Ly>z;R^fn06=Ig@Z0CaU-CkDv{1I=So=OV+3SWph7+X@kLqObWsf&(-CfO;#umq(TSIEtJ+tOI-)GeQPyUrw$Vduez*9un*Ozz+_z10UBrI9`PUDb!kzI1DN39M)#GBjT&5v{&=`=any6?<`n>xmzMmuz*$}?!QbK&mLUUj6>KX{}F^_@d%OkeKn;{U@R2p57g3%qO7XjeU$*Sc&z*Z4-!fbN5&V2<|ai&g$pHFBA>h#uAlor0&-QzhP)2MB-#5uZ6%6_gKpM4#W+2UEf2EwgaPJ<LbMYq<whyX+;V)HIelNLj7BlFcN+^n#n7;VNF8^F76eWn^jx;~%Km_WS^*NYO~4R>_#slQgUXZFVNPM+GnUZkdNA0^1_pDpAE1978fi)bEtG6h)1#whP|I9}Cb}FRR!)M^Ng@64_-swkz!*!}1{<_Ud0QuK8`BQY@VZlMXQW|9ZMMRKbyo@)qbeBdp3Rn?OPx?|z$BcVmu^2@GGXC&_wDJ^_?7<NDf@a|iHtGtphLCV+jzy}|4@ptu8D@HaUz7sh(2)&87>o{mI+AeG*Sdr%1k*M4~u&_V$_!=K?1=21cFcLC(h{XsieBQ5&%R%8bmkK8>an}?E?_m7IZX;Uwjr)=M>29s!j`{xB4<Vt>(nTg^4Q8!a7emJ0I*%hw~x?;O+>wwiLwIZ@N27DxUSz+A3io?wtuAQ}@h_Ec$bJ4Z(s^>;*QI!g%s@lLKtY@xXvnqYRT4W05v{F$$DN-b_5XrP)0Q@xgLQmWP?ghn5AXMeiLdC-0oG+a=I@1(lsgl>raQGb&yyjnxE?m=!W~>u~_o>_}gBXcay#5^>uy&}w8of+pEeWSJG=2kMKi`W~AO>7m_W0T(8->#@ZqRO8er%dA^yJ(psLje{m4Sj435dG1%#t2b!INLEw#@qA_^Aq%!^C)XA3y>mX%<QV!E1NYT8OKF32rGCLj71pa=VWoj<T7RT^q7*N<e;8fs;zm9<{8rVXtk=1e?8M&?V};?24|)Z@K@FO+g)L(T(Y!LTx8n?fDQ%#F=#@wD=hn^t9|)h=h>E=w{ArLMd_t5y8s~|;nT{M3cvL_%cV_hVe)mGf>B#*(6+tW={}H+Y?d78y!JJQ(`E;=OOUy7`kLaRY40Zok@FA;Hh#PsdQcXUWuM{&`KoG04ToK`Sm-0#!lF1s5l6r6sG``HC$@+YkMaZ~E=OoJ4JYRHvOE_ke8Yx1;a@90zgO~pXmVDGIXpeVZ(#T5gPTACV&>L59DsJ`Ym`HXyscZ%>N<QCKCW+c}Gz{86;vh5&kY(AnfF9vjEd#+I#tFB%?zO-0EQvM2$T+kJl%+74l*FnI$7xvqsaOUcoKimZkl;)T18#Gt3YXFmm*kXHp<N|#eD)w8ARb>H3Vafb`<coRNMlF|&vtu6sN<iIWbxYtqoH7#SRL@;>gnRN%m+m2)?Gn384^erRh_)}YkOro>1IAA^C;S!J75Nk;VJCiPgVu~5ct82z23&76$%B<&jA{hUjN!+M2B@sCb$Z)n5L>*%pifUbF;?yw>_k*b!IOZ!^_?u^uu}%{4bZ<2xqknB(-4W7GvO7ymk>Y?WB{^L^b31bk?aey=pwSy#t+&^HlLe-uqSV9MO+IwxW3ar_wug)qgZfYWUgH1Jg{1jL4VMRUua`r#N_D!<q5=C5E3EoJ#zor5i))4Su(y^MY)A1IGJ3MMMUlDs#?05#$x)xQ)KlZk9mY`)$Oc$JwB6<cx!bN(@Sc7KzLU-ij{Vxp-M~AZL78aN*?!Tpev<8Lcux47(3YO8RBy5R|WqrF@o#%s9zunV8zyX&g!!KRzG_%S<%FZE?39>LmAD^Q1@Db0ziOojj}<026&5``|im-@yNM+gMN2y8j(BbQ$VEWupaf0{X4MNZiK`X~p6&6@lNugOU$GFT`o)FYqf|CmBftewYl`?G(Yu+H@ZFLyKqKz&<Q<Ws{26dBOP$(4XFu=M=2JDin%qTaad~I*`rWbeLZM(rLp#CaPx(BRYP%;@{X4qNS|U#e@kk$u&6YF;ecfYaVBTzU3I$=JMKrLK&3fq(ew<t^x%E%mlY3^Wkq^!hx_GEz>&*5_E|+NB3-ZnN6ND9EYm4bZC^lY?zO_n6D^6G85v7m}L7PznOn;V|`0CoVOD@bMB=)WUN!WzjM$ibQUuFQA<uY6^5gG2u0YI@EGlq6$LhhppjSaK9s9MJ^=*%@7@$y0kzOZ5QS@9I`fdQ*^_4#zlLkvjBhT%gw}d7zuHne5~b*<G*m}d!R=oNg~PWm`YFQYR%!jeaVRcx_Hvc4+ClQ|aP%JkfT#{yHR>yNy-StSOsezoX9gacV7?UEV?keLxMSZg<j=ZZPZ$TWslf3Nf6VoaBFv*ZWFCq9Imc^UNhO}7K_9VX$5FU}N8xXi$!zUNjrk>%mW93;;czbL5R?#E%M%He_;U99Q{f7qLoYv%jFm@aoDb&@HBD3xd<JyRA&Xm9#pgG7L0IZt4xd|kc{FLME_7EB)$0c9SOkwMrt7O~n}TLW<OYt-xiR^Ug(OS$jmR;XwrMr4bOtzmiu)o+^KBTSYzqy_$lzkd{z)|;e2U^c{~R+esqq@st4)x{p7dG0#PkAAd!y~Mqdt|IDyV?a11p?~n1sOQE<(Ya-oBU9y1B>XA%D+~Y$NfspPlxl+We{VC{msdsIz6+kWGF@PL6uwP&jx-yPbchcu<S*cFsdM6>&+3kQ5^ZvHAVs#$|x?+DOw%Mwgrv@q(Cc$!Q{5FL>&8m)iSrR9?G};t3P5{_p#O0#Xc0<O^VHDjUR&l3JD}UK&aKYvl(xr{%<vF=`$Th{T`fEPb6+N;n#sPL|1d@-8HCxJCI(+UNs&n1NgGSe3AKNFmx8JDQI!aO?C~;W{D?{6h~m@6lVON*L=Yz7m(X-0Hi}-~mYwzJ|t&d%kkfJXZ#|7#BI0@kG1m`P6Fn@fxdAO*-k1>%_b3uK)k*QeSvyJ_+pPcER41+;bJ``i!hb>ni0NdeFPZEzOtQQ;j%Ig=a&)8PQlp7R5fx)q++l6`EWwL4-T8O)7^sFgiFdC!QCn$30!OX~KvAf*NbfP9@nq9cS%}n58iYi;RLRrPr3pSC5y2`Pogt$n+k`{Yq(B&_1Y@4{lE=;uo4&^1HFOBbvm|L`n9rmx(3WK>`5%ItZSTapS>ix*PaqTA!bT%9a5b^RL-Rkt==vF!(%FKu|RkR8(p5&fT;$3Pd_fYBOwbMIco?=m9k6VUhX%)wi~tXyhR)5FnrxMYIwY3YCPdL-X{07I7(Rn4EWblcwpX;_&`9is?XWpTX9RZ?bOk0nJc}&%y<=!oaMao6A)(^0BqUCYXp}^+@b5!T63oH$EyKPN`n_`QK0^FpV%OPAmkMTNZc59j+j{IG=^FKZ}uy2P>t#1h#XbCKYytXtAb*`d(oy!+YWhF{BQj=bp*79<9MbWu!3Yq9EGo!%UIHhJ&)2r?YwPv<1RursTbSY%W8=l=N9q0ghO+et;*byTtBK-=6;S793(hJSn`qfG_U@>y{>cOXmSltH!#@-OR@|R*u?B-js056sRu?f}MY0UKZJ}Ps6e0WZGgO+)Kl<|JQ2V=Ndb_6?2&PM9&xxMpayQ?{Ux+WO->k(@WEXFTIbzL^RfD4yy<YU1&TkTbcZwLagpPQ3|WUHsA+?{~_*!o9QD+M4T8Qjqzss#;+JmLo@x)i0zhXBHTASv>PLL%+|<9_=HMg?}oCGDewyTtz=S38L!8kHAeYBA-=a8%j)-wED4<?VY&>7)wil&;u({nhyqC)WZ5&-)0(n*y;t{LNZ<eclzg}T2cjPVnNR_z2qror(cHS)?P5cc#=L(Auj?G3`xVj>A95^S=r4a7q(+*4Q<|NjaZD;2PD{ILb1rK^+#wx3S!qCz6YDkuMR;~z0O<NcVwV%<JO0#Sw~RQhILtFjXFov^t__edJCHnECaiXXRyMk_BunZBZsWeZX&;kb#|mf!sK<lTneN=m(ftZq5Xc0T8JKwa_VPj9il@81yFciVm{rnb2;f-eWxN}_%IA7{1Pr%H0rd5X_Uz+wR6_#}?fZ|$HoC&@w$QaD&)2XK@_L&qEQry99d?mIUPw@ag2Q`V@Qgfp1u>0R+RIrIKF%J<kmg)OBWtNisbzvhZI*5P;r3n8m3*RVbX=h#AX4p(VAbg^aMA#8Za+W0wDQ$U6)U`_hLb@oq$#C>^EUP{Vf&V)JtM0@sp?E+haoh&w5Hr?jZ_y&7SY&u1aJv0pNG|jx-hzE+Jf(a8}++DgE!&E#*jzZ#Oc8<sXLYyxM0W!+b-vH{4nZSd)Tb*%N|dr#UD6Ay9c}P0Hr~M^aIJa*LDmDJDe0?D|H&kf6Fv>@KfL=Z`lt3du7K~BP43mB;b&KPURjCl@#7xo&2VgLicKbJN6vgnzJV5x7pep3Z~s}jRt@&?iM@)UukLeU%7%2psX%ht&H_eZUXkyC+h?pQ}q9E{>N}wUHu{WI-uG)sA$idBP)=F`1=^z6r}FX8NmDyKwU1-W{f@#AA-sPb-tCGNT!3ff(|LYtv3X~u+4gNNOULW&EHECWo-yCby;-q=C5d}OqA~g<eY=QATqIcnHL;B<wXI3aMv*^yQW^YsKC*jJ{0t;#nvhuuqe??{~1y|<>Y)2E3%~?UVa6RT8Lf01dTO4TIC~lt;TiOL3w9o{>MRwF?)uOX)muUV<9k01^Y0o$ONoIf-WqO+PlHK?b~*DMGcQJep`8cjyhpofO~iyWbmiGQ$5GYKrBa+)GvCkqe2*n@#(6^Kc#1HSft>o9j_N7@=i>z<}=fyI|!xuIuitcg81|faYxYfe^&ON=^$fktBsep1|5!p|4HZ-pZtt1tuVKlnl7TLfz<Ly4^H*s&nzjRkn?n(e^NM^Fzp@tWvnt9PsAtV!~c^1fj~(cBZqW>=9TdTPfh>!HnUSd@8b-WWl0PIs&Avs;Co`Js6ens&}&fcT{9-L>XX~Ax%@}(Lc4g{;xfJuuRM{;U0vb&5FXvm`3;nfqW?R3M5T+2W|{s=K|Fmhct4x823pDpgS*miGucY08PJ3v|K+25i`j;j@NC!~?$gl$<%%?f*ivqSJP9aHY-0kRlD|>38gF=#uysx1mqk^2tz}_a({Q__;;rH<jgZvy9gOLekBfw&^=|A8>yRd)lcHmBR7Yba$BYoPZg+SVo{L^hKhO2_3cW?Tz*Yb2h&zoCaT@cpd0ZVRSdjW!!NF2%1UI5YJWQ3{!D+R#AN+$ez!6Z@&CP(O;ho$o3)O42ONp8Bz3d%-7_{r!*&IWg)HUhHT_hYlV&2lZe|fJxao+?0IVZZUvR$N)Ie+FMeLzlcd_(4xGJMJniM*$NXC$m0zTqXf%Yw(LXX}A$tCA(t6PIvGc+;gySe~;XxH1&l>sysitb>8y8cnSPkj2VwC#*y<OzCrOOLwMz5kF&QdT(XwyKv0s-^c&B=B_al=aaX<-8b3#q>7b1MG=J6hj=3vJzCYKL<!#big;7`3cY|t-;<6o`Ln&S$P-)_Zr4B5hE%6Qj<rgmSrYu8-R!>0l~3*pQ$$o?5oW2rX%aYckdR(Jf=Iaiu7=32sf!mstFF|}A=Zk(0UJ2mikV2wtq6l2>(dYB1YQhVitC~}3K1WrP?-3?B#l>JefaNVALf&|KAKH&S%T~B>Jp5oRUUV(5f!nc03ezBfUap4012@vgkOeGGcR}W@^!+e{N=s8^>g57m3n4et34;D+<aQ^6sfC{#A)gyCS+1Z@avR22HWI@$Y<>tQ72n7y*QYSPJ&DEi|HiwXgK6*eR&oOl#eFRRY7;JX%IuFn|Ts)I^0ZQd&>&G!j0QNT)Bj31wAitYGv<EkrS&o=cWE|U&wsiK%1zY9Y%%~m&TrIp0X>$EStc~CuAK%tQ?9BXFk37tZVZEw>D+<F5S(e4K+T}-pjbG&fv(_CFH@l#NE&Bx!}1$F!(%k3jCO*WJ}{~9w|(+%_R{gERd21RA<`aLW4gTtGrEKu&-^~gaE?r<631tQt&r}Kn8vJ4JUtZ$ws;Hn%FO?*IgZdY)Vebag&<4GJS>vN<hps01|OE77Ju~T2ej&7fV<g<M$l1FH@zfp}L!I#zfwAUWXzZ7M4cYsM~jtL9MPn`gcZOV0a`)Mr1}wPPa}6tEoL?#Cv4TA=>R&)(1z_i|p#Wo)C))T@q<+9l&ze@@abHpfsxz(?${wXX_#5;Ham;;78@)bfd1`I<~<_7+c-%db>wBIIc$W?zIsH?A*Gxbv|m)e0Kn^C@3K?b3b3Gu~=RG3OpLIewwNAPa3v4(4eeFb&h0QTh_?^448#yP+0?nj{WF4I1zxrY=q9efujHi#B}8zNEzJ2jksPQsDEXmQ^11uKX`+xl8;GYEyo=#Tz9f5n&&s>j6S7dzjV)TciQkKW8^uWv7{oiWnb;fe!+$icyQO0FYo;wKA|qE@%h$j`T9*(nj44fb+uRsiTT7U43KtN(O|-H4sCi(eNHzrR@r;-9t^8H{kFs&9GvvahkQ<Vea~JtI$>g$C&CYUZP@cpY5h_PPgtX%orZs5c{I)=UwI*#h6P2P(Q>vNA!%B9qWG<b<~7cG+k7(HihOPV5tr1}&;zH(ff4u{uZo03r?;J<WROss-B8u>j<!V=>2~RRC`jDoq6s?+=1#aA*Xz~d9A=kU_M8E>do@9OWY)@#NJmq<N65=J)%?TRZDZ%xo2lskj+vcaCuEaI@pa2NmNAW3m}q0l!|Wg+ZBvmzIM9jgyrFXsdeb!}x;pWj3W+tvn_m6&Kfl=QP4sM#A_M+}c(cKkFNe+HW-1@O<EnkO-c`1i(O>U7Cht_hFWpn|74(Zr@Tj<hLW~~uZ3g<u8#L<44YTRyii*VNF_{xx1e#QKXy83~Z}<NMJ@g|t#%9G4_>a&@IbjM8cgIxFhrhSIQl)Bk2^Na3?WAeO4Ev$8cSIVOl?>8k9gdD)ZPkq+feR#IHzOQga)mHGIw<UMv-?`KRdH>?btUO|`aS0@LUC(BDi{^Q9}(Cw#TcCtamJfVLqD8}2`M{JiN`qR=9L`1`$$mV<oN@s;FTFXyjmowOKMVv3@;>f;;A80{L<ecZ@1E3aEtm9;t9^{?3Ph(ygSef6o|KYuBCU^*IWR9-0AvZcNh&vv0Sy6RSTY;ZT<vrCa9BmieU(x>6TT<qZ4?ZSX!XKi)j)^dhbtDTPc8y(ut+7tq!Hf%YafZSx?X8Rm$^w=NCxTUucqcT!hEi0OXaH1DzNL_W_t(3R9gBp)ypz&j9IVK9`*W87i0);DIQEBDmW0&+Jf2Nm`&sum28)hV}8$lrYrYORU#*AzM15QmD3dhqctF+g!9$;iSv~xD2Ei+bbVwB8bYf1-(d~R0<BNlU9{;2@StZHYtaIUvkHzA2lW);IKel{$YAKf=rx?bp~xjIIO0p#kK^y*C2Z&=s3$OYy8g<0xT|SD%wxb?k?8n%2zMuZ#v!Cr>LIuEbGjUAS=$w1SNxG@xGU?lyu{BaMGpl5?mfYOEC$vL~^fGdZ{e|foi%U=F*rAZ}1A;&r40(#@CfY`vv%G6vqzr!u)MgPO%xQXJ8_Zs@XSimDO~$M*J$aIiNi+s;`N%G7mWO&?c^yi#qgUvL=IH?0VcW>*Ym6$U@Xq8iEzT=P@*%s-@*FBlQ*x>^5A(*B(Eg71;#FUw`ZnGj4ak_IkiQ&61WD)0QoF$o<O>L^-@KMU~skOw3#yl4?REIiM{E1PBgCd2Zs3FNiQ{MUi)*ADE^FbY2}K=x}n)Q>!=H2NylPnkgsgIl8CJ7;?DVs?FSf=>e8XKwU`TE?pm9q-45+b5-7sI?wrNko#slF(PAK5eMpO6k>}-21}RA7gGFXU-FBC(()9Y=&=pg{>*#Y&R04;`OZe*q2mZ-A8Q==66Z{(0g1S08>3sA{9QyGl<qMI0J`P?BlpQF>_e5~xMm&9GGeB2%sxk;%{vMNBE({$LC&#p*0&QsmlzxU3QstxOuprrBzl-^8WR~&usFL>_-Remo4hj~u7JF$D4XP<c!xT6`o-2ury6*BOu};CMEm<|wkqIxSfnHbQdG!h<^zO!ol4kH=N3G%wmG+W+UEq}vw4A}E!QX1amMY1WUEmS10?1TDEsYY7OA<V6suAHbK)21a*zo$B61}5gexw!>*`I!su*MaEg2H>;ESy8U;*M)j{|Xs@DjDTLiMbY2$tP=EVeL`*VOFP@E5MY<u~#N)d#EsP&m5U=1ylvXkXb*iWDKZUCZkMn+1Ql>ucLC^CG&Q7V`RW1y`ZSmDAkzMNrzuGmAv?3H}yp)VZ-(qh&Rb2?C0#t~)YVW2;_Fu(KH-1eHxHx;AwrCSnl%uD`{sHSZQNP|)o>wrv=m_=R$wVDSXV+a0%XT}adOi;0TisN>?vetli5{4vPfPH!X}&mmyN_^b{-O<X9MT>wJCpY_pd!0s)h`5pX|c>WnITn6H;F$rCfIkf@{^hvJYW#%SHrP?MmK~%A`>Zl&&5y8?8|6GdzyqGfGJ&Hu=FEn8H00PY{uWC?bDfOaxYRy!StS4Pxrx1H|63lyQcZ?7Pg#M!U*ZpKwt;IvbHqk)j3e<`?mmO;YZs-7SQB7@8>kX;3t!S<mcP=Ed(S=$LcweAYs6U2F9<4EpidndEQ*UA2xPjsoi<>7W0(JbdZ3`?#44RZkhnxpQjX4?WsAB8Y!|p7bCvCJeZEU*U7C0jMKUZbGkbJl<i*XOiiT-!#n%VmF>r12nr(k=uf#aSbkgTq)$F=*w)!z;khRIt{=a=mCaRJL|*~3FQu^@kxpKxHlI3V&wl-nrqEk!G+KJ(^G<xoC@h3*x7KyBIIuoZA*IZ;mqio#^*H8M$qQ6s0n&Px7U_X*|gR7{>Qf6tj4s;q8CcvmH@gUj-bqLNA)&B8x?U9!bwK*|z#6GQ}k=bsNSWkmSm>|7F`Ct*o}wRES)B6%hUYLgNrG0bGrgs-Bq>JJsm3}@KFBx-|W%2Kue-{Fle?*Pd^O1y9X_z9M-d+ncXX)X+ggqWLHE5G#6e!Fw<AV<}0T+p<#oaz=EaWN8vWaN*Pg7athFwm4MLmo8<gS+tssjQ#tRm%kacnrLO<h`T)yBlWCZWIzuWsLxbb7jy?DwWU0=u#+Wwv_v;@S2{WNdpG@lHw`xeLqv2L@+7NAnQ}kty;lg?~~uhA*GSs7{@jF*^iq34T2tP7C(YN>n>70g>iOtBaS+&40C>{H#qC5n(zoC;G^9)V>mY}&hT$+7(TKjqGe}HqiJ>6evT=SYX@IpSa@=HOdD(tNT_hcs`6L#)dAgB;{QYqVv$r|yu`ejY3WBt0>iR)+Yq^xrmoY5&d+h#$H8|N;3XC5tV%RM=plFzBz1}-e-MOd0H-UmAO1X5xV8mjkfG=c_Hnmymh&Or#5rz_Cubo@E_L=oN1qx|c;vsGy&3i!kjwIWO_x57t-|}_)?ROV!1B==_`x~#yi(u=gP*0(x)c+mqc>FBBNm3|=(gFWvwF$R_c&(|GWi`Z7JP>fen2NBPZw;0J6JkdwPD%KK4UvsE}m>KEk?I_Sn|SBk~{&@-J1M`r1-%urY;I{h=u4j=dq#@^;7-r5GugS#g`K3DcGIW{WMqa+wgmF@DcMeY9owQhxY8+X1-xKIp04o)DEJNytZrNwoA<(7bTW(bRd)(JG@EzWi38S8JPhd=sbg<HCdkvqjR8YJtuS}o1@<1REgdudJ<P7&U_#bs7OeH%$m0+J+9fLYV^BkO_7n~kksAQuMdcTbbaI5g9mW2A0F;XM8sJCv~HHUg=^H|d{2UnEWJb>Q%c@BN`!LT;U0@eTEn@vRw#>pW5k7{x;A-o%mD!^H{5U68x0f(omIf7M<oKKQQIR2+F=Qsn1X$u6ObC5SM6SAWSrrHsZ=zF9m{zDbif5{iQ2eSA03A+vtndCi9-iL`1O+E3^;c+wBDc?5PqOHRwa4ic0{A=fNJdrSn-vlY5&~e=Hr2t9O?{$*{dGRs1emIk6En^9nEvTDiK*guLvE<V*f*&O_UjVWjO^mZauhORGlI{@^BB+<x0dh!aF@uC@0AIk~`oITq|B?Iz+@9k+dOhNFUo$o(sOaBCLvnio^OsE!lL<{U+N(gJ}+>i-Dk%9ag{6_O@kb^7~Ib=pUirsN&EhHmp+}TC)iHY%GL)SKX`aO49<#;oMq*LiR&N(PkZwV_`{B9FC+}un$8{;4SQilw@|7WH8v7T(oVt#^qB7mSDc3ZMiT)PPc|_J0>!&1kO|an>u3>HlGPD<(1sw`+&GE$VfR`>M~gQ!yU`P>q^Mj9i{%rB-j{$yD6f^qm=>Ltv<Bm-cH+vLQZZuK5tLD+7b{VBc?2ouL#v7QgD}vDe0wcz<^wdQ{KvI)=CLSt_l2}i|Em}U!#O{Io|L0Ttn*lYhml9EWc|*MKn^Wr!b8#61$0#>4&{kLrw+?u4R;BHR6|E>^Nb4YUKLQ*$MT;c1c|C?<`iM7oF-_6<F3oqRWUPGqTB~!aYB|()<*-F=E=Gq$CY^P%XD>THKp?_s5LF1yKij8A67US)~+?m-)0#VSZtiC0OWmRe;UFZYl#2Vk-gXE%2I%n`_3)U1E8HnLG2YZ$L^dn`f|hho*0XTBzr&`v2do3&PRibr^jOxb)cJN;Ck~_(DPlALl#XRWBGQ)uOmK`O-Ks_r1qB=$I`p`YZ9^WsiIIH%^&i{eqrp(fUqh<_^<f7i+*Vq*RriKT_<mj+wzXX1c>2EbBVocV$qz2M+{c8s@B)EA|g*48hVu=yGdo)pSf8IAV-xdQW@3kx;-_wTN%1*8lBBAm_ceLP)5sn_>G9>s58%?@O~de?wD8gK)r7!5<HXo&bn>>NA&L0sBdNKaeL4jU%6#agLHnaiqX|3^Eo%Yz04t3PY9k_oM`x`J(PHIN?C{gm_ziw1Ru9@s~QuxtmiB^VgXK6v&yuZyM&@x+qxsma58+wdgaKCg9_7YH0!DUu|usfDoQ;$ygO#Wm36I-U&{=GmwL70XxN@v3@}mPteuH(uxV?75hw$`q=f===Ix2ww1Z0(XrouiIEjken@N-%y}LlPbPPCZYW<)jgvu2@n-eB?u$Hz;sB~jQPlo2EVB#7_+(iaTJMG`%{nM_I+UJ(T#z_P4d18GuG~|R;mBJE_poPt-3)$>>N*4iuUIe)%RO|+y=pErY*7<|7^k`Fdqdhc%b<Yn?eL;DUdWWytEBNzhI=r?UXK6&%4&6mzh|xTz@o7D1(q}tzqri&6xYIwKI2#Zh*(%z<hkMCeTXieSFZ!Y8S=ruw3E^-P1ir6twTqx+JPuI1Oe23kMq6^g-PMKjmb5#t|T^pZ_iF_B}DKy6=PjRJ*MdKv^jvUR_#7+c4RC6WZyW-{&~Xiul8&%RXU3~-VdeP<k#<w42V)G4>3&y5PLlofFJBiZ20Ywi`U?AY{3BJA7$9Ho%UljmW8l~T{NVa%BU9-v$=e?&zLBT!^OtQS5pQ8*BCZ?HoBJKlVm6~A_W0mZB4f@>6#vnD)YY$W8wD4((GV;S(7;1;iE0s!tJLfqyTMG@QDg5+|zcPMH<aa#}%L3dH!8F^Xt9GzGM!_4*&SApw}0R*>Mlj*RBp<W~6=;)1}sX8AAVj9Kz0&17@-~iq&m%m=Er`PCKq$Ek>ejzN&ap6l1tP{*aoLj~-$C7;IN?99t5@zP(zN-O9`tDYEq*b94IER>}=nbQ`;$_?tQ?l}r4;ox(ArZv+d@lDTLXb6Q6JV*l;~>}rt+^0d~z#8Bk5vTM=^!l9>@1!IgT5cMuoDtM%zw2V69SCbiW$8hv?X2s<x9%EBRT#1ZgW25e33Do@0Bb?$8K{RQ~NZqLlg4>t;_996<RpZqhy?+fLCdikfY-7}@eRPO9He-bNtUMR_wj!2?AEZpU64i1=g6U!v3z9Ng(Pr0bB$7?=4MHL5g_p`;Lg;nJOn0WN6vNgl3{VI$J$Hz`H}+9;K)-M>NPuHxhB+`iq)ELgx3WPf2aCcPn2aNC_xn1d0Ruy*4pFsfFnK@kOBd-n)uVPhrOF%-@Y*XxD1~Q*>i_YlG*1oTPC_<8WaVk$4sXi(Dw26ASsrV!F73zEy!52A%d2e8YqIWAZIiKwB`yl31#)Th7`joC4fW`K$aKQ@4$3RhiWB=Zc{G5U@Bx1{6MkfFi9>`D&&#zTwl1ON_sDH749~<y!(ar2v+l>cCH`Ycj+B5?R~ii>fe9lLx@?rCrR|=4-E5m0#U)}y){C74m*k#xTwU&L<%>$L!1$97SA5CT+%p7>3Dj=FKfnS{)-RCfS#HvEY;yr#y7rV0ddfDh!_Y27gSHr<J|tsX31zbqojba1UmwYf@{Sx^E@&6Rcj|p-=pWB}(tfr?Vf_tPvMl{I;`U$?3<wLmB=b>&jG3@A#Nscc#KkOwoP=`ndQ|GNc%C8!Mt(0;b+V7!1<|?<V(ToquL->-1WVqzjFgEZ$b#3wtHWZMsp)VIR8L^vPR}fcX7oQN!UL}z%%j0;72G5<RR{w9yXpK$6li9}86bdcCvBRof6~>}p2^m|FLoAXuPG*6xZ?p|E=Ldx-qqbkXIzSsQ;CDZ#=+ZilEvQyGCrZZflkUhno~lT->p>Lg7RWr_2+n6rPG1`ZU&Nkd@`8DvQ)X%9&6R2ec^iro;;G!d=>+fIh#($#b^7MPZ6>9|LizkrlG=VC(r?mP3!1LH^qqzpXWM?`2yM0ox7Imtv$Hm#cW^#?NbK*142P*9ptiLI;?7FOnzRO(*6VL1)9s2adYCCLR%t5pIjD~@bDv4rc{`_7$bK%Hh(y^ohb7K$?~5ds!h02p>mSL+4BUTp3(74_lnQL&E_@Hz0=^S1z0+{;(cRl!??F^4fYD<9f65G*nqs}lr1Q(ea}*dri<CTBLfO5O-K|e*g_igmXbCnCgoEN;{p2qV)0QW3CM*0n{zF6C)Z{~IaOuD^2^kv&<RR>nZqkaP-?K-rZIo2@TYRcyx^64L8d&gjnS`k@;ueWFpP6%oG}`0nDL8dx}yEOp37u+TP3=)fJlq4;$%%N3lnC%gozO-P`!<WFaJ2r-Lzr`G2HGw#87!b>wJjK-u}^52kud;1lz`dT&dAjvkH#WTR>^dc;6`COke5?Fimc3{w#Ecl{OwMa}XbLN}^w*1ip}`e~ta8@*|xxl!^t;b)q;8N=5kS>E#XQ(%23!cGczRoA$XN7G#^|_1*vv?uzcI!m2)hJtcKUodQa_Wa&GSy1`%KdZ1uhz`v1?p`6b-B|(`HL2v!k!%YnMqR-q8y#NK4|E?_4^m06srL5`7rSb0Lo7NO|3TZtX<eyjAI}<oOyKV_y`No`ROhl^p{fL;r#VN?ajG!PL7YTI|usr`%jujztaYpG7M*%@1&n>;Y4i;^Z&@~1$dz?&H%s#;KGmUmiJ4Z9Wejmky+C6(DU>Bg-V|UQ{fzWb#sz8D6M9{bE#KNwqq6*TZjzz0%KXyQbXT}pVpPqJ+sRI1)jmUl+75Up)WVm0(#V>i5mFBM=0`9DP&B0p7R~dQ~l9IU<Gxp&L{}Qq*YL2|aGBtwied6p$4UPvGAHh`&P#N?j%jzA-`HREi*e^CP|1Y8*i4LwP8E)^Zg(k&c1`@JdEOw|5<}S+k)r3VNN?ogEb2wlF(n1RYADs_`zPice#~#$XWn+l+$dK^8=<s{f%5MRR)MAZmIh{`#7LM))knAo2*S$-}a$lyBM_`R+giLXub4^GDd0Iq#YCf7D^5O;ZE}SP~IzZCd`(r+6wHf`X8!Y*Pa3j@0YHj}b+oVaFOmlF7v$<NY5gnz3<jpCb!OdE1W|t}6n;_{`78~0=PyN`KWklNp<;jFmVw<LS?D?Np`~&EUMklU;DsOSZIkwDExJp&R-T^gsIN&`N9b<MBy~Pw~5lQ~uyie5oU_*!6ISu3)G!*R+{8ZW~j_HbwLbjtLLMooH9~LrA8~%k0Vf2;!xT-@R;H6pZv)z|0JqRsOwNq2&ZW*g^>k7NqCuozk0i*PMN8S>%zSe<S-!^ll(fJw~wXC0v>TL1WY7g5GY^YKE_H_%#-E`!9#uvXx34w8AJ`az5$#!yOihy9O#c_50zdY(XY!dE^V?=RTvzrt9v8!u4!?NW=6y)Z+HIh(h`YmnXAT&E<hj_*NZT7BDm=Nrcz;!Pv2RVY=9=|~V!|Ig@tuYUir5gtD@h+weH{FM2je*+Bxw2=NF%EcFjya!aej*4H0-svLWDEVCFFkL3nQh^{`@93dX|1|P6nH!Ixj1f^kfV~7hCnK~1$j_-eCivV@On`_mxKG70^L{u1)%cmpA$24@@mY+AKYtQb3+N{7H@uJU(b}?y5Z!|8%_5tm5~CAJ*VCZw&AerJQ~1oX@P*Dv)(g1)@_TE-5Pudfnb{*Js4ihDUj`7(A#S&DzwRs*w&?SL!?KxVs-!_jZ1||@qjL_e=~<-7vkpx6MXswc4j1a4zH86^3f1#adC*9=6LOWe>RoB?Q(B0^7KvTc@y{Du0*<rL<3-pb(jmiKb{kn(p;`Wf~7=y6f>$dA49x-B5(JEy-65y#C!<LK?u?;7XsT14{9w6!d~~@#jt|p=s{}`r*>?gIRbVm(92O|<M0GG%dooLF_+{u41aam^09{kA?$^u&`dMwc|P<YAkjb+<cWLwkr-yUI(E=T(%Azs*FHRHJCoJWxP(raTa<N=)61es*i0*#9%UYE#NWNzvh#CZ0-5TZo|1cqv7+QM=Sl5BeaupQ!J4}sd}5YUYIcF*k&VhW<s}>q1AF4{Np;w;oRIN-iHg0NOIO8c@wYaQI=+Cjwq$J?NO7m<v$ziO*>)k$a+rL0(r+fu?QGjmPO#bZ-4PGZ-HWsgK$LeIzg&WnFsfNb1dU2a^%Z0`r^uqWrZU>@96I?*RA(UF$X@MEz^KkeljvRNwwR$;{QYB7d_vNA1D9o!Gj!zbX8x5zVQi9BXJK+AK23rpZSDdO%s$=@f6VtI54M+W$FaSA=2%pC(%_`)J?CYRx}3dZQk03nlL+;u5-y!9+$f82k|TCgVcodX#)BEMn3dK1o@UY8(VOf)YIGt_2wToY4<ZDPTjzj^TJ?07wKBMqzEL!7JH0<+jzyhtgQZ5f8L?=&%KjvBobpgf|F>pC9mQ4`9<3u$GP3*R37_wqX)7Q6EMVDTRo`2Tm|dgT+w<^fV{qvL<wSlhJ>{EEWgZ4~q<3QMa(Ci6r$$%XI>}`SrBT2$SPMG&rzg(g;FRt5uS3|bJ$2n+K5JmpqrdhD5YV2*+;;RSKsI!UH@@)L;c2H-I-iX{QI}ez#xlVOD<p2PbrkDDIGcy<pbNFG1P6C!+JQP-EZ{4x#>e!sLxull>7{7}<h3Pe!&xW&n4pTj>q=-EZ_a}ikaZMtjvKQ{JaU(5F=gFuutjPPH0DnJt+q5&?rZm%b0Ht6J%4rkgc^mqyLNl9N^S87#@SVnB0IPBhL4@5AwG!{%}+jD$00ZG=)tqKhWMb6yD)-=yQ3C&=W|1|7A-ZS$%GRd;uzPZ!=?u(qVt{1t$8Z%eVu~o_*|VGaPQPFzhrV|ub}eq!79rz@4~HhH^+GvZAco`c#3}ljB@q|c(#&^?;n)<83Y<`7Z_OI?8CzjZhi-HthSdoef={O@*w|<@IN7t^2+R=#;|pLUN=9S(@>Az8&W9wqfAPBlbvPW%oiR0;>p9T`r0HwoaG|qW$l`9lESIu(*_yRZqh4p%DH1dJ+O+rAY3C*KABH>+I%lFMUKJkddRDqVob5TcxnM1vApIRds4?y%*1Bqvar`R=V}pQ*d*n(F&Fns1ubxKqQj1Uy!{8um|RGsn|;0g$UHZo{zL_=v(YSKWnu`ay`XTg_*-sKaEXS?HJi?O^HpcfEXW;^iCy@D(yFwHc3a<cAQ2{o8ph=L^0OyP>5Z8z19fWhvLgXVNyd|ndjXWnDOd15!=*btMUs8^)@4oaoq|_a3Fw7yoEtS4$SGi3`=IZPRsdQC+E32c$aPjoH}wwC`)Vn4kR~y1P~l>4qqPap*PX1L&?2#mjM|X8ZyIV#_jor5Dhad&^#tkH%?rCc8FMRVZDyLf=N=xOSX_T1dO{l6DU#C{DaO=!9a1nrw8Z}t+4}qVe2nu`bP9V~1RwfVdjI_~=k;8xb}~FtiO*IkV#Krvlk2K}xF3h6MQkqC*pO={P~Rr}B#!F?m7<{}=WvT-*yN1?Lc?}jP6QWEa^bjNzpQd=Fz(33Ibr(q`D=8k>l0+m=a1pY`!vh^S%vMn&CB)SY#s=QKFDt*i5-oYvek}sD%2gStf!d>QyNxjq^}v_3_cs*j-lSQPm4x~x&v}QQ}!I-iZveIiYP=qpQZ0#`;V3m8q`!FO+AVO()i*&Fz>dCI?QmL@wd9%u!su5Y)JK4j~yGYViBME8Nr4`1VD|8qiE-_gij2gewrf>80IQGNAZG9kh4NWeZ^HTG}?;ZdOwFUv*^Q#xonKsQjc!}=am*t&}|Paq#(g7I*`;RvfN<^D_lH2V7P*^A@Us7S~bsV(YTr<>k0y92sd+jJ$LX&F2MMZhin1PZ}J{;Bjrb}P69S>>=2M^hL^NXno8f3e-M=6R|i{9Q<^YF?T9qXz^v{b>hL#{4{MS76!2{+dq(pjl$J0(91<ZcC051xvx8WJlsU?&Lo0+=N{8QP{rOE|`Gde;V)K*osJ|EUuA>n>2t4q3Nw?f2Lw;Z$IY&huc$@Jn&O079$4Pnwzt1Qjee~%M=1=nRWtEdUgSjWt4C!7s()480D^oC|cABg8`Y}%NJG*~Ji}+A#O!<P3zK2&A<E465q#GRUi;@;f)8(1}2^l~LN-7TI@(h;BQewW##IYjXv0`SdEGIl2ap)8zWwblAK3dp%8;BoL)0vZW8A!6(huJ0uc@T6Ye^lmc--8&sedv7K@^KvqcMUcfHnESMp!BWmP%GcrDL)dj(&H*-9`HL;hK%1WMOW0eXnvntpy}~z@`&crE{lw?;&*Gch8UUme!d@+FEWyyT!w?sC@f3s%HuX00G6>+)mR|=yz@|Zrt-qu`X)3QCPv%vum^3umtIghHa{EHxqt2BlGv(XN3airt=gO*A+t>+l|#EG;|`($QNKD^*5sq#nSf^0+`?6$sOEs<T^Lju#l}O(EPD@9MBz5E^8)n-9{B?CzDU_ju+nFqU;hVGQxsTA0&Ry`lv}gi#lzQ}j@|=&POcCJaXH(Tz@qMHU7{v}qgDC%YV8!L0n7?gaets1Z^b>vq%=a5QYm6km}LfHe)QBiS8j_G)!9!Glq~PfAq79<s<Mcq-oI)<AAi<*{9Zx+<Czu|QBG7EXH?IEv4FbV{wT>lQxR8@DtMFTOT%i!$3+Ch_!3cN^+$;f6DdgKV54oax$gXb>k{F=yUI3MOvLFnni_7NU#-0?Nk5bWdiVfAto;JMHnK`JZiue4@>+r1YnvnmiszxrQ$IcQp?dQ1)nYXy5XMP@bH5e^loaM3E0~Pu4owPqhu|{i5sZq3%MP#RG7x?oyzDWM2}W?DZET*-r9q<##H3TeEIiBO%mNu=m+MzAU<5sds71*r=F$V@$6V2fX}@h5%$ZI!ylSyF)b6{{;?D2WY`B?5Q{;#CYpdGZkHom^dLtK%HtMWs>>txs%w?x7cavJ>D78J>HY>lgy{4sN+#ttDblw?I1b5d=o#)-z76QJT{{wbv4leqCG_<J*Kq`Phvk)ZOT$WioK3Ym6w2+%2zDkqpI*P<`{IdN6LL~zV%~qqn{u8R<BA8+7cmUC@$dw*%)_T;YMdWG7E1@tR$HTAKpV|2ua=-8v&KT<2f9S2*SXaFnv$NiL>YkP;yk6Q+e!EHu6VD85=bU`0S&Mngd@EEMjmGGC%|1ZK^7*Vi#?AGavhK`Ppl*V$GQR~M1n6yp9T04u0)tHswCjkJtL4kQM`xwr`gSy#GLEp?UqrgGR>ZLnOPPVY8Y3;0IX@yV?fr-;_AlOQ#d%|fK35t#=RM*FHU5v^@doe@w!R_4g!2n%t3<%{UCW@smNJi>MupcLi7eGf1>J+n0ko6iK6c+o=F!Q+`RT(<aes+6>1bVK<g>f$1O69T=D~3?&7>L>P7ccud()6>CH7(y5WOo=q8d~LF|ze3@aBkf7&vc=bW{+i1V8Yal`&}Mr|hDB&b<#nR^c)Qt`EYytt@Z@54)5>r#SxJLh;xVldINYZ4<AeHgS|UTzQbUiYKc`8_XO(=~p$j!N|h=hqju272&p7EF2-rI=gltC3<x#mVkYugNDG3vZ6ptWk^1;Er+EGRmc(PT1l?~!e_V;B6$lOz)o^x`0fI{1tUT#cX?Wq9>K5#j<d12Lu!~Ck9k_FmVsY*Be}q-6Ir_V^Zcq8#!VHS0JaFy>}PskeK%#7W>wVKumFC7L3PVd!D)1xuvTQ0{9}L9u?0${`o1%ZWvlW%Wv-UR%(f%$gQzT`?lG<=;GA=p%;K@#eTcLPc7J!n(Bp39A0%4d6)hTG?QYH$=}<rjBXM#${*nR^{Sx^h0NKzQi~O);l!dfS5B)dpYpAg2rk`M*EsMW*+x}XUnl`EDVHoQ9QJY;bbxKV8X=jm1u^p1#<4z|<V5XBu`4ZKG#PvHR5}IPj4uXQ9QeLaoYq+B3H`~KoIM{F6cliIkT7ySK2=o|ZfM1;PbNWEm`gj*lBYad+1=|?9GfZUE*tmo}HzYAxy&vfafii<w*3a_}QUvFnuG<P~0U2MqsnpHnlc;pFGWnC1G#*=z;mvd;UL>OIFa>pBS9M5%B1Cn*(e$Fh=M#`SE_3f1Rua*>aYilAQ-TO(1tMMm1cx{)XNi3YSa}1pwHwvsoycCO^w0ljTl!9ya(s?YFZ2%0U$p73T?dx5IYzlUTn8JThW$`F$UHBd^5>jVO_P1Psq7mzD|(>5Y{%#+cr~W2i%4Io8~{fci!<L^FYsZA$`#C|JK}a``Cc=Asi(p9<SQ$~f>&*6nNedp91#YPE-{5IGRL8YW;b&417@kIapedEH@D!ArB=s6mG%4;TU_0gLxz%DF#uSM=g^sFW50UaU!0*WJ8CWbZiiytLDl1j>$6X?O?~}Z`8;jl-6l3^8Ar+0lUv_Z40iNQ0>Mw5FmwcvzptMp4rCr$nP4Znt*w%tLy}JP0V0&0SvvPXk;JkzsEN$GguA5N8@P@fqyoteQ@D>gAo(*Ix)JLL?3_2Qq=p@eh|U0t${q5}T&)9cO!8Jc0>tEwNxjR*mcnn-bQ=fmuE1v8yD|=%{tUZ;j6AeSiR<rE;gkx$^@=8zwj>b5(W63eaZk3C{UunqVS05BS05zFDKq387hiZBT4S{a)@O>je+Xk}HmMH<cC;k0FV)e!TV17ycj(}6?`|s=^sfMSj+#Hd1CES$Gs8_q6Z3{0$F(f@E-O`{YqU}a$N$jaQ~cU;2zT(H6ND_3&Fd1|7N&rP*&FEd6(~!y?+3%g!Ap8@mbNuj*C25nH^GTayN~lp{E@u%_h_o1pjp-x$m5Ps3@<l}`Ye8kk?q6csBPam*uwn$PV$}#7mK6L#JG`4kS?F+NS+FPC%s1k*|rn&EqrlHVO9$tMV}_;XgYX!VDfk4tNm(`-w~6I%7%pVy}=Wo+$T?%1!(k9%fmIo3LS&Yy1~52>-h74mvcnV>a@Q@zfh)%dIcV?&~!)*taN-&#DGi3pPf5am;%cgM@3U&J4*gn(h0e-6)IHc)9SuAL5P`5H@@IKJi|Q)DwN*q!XpO2bsv`-V;YoAFY_>24xC8$6-brXMF^Uh(bwnfHSh={(_B)kC);{~)1=m<&7ynRnAaBJBaMWnz`XfhcB391N^qn#^S%f8#Z3)Q!6|I7olW?)!nyX|R(#>lj$<<&E)aS+Jm3aNIY$0X%X60SCe2<}^M(}2z#!8Hsz=<@3)=OOQ66$ik>JaGV%9c?Q(-kN>hY>U^c2hyZF3UnhkX;D9`vu@GKQLUSRQkypFQ$I1?+-KX559UW{giVw39qs$3X4wzvr($`{oWFHm5B|M=q{jroni1jFlcO@k1r0BCAdq$3q2VS}W4$pkI^BlNEZ;w@uyT=-thwRIQSzG<%Eo3n!JV22>0Y&c%W7TQlc%VD<Lk)hx~Ly=BUbplj24Nxh^ynteXqg#=5LCkPqF+FO8)<s(xL%R#?jzjqg4TadiE($7QUdvj{5>Ni<Y?o7*fHO$JdTnKmh0~;Ya4#GS;W7mM%e_&zleaB5{6vu=sp>Z&NlgT|gt@~M~Uv}KOlKB=5$xDqTzcgv>xPN2?fl~uA78A%TNfA&7l;QnGOunp*wpQI`E{w)BiE*w-8y*xlG6Iad)~!x?C%W3AcMZ_)$nSZ0Rk@*gqyN%C_B;VcKfe$_w6(F(Zp6MvtNSo3M6xyjE^HG)YVY8LF15p!6)Vn3`<yF_MouFSi{Je(`q6|>wxF>#I>2>)_<vTn4oB<vg1>eDATmEtK_$37XBW%&L`6uK?I9~~8~7#=8gm5V0?m4|zH?kgae0=G!iG%Jk#)dp@BSb+$vj93VzN?C-=*BztCmWK#CEsEU2hb7Xh*UH6A>>djQ?5yflYD!ZZy}moi23MnWh7G{ddF--0>%afoUsW9O@A;uAcb;hC`>RAoJbnt2t9rE|ngFqRyOsg}O2^LZANGwyUp;c`Q(OqXuI>sO&yrE`|`;gpBXH4!w@Ox;Q}<uR%qQ1uteFf@<f7=uxS2Nn}>i+N;>I^-)(mRYlG1d{V0-50p*UKD$%;GLIO=$3U$-VFADsX;<L$uzM#ft~+`rfATDpdofiNxtC!+*?{)jIPG{DM*wNY{<?{cr+G~)Og#{qV-t+l)sJ?|Ngq3Hr$}d~1OX>jXl~d0gdxHu(Z>3SRG}z<tz7B_qaT+c2N^Uh1<1P25WvF=#h_6H`N>|*8Vy;x&R(8C&FwR06%5<beqh;|L(pnL%U3unF$QA<N672O*HIj-GHfgIM|CE3hqaj4$c>=sfSMMYlgL2i-&!B*E)>#;A+oO<9V?VoT38@b5aSl^t>&faD`<9sM^-r6n0Z)ul2t4%5&1T5hYPaS8WWARZyvhAj(kxJqbHF)B+iHY+%)SHCw$1sXnqCX(d<vtr-`A-R)W^~<6pG3&QE@72X5W&^5pjc72FPxIYzk>yvNzPe;&<IBJRf*(C;rYp>2p9Op&NGR5%`<>VG0)AF0DpWzx3;%?@g=Q+l!p2Kk~a>myr~0xEXlvq3<9PV0PsqNX<b&zxB*Fs6I8eUqGFXIN;bVGRkxe*B%F87B~D<jO&Wg~H~szzj*~d~AjUlOfe(8ulRL)T3o!vas#514Da-K~3lL&yx)K<(yH978Z3+`BR8&y+~IZGq^cLAGP5=EbSaliP(s!T);l>n)AW)#m3(}OMO$^#)nWE>0s1o7JHO3?iss&wY6H8ZKOI@YTGHOc^MyJY~+lI?7nIyXgd?b@t;kJ4>Xod$??k8;HjGp3_+Go8iPc^=6(IhHxA!Vg|q(ABEcR*!sq7t`^wrXh1i^u$26l;fB^vgC+LP57;*|51#*I)7R^pD&PhldXyu0n1fsGt=lT6N-&%m*vi&v=IQMX>x+ZW%w(s4WtX^;qC$k^1uy#L`y#$l%^QecE-MT}#|2tU|59!#X*J~Xjs8PBbvW`9wfOi9pgCW$IQ|%x|?X6(Oj^Y}J_lL6VMC0hF`(4L#XZWcE`|OEWV@M95r!lW?28;0#^z1;PW#Do(1NZH1#QP$k$(>o0TuK#N=S_YIgNuvJ^uuHQlq$FY<RF0QlQl|hTd=%k36=zldC_^AJq<LypbPEOy6+*^`y=X1Icl@mIH8Ndj<?mzn0Ae91NeF>OmC-oI4;VreUOwdT(!Q3f}9*guFfnI?+tZi2PR8gdVfN&5x9U4wF8fPn*OF&2HUAU#@Lqn7pF0_`na~n>_?3&v^Ks0X$_}YL$O6DdV2UB+4~bDa~;{~P8dMIjGmq68Kj{jwUXr_6>%y2AdHiucd^$0D8^Q4dcy!+wOm{n6-9TZ@--hO^D&#tW`^GDy{e!airrnJzanMD@zUvZd690A*bod^rfh=I77mDhXthP0<FIxwXY4jlmA{(yjWePDN|h=&rFzP2c~Nm3+}ix`#G?cH6$f+-UPWLNG>%f_w`<ZW7S*Oly?c&*k=4q^Uw4qfaKp-=bty$QE0NHb(I}<a{!vaUg0w|9Z`TQ+B-Jk|Rw=@%v_Tv04*oiIIZ`LIZPEmO1_ISQHNwL-a`;S1b|)WX_EVB`vd~n-CW={xr7f5?f0;qzcA8_S$4wdENh)J|ry1Of;{vq!bhb^x-KoHJLqgf?hX;CDnmNrc*PL%|y%JXRGq7v@7bj$1^d;Ykr+vH^Wh8hi_|Zg?28@GV?=!jRd#%Idrvzro>h$@H_gF|X@j57hdP*DE{9B%$ivdo7YA-n>)}VRq#5`X-sMzgiI>jr6RJm)Ni7F7^hxG0P`DSQ4yt9@hB8p-E?lV*Y^yA&CvdfZYU+=hR8p5A}#8~b)Q}*Jp1mVPZ9XYklA=<c0UJliH0EN}4AV=G-#Wt3)#^Dei*yAX&4gj7<!V;cOrG*qe!_E#*1BODgk-Fh=md$MwfrDudVnV|VDkV}#r6luvY|4@fW0^wnrcuk!1>=-Vx-*t%L|TwA_=xBQUUuMJ7*M#_-64a3?gPXE--iWH^3pW#NXnD?$|nn;C<OooxKG?Q04PcZ*SSF_%cK`(P>30$AiF9nyPgQLsHAh#RDW3(DApzI2bM^nMcRhY5I)C|dcfZ2XN)@oAjU`+*JvUpek|a!igmqL60tkfs!7t}e^6LSwc((ip@je`@?#JH$c@O*vJr)m3H7A>lf!P(4^mGy)l_h5fq1c0)16=L4NQ>TEivs_vfzBBhd2o;>F=AdCBD{Z;x*G@N3b)ffMUez|A@^5D~&Uw=4l;9>Dz@Ky2HyWL@F<B(&M{G9pW~2FliiyN@Gv{;O(j{$%v8LbrM#xpkh<Sk=J9DYrb=E=XnWwjmQvLdnbM820STEo_%z7B_*-ExAuS{-W@c19(A#6lkGr@4QJjhnm3ZrO52`7cL?_9Hx<J`!FFRpH42OsGVxn4Uq;Wt)hw{?^yhKlH6LEq%RUs(27Do%2i-rTlaPv$o~%a#hL|OT3^>6TMT-M$#{Gwcg|BSYg#VpvRduPwRTcns=h`Fd*$9KmX73Mn({+Ko*s6-y=Q_vHBgPHydcF6-@gX?3AIQL)_D0a|>72BMf?q$uv~XH+>xXv#5EoM;hFusIYLs_Uvr!h}gAq+P6RF?2V%90Cba<xnVGrf3n&-j9m`4vym|#q{2-c=;_TZ#?-QhD?LlbK3SiVCN9l;B_3_q1@!kdPfFjo@E&@9A|g=?}E_+6=3^6web&+l{!C)_fPEC-^=9hPMvHvBh*%#o_JXx108+UmS-bPi*B9mrhq0-nq(shJ+48yrnGFg}n)E2-bYLmVUk>h{+`H_Bj%NXkuW^o%s0SN~9OD|KlnY|ULbeh4g}OT{OV&H-yWTnnjO3sp^d7~{+Ag}K`~^D$Ay<^e%Q95)FGD+_R@9C<W2Kb{o(n3+>=vViqf0{lYV073Inw6OQA(tHL)BqfP>e_q(?r_G-<OY!Hr!-O{SYKD8G6MDpZFGKRj=Hez~t*jsc4KK{?h;XzDTsuO%nJ6e6gQWv;V1Ra++_Hf3BZLdRl{V4+7@GK4Bp)U~1aX{ChAQnvAT^-~e6V^=Gz&gtQR0q__K<A=urXG=JA5k-79A})l;yzqkVm={9IM|y+8Spp1Y>J$(Dwf>H`MT(g2Ra={59AmrAsJxLyW%jh&AljkO;SE5y1YlPKL4*9%VMoXjt4k=~=8lhQ;8|bYgLw9+O}TPa#}zu;lI_GT#3)+ZzIkU0*5MZB+k8PWe<Mk&?i4E*L#Im)ZHY*`}?$R;MXR&wN$8FCKEmqxA|Cs*#Y>)<DK4s747tYE%IY)yj3CAg1stC%9cCUVIMNM1lJ{82nRLCTJNSW58ENhuz?#R~8$b^is@kc4MUO=IHOasJ;a32RXamnV#g)%uP<BKS|6SdBX$;uQ+j3QkB4E`%L8+?IR)O(+$oN+34$=@P~7~-a9TC*Py}XpQ!vLp8c%tt3DN8?G4T?#P~29NKj+<@LWBg=(_#T_My+0En%LvWmfEXVNu^%)MLkLgDMp=kLACg(4ppf0;@_;|G#kTIqUj<tFqxjz6$|&I5%fmD370)v;F<o=^qZx1K&U>!B`4^w1R}o8*LAjIhYmX2ePnMwevP^aZwMg?;&5<l{t+>YA$A%;0C*|;Onu7!~4asqEU1qQGwAjNKr1Ubr6B`0wb3<lExiK6v0d~0Rexi@H~?M0(j~h<-K{gYx9|NgnSKT$%x}imp4x6bXYJzb&fi$+mV$cz@4*gfKJHY2h+bB9aisq&4Euz3<v^4klGVN9D}t>U$pKw!(5&N6XP+^P}wBC@cUT5uQHnr^8$~NMF1x@i<uLxfN^|;vUY2=rt!(2cDdxCPT_vRqPtS2o%dJPZ6swj?3SK@ACm);jr388i~%?!BfYo0vf=w_fj7I8fbI!AQb|eLCzu#G&k1r*@X024*CH`N$U2n%&~jVyR@yUNLeh-{RHC&xdo0srBZyT&*AE5(HRM~VrAe0GM_cn~N<MMT+XQdLEN>xryos9Z^@ULT>m`CZsck8m9K6(i6`ZW({2rq7!0TZq=ieqGcCyz~F8Icdq-=fl&Aj#jVMcg^aDJ-G*ZPa+PTDVyWV-ftMz{QmI<@a}N)tZ9vZ`knaarwV8PJ`ZJvWr4@8H9MR8&8yweh#Ts;ng1&rhfE98xi@q(y$zA6lr%_`uKu2${zLykD_FQ$9pN_d^K58HE)LTg&M@YY48UTTF42w@#EqPdJ7!ZwbS9B~<Q8oF`e>F!efP2;h90A#k@Z$x#aY8v^!1^)peF^<zKaxQ5)6#T~WX++dova>BQ0FV}~B$hkSHOR`TX%@!K8NXEarTpkgF?C?biBMrFK8X~<4`HY_;z*Zdz?@G{VmJq7VV_7C`sv+<fhFi=m*6x?AFmTBKlsAyu^$baZwYVG0Nyel#q+i<b2u>NhdSk^7?!!_9XBtBfom-{4NJoj%fVyS6*2(rr@%b_`B*o2N%HNgGz@3|+pd`_1!wj2BuAb$V%dt76%nAfVw>FF`v7w>0BlKvp*1g1bS3eN!xEfi4A_wSiAZp3#v;@*yu!x1b)@C`0zTM3uqT~=CQ8mMZ%}hvbr2rnzgsKikGJlPI=2K@dMPD{ojEW#SWN<nJJURSm)b0{GO*Hrfud;g3Xf1f^tk^*C(djy)W~3{Xo(q7oPjQ40bfYO<Zsy~;We!b$BOvUo0qnc#7=ZJs@KTcMzj5ts|1yoqa8)HsO5tstV3f%Kf<uK!B1@&|aCQk)+Dp=n1P}idUzUhq83K9#Z?etRJX0QZ+wM?Eacj$po{9SlVJwO^sXP2@7q#b=4hi^bn-nQfg)|WJVj}D2)!hkTqmXK6>P`VX_c8A$@B{SjTi^&lIg{L(xx#N=%5ul!5-&3v4>eeUUOlQG$v-R{mJm$zSHs&Bek2P0C%iow`!OjOC|yumK(M<MciM0@0m2FpEt$U+2P8r1qR)N+E(fi84;oa?PzLXptfNE92cLuy*65$YifU>Not`dlR`rl$RJ>--;l+ThXG%es3$CPNz<M=n?GRbAgi#`$$``ZjMuc)j0RbhoWD{KV<l__E3Q5GuSAQ)sg}!h9op<=V*0DYR?h^F2bzH-CTyo+r?7d3uo6-?O%bu{OV;<u#%3q-j8|p577FW9U6#<{vjfS+LoYohEDS<Z)wjL=hzohrU)g~;kjDD9QZP`4WvT>Nb+fi=OFp_6{+O~!r%NCFHqn=Fyy#6d!5CR?(4<KrBQn}kMD)+8r-e83cSs6-bl^oeX$q&gZC(Pn-$7ias@-`H0P>OjwMNK@5bM?I!^1X#!2WCWE4geYC&M`6z)_2~6lrJ|-MOYjoLUL{3Om%iklszT0y*~Rd4!EBBsNG~tBI89k9-yDHcqh#`5sS|l<%>B5?=1AxiN<itPIgSmM^Z%=+>M|2H55}iCob>uSIr(Msz=p$m@=yd1=>dVW)&cW6OSAS4F+~4^?i{Q-mxk&1St=zBNrzIKKxp?mV-);VQ;BS=54}E^xJh*`1QCqH(T-31v_iO@{nkwB;SDrh5ST`SifV3`YvU_;AYI?QE5ItrgRp^x#CK-tc=UP4osFQcHKAyAP?$tcCe!mus?ag(N*>jyZe2E(-F+0ahFcmGu}ou*&(phHcG+q!u}L<r!NeV=3LnnhKJul&!v+}9qK3!dH;+45q6LSn~VP({MKI(=<A|Lt&dHeBAc1UzJS8jE`zPd9FMBkr|%7S$6f|kUo8iN&ql26(4eLR$M|U79<3iU%qn`nM`eWe|3Q?)u`y^_eq&cC9D*iX@WM`h>C;?|;72K3b-~l-+2wLm{Y`K-{i54qK8}zj0*^XfK8*d7nLu=yi7lkt;uB_QtW;^>URCIN<4S@-ueRcN-qXwpR#MW%0)w7(R)Ffrq_?yCg-9J>Sg)X_B}_HM(;q)M2_Hsd4#ORG_HDHV8n!y|4CU0-_StNkx4r;!sV3dBbSGJVd8O}&XbwLZD^v>~^Yj9+4QBD<q2s`850I->01i*a*N8GL%-eC+7z|bwCBxm4t*S9;N#sZ!y7tKmOPI~Ep@5(S&7J~cuvk}clXt5yv^yle@zAd>wFCD{k3kqXSU6{+I~mXGW#odl{6yxd0}vN~{RM2}GXvOjm?2F!RqjFCDo|?C*fgzRb1Ny5Vqwj82lmDP^qw4;qX)H}{>DA`gzt$U%<>(~WF4y?IMx@V*I;T0v&%jNt>Z^E`FmYnaHMRO%=Z&`AoWU|4@+Q+3m|F3jaWI2U8rAB&%n$e0FXWHE_2Xp5_XWf@W4I?6STN(^bEkae}c?_ycR84l}$jXO|Hl&ork?<64P~MXu0ZU%OLTkZ4$15SySjvrhp*)$D|3F$Gix8G8^Xd2BHROJB;U(F-n)nXIg&9N@UvQ81UIKa@w0zc2=-ptx3b8{^2|S^{nI!(~e`uCBgYsjWP3DaZl1%b3_>EujSIcRB<GfynU;MRR4tYzOLPgaQgIp3Jj6wR9niT$Zuek3W=*gjuRI{x5I6gxJQLN_10;78`%BOHFTF(b$hS%#&c@HG0zyCquv~yKB~ACLi_y>g|++RaG5Jx2!lXzi-!km8N=W6PMFaQ>qDW-yP=k++lw4k6}}bSqy#Ve>`&l)^dQiuD`tSkaJD@k{ADq08%*Bm9-J|<%p)XyRhS{iHtPJt;gg~q#vOjm+V^$Nk$e(jOdPf13<$N6%{vuoRvx)3WzukrMY%;9to!=PeL`s3%-5iUnh0yAmd5eU{e=gcJ~S9NF;ZWqPUP+{<(Ss15g>%xY%PKP(*#!nTN0ltdEuYihwd2AjW}{@146p=d`t0J)ieUPxKNd{JWyg3%=UznEYRswoFt26@Q8>k$8E+p81_r%Hnkd+avVpa!Yg=}`C`UcR6+3}{K=d(P;~=4XgFh5$8&P8-LcA>hmhNNuv*Ap1RSu;tFrTjU2GlCo-!YYI`4wVk%y=u(rR-267Tci){$K?7^}rgk9DXZ9_4@%&Y<H5&;Wr##IBA(#`++YK4Xdp>efxP3TLup!KC^CNvb(YgMjT}9`WzP3n?A>_SJe&Jy=n-vKnaEv#4G8-|*yXk*-58z@gOlz_F$mLICPRDY!rl^^oeZ(qC~K<I}G8(zES~Js-zsk|M)nGi))f?zyH`O|5j^w>@6e3sPOd9UVtr3~<Ra@WQGLMDOc5828^_ETy=GeO%I{B$kgY?#ZPdVAVG~s>+j4n|(xR?Q25+{LS$o_j$l6MZ#H-BM@e>i9;ZzFfjM!*JK$E+6t<ly0DfuCF)+!3gZ9D(AMGT#7es2Us;tfvEl_r&J}eqATHlMt^eIH6x%LGTvAe932lxFUWNcE2+rKJ4StS7ERR8Bpe6CLZq6r|wA>6FT-N-F9%Le2dnubBKZ8zI%1wnNj@ojIV>q!hP|O<ntN@%{RM>bxd*vEz)<p}g9}Z9|8U5JOGd<~Ud&vR1$A0bK4E4<mno{%ba&xv>ah8%5;BBgXk4mXh!^9M$>P<Q#oy%&`3+^=Sl%~9)wMnj{PJ*TboUp$<!yeXiJ*&)H2z>Q6ntzu}=EY-r$d~QHm@JA!1Q0PC>Wa^5*FIgxZpMDU0jz}KvyG@CO3E2nP}Lf-v3uFdaBYI}sV5-$CIfi@Xh>G+(nlvzyI7&c$?LP^K~hTGd|kXp-Z{c>9zbr$E6M4=*a!M~mj+;f;tqe~61U`5A#B=of_|$>ttS`aSlsy?@nD>5)-apMP(Et(t7)a|0f%fSGpI^LsH=_$P5@VwcdEn_9y&M45J4_PB2XPeC{E@nJ~1eo(Oe=xz^{b`X;Des+uuFk2EjFS(vLW&>ZHuzoV(~-q^8R|>=Zn_{W`a-p_^tVa)8=)h=#~BTXB#ho4B4HGq-&Tx#v^Z>VY(EvRVh_I_5X>{$kEbI7bh<69Y#ya#o43{Lz`A;Z0sOS2D-)5Tl%_bi<X6I(<K1pxcKqMKS$I;2-Rx3F$OvgJ+df`y^Kq1Tg2auX^W~|7RU!tj{YzXk`gw7Zg?q;KM>xwMPl6AuCtbGor8+>y^Qimx`RcGcQ6?ZY2RI5e3Asp5+i8L$xOob%Qzs{70nkX;<P;o7ldFDYt%3(McvyG8YX4`q%%hAG&hE12+x9zj5pjO@}G3a;qo`QJ!0DptHL&S=o!-$&jFBQU0)h{AMD|qlyUdy_|@=Q!E;C015p0=#70Pn`d^M;I5o5#NyM<9ky=udm|1^#kCRg2dG98lC=sP5Dl`PQ?&5@$L<Nh%D^cZ-xM*6Cr7K&!+iQ|lyn0%u`@n5{ybfDm6|e31^=)@Qt;9o%OjcwGM)1Kv4Q$?CX@B1Ee(%Ou(mYWM|)%TZR-!23oOvBMM0EJSSy<`gLeh3QinWhRvLLHDZ1!}XlAYJya}^BEDBoV@dM^F+{K-129LPoI+drJEVS;ahbWsd^Jsw>=lGE8Xr`XQ*VV2rV$_~@?DCv|pT?VMY-%e<{PP)&=i$4h=57<eP4QfWkB-04Wi`gAy=zy(Mgk1+p8H2`Z$Xj}B!lM~fQGG$atR1wmPbs=(Hmw>#Mn_*XLsb{?m!zoA#UV-zQiq`4fN3o<@xToV94pd9$6i$({JaRql{*{3>*+qy#ItpqdELBw)!f(df!F9ax6x2)TNd1br_D$d7zGy?W%iW<j$%O6~GNDL{#5p*tUJK)SW#TT7n#&I=i{I^C+Sgtxhn<RoLD#7Kc>-6c{SVa;29W)1d&cTH6zfS2@_lbuo7XkQ;ggOdJli*V?VEhD3=L!)C=72Bj!VxKT>(Uo+Tqs%hEu6Tms1=Za+)OH`lV<$BFtRa;f(p|&^V^~s5Wt$KN<dXsH8TFVe^TUQoleOXOM2kKuGB3@&{IvBR@fNjuV#wCP^RGkgEX`V)-oW)u7m43o3mjLb7Nc4;s+8#GPyp|;bL)Mfu-%c|WXibsN#Z>TcA^)z9_=B}gF3?yh%46`7!}kH0s5Sz#s-d8SzM7j%gYoe`)hWW4@i@kI<N8C0`abaZ>cnGCm*jI#WK;Vxw__N}ZnrlmLt!a5*vzqd1Jg6R@f-!-M4d%KhjZ<${ZF%UBV0^}HrU;H&~pl)i2!nJi<uGX^WM)g#EiKo1k1^3dCX&<73BMZ(2cm02f+=GY+Q<`g5THPx=xxDYiD=#SS``bBij60G-Y>qQ;-(#_V)8VSyb2JT4cIRrMeXM!x1BXQ+MOU;vCzXmQQo@&?Xa7XXT4xec#vAPBpA?rSL;XcyGE<S3po}O#Qy$k$!<z7ZzHt)NzE1ro~FCH93XxXX0doNmpMdc_|Q%pGCq7TcblMB*Q%f)$*8F6-IgI(C1DCVidK3he-eKJml<7<AMMjk-2)H^xrL$DJp6qop3zklOI63WhL$!qsyVEUW2_biSWSzk&_8v@_Etk(+OubFY?<NSXGT^Y2J_aH)>r+Y2+$;cO9-)DgRjV=G3JGrK`qIOMT&*B4U=o<F6bozKg)TJ0ZqJLe7CnILZeCTL9EMp4>^+d<Hs)g_SgEBu=g;wxK}yY&ji$AeWFJUV?5M!=`;{L{iPuh%1K+XH%?fG0iTDy$)Y)BX_**T2=yVn_X^vO}vw+x%rC`yPh%w<2Nz=u6;B1yraNJppvhRYdv=Z%vMHhI!`R%frjUq*rY~!tAr{bfUYKyE*mf0OsMB<r=H$U-vnt%BOp1*aN{2nt8@SW39|Tea`n)w00EDr0<W$FF3iy8vBYQl0ssI200dcD'
def certificates():
 decoder=lzma.LZMADecompressor(format=lzma.FORMAT_XZ,memlimit=64*1024*1024)
 raw=decoder.decompress(base64.b85decode(PACKED),max_length=200000)
 assert decoder.eof and not decoder.unused_data and len(raw)>=4
 n=int.from_bytes(raw[:4],'little');assert 0<n<len(raw)-4
 def wrap(data,version,encoding):return dict(version=version,encoding=encoding,data=base64.b64encode(zlib.compress(data,9)).decode())
 return wrap(raw[4:4+n],2,'zlib-base64-rup-ramsey84-nohints'),wrap(raw[4+n:],1,'zlib-base64-readonce-seven-two-pairs')
def midpoint_certificate():
 from math import isqrt
 a=F(6249,10000);D=F(111863,40000);wm=D/2+F(1,10**6);w0=F(25419,10000);X=D-w0;h=F(119340301699,10**12)
 pc=a/2-F(102,10**8)*w0-F(1,10**13);assert f.analytic_prefilter((h,a,wm,w0,X,pc))and D==w0+X<F(139829,50000)
 M=262144;scale=10**15;yscale=18*10**12;step=(w0-wm)/M
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
 assert total==242318959769412276707 and sum(right)==242319468752867876256
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
 assert suffix==0 and Jsum==549999513468646731662
 ratio=F(864,109)*h*step*F(total,scale);confidence=2*h*(F(17,18)/X-step*F(Jsum,scale))
 assert ratio==F(1785857496442746346634861447341648444053,1785856000000000000000000000000000000000)>1
 assert confidence==F(3754643731142310263976092824071804383363479,12017074176000000000000000000000000000000000)<pc
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
def decode_terminal_blob(envelope,version,encoding,maximum):
 assert set(envelope)=={'version','encoding','data'}
 assert envelope['version']==version and envelope['encoding']==encoding
 compressed=base64.b64decode(envelope['data'],validate=True)
 decoder=zlib.decompressobj();raw=decoder.decompress(compressed,maximum)
 assert decoder.eof and not decoder.unused_data and not decoder.unconsumed_tail
 return raw
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
class ReverseUnitPropagation:
 def __init__(self,roots,nv=95):
  self.nv=nv;self.clauses=[]
  self.watches={v:[] for v in range(-nv,nv+1) if v}
  self.units=set();self.visits=0
  for clause in roots:self.add(clause)
 def add(self,clause):
  assert clause and all(1<=abs(v)<=self.nv for v in clause)
  assert len({abs(v) for v in clause})==len(clause)
  clause=list(clause);index=len(self.clauses);self.clauses.append(clause)
  if len(clause)==1:self.units.add(clause[0])
  else:
   self.watches[clause[0]].append(index)
   self.watches[clause[1]].append(index)
 def proves(self,claim):
  assert all(1<=abs(v)<=self.nv for v in claim)
  assert len({abs(v) for v in claim})==len(claim)
  assignment=[0]*(self.nv+1);queue=[]
  def assign(literal):
   variable=abs(literal);sign=1 if literal>0 else -1
   if assignment[variable]:return assignment[variable]==sign
   assignment[variable]=sign;queue.append(literal);return True
  for literal in self.units:
   if not assign(literal):return True
  for literal in claim:
   if not assign(-literal):return True
  q=0
  while q<len(queue):
   false=-queue[q];q+=1;watch=self.watches[false];i=0
   while i<len(watch):
    self.visits+=1;index=watch[i];clause=self.clauses[index]
    if clause[0]==false:clause[0],clause[1]=clause[1],clause[0]
    assert clause[1]==false
    other=clause[0]
    value=assignment[abs(other)]*(1 if other>0 else -1)
    if value==1:i+=1;continue
    for k in range(2,len(clause)):
     literal=clause[k]
     if assignment[abs(literal)]*(1 if literal>0 else -1)!=-1:break
    else:
     if value==-1:return True
     assert assign(other);i+=1;continue
    clause[1],clause[k]=clause[k],clause[1]
    watch[i]=watch[-1];watch.pop()
    self.watches[clause[1]].append(index)
  return False
def ramsey84_rup_audit(envelope):
 started=time.monotonic()
 raw=decode_terminal_blob(envelope,2,'zlib-base64-rup-ramsey84-nohints',200000)
 assert raw[:5]==b'RU84\x04';position=5
 def integer():
  nonlocal position
  value=0
  for shift in range(0,35,7):
   assert position<len(raw)
   b=raw[position];position+=1;value |= (b&127)<<shift
   if b<128:return value
  raise AssertionError('oversized integer')
 cases=[]
 for code,roots,canonical_A,unit_count in ramsey84_normalized_roots():
  assert position+2<=len(raw) and struct.unpack_from('<H',raw,position)[0]==code
  position+=2;count=integer();assert 1<=count<=len(roots)
  selected=[];index=0
  for _ in range(count):
   delta=integer();assert delta>0;index+=delta
   assert index<=len(roots);selected.append(roots[index-1])
  checker=ReverseUnitPropagation(selected)
  lemmas=integer();assert 1<=lemmas<=50000
  for j in range(lemmas):
   size=integer();assert 0<=size<=95 and position+size<=len(raw)
   encoded=raw[position:position+size];position+=size
   assert all(v<190 for v in encoded)
   claim=[((v>>1)+1)*(-1 if v&1 else 1) for v in encoded]
   assert checker.proves(claim)
   if claim:
    assert j<lemmas-1;checker.add(claim)
   else:assert j==lemmas-1
  assert not claim
  cases.append({'code':code,'root_clauses':len(roots),'selected_roots':count,
     'derived_clauses':lemmas,'watch_visits':checker.visits,
     'canonical_A_values':canonical_A,'forced_units':unit_count})
 assert position==len(raw)
 return {'passed':True,'cases':cases,'raw_bytes':len(raw),
   'derived_clauses':sum(x['derived_clauses'] for x in cases),
   'elapsed_seconds':time.monotonic()-started}
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
 raw=decode_terminal_blob(envelope,1,'zlib-base64-readonce-seven-two-pairs',64000)
 assert raw[:4]==b'RO73' and struct.unpack_from('<H',raw,4)[0]==350
 position=6;covers=[]
 for _ in representatives:
  assert position<len(raw);count=raw[position];position+=1
  assert 1<=count<=255 and position+2*count<=len(raw)
  ids=list(struct.unpack_from('<'+'H'*count,raw,position));position+=2*count
  assert all(0<=i<len(planes) for i in ids);covers.append(ids)
 assert position==len(raw)
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
  for j in {j for i in ids for j in planes[i]}:
   a=legal[j]
   v=[signs[i] if not(a>>i&1) else full^signs[i] for i in range(4)]
   v.extend(full if a>>i&1 else 0 for i in (4,5,6))
   for (left,right),(aa,oo,xx) in zip(operations,op_masks):
    x,y=v[left],v[right]
    v.append(((x&y)&aa)|((x|y)&oo)|((x^y)&xx))
   values[j]=v[root]
  covered=0
  for i in ids:
   zero=one=full
   for j in planes[i]:one &=values[j];zero &=full^values[j]
   covered |=zero|one
  assert covered==full
 return {'passed':True,'shape_orbits':350,'assignments_each':assignments,
   'assignments_verified':350*assignments,
   'plane_witnesses':sum(map(len,covers)),'raw_bytes':len(raw),
   'elapsed_seconds':time.monotonic()-started}

def run():
 m=midpoint_certificate();p=parameter_certificate();r,q=certificates()
 return dict(all_passed=True,midpoint=m,parameters=p,ramsey84=ramsey84_rup_audit(r),readonce7=readonce_seven_audit(q),scope='Exact energy/confidence mesh, actual-module parameter probes, COMPLETE four-case Ramsey RUP refutations and complete read-once7 cover. No field tables, raw extractor outputs or accepted tree executed. Universal reductions supplied.')
if __name__=='__main__':
 request=json.loads(input())
 if request!={}:raise ValueError('empty object required')
 print(json.dumps(run(),sort_keys=True))
