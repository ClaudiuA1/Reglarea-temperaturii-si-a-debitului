close all
clc
clear
n=5.01;
% A7
% amplificator de putere
Kap=20+0.25*n;
Tap=1e-2;
% motorul de antrenare
K1=(0.25+Tap*n);
K2=5+0.1*n;
Tm1=0.05+2*1e-3*n;
Tm2=0.5+Tap*n;
% tahogeneratorul de masurat
Ktomega=0.1;
Ttomega=1e-2;
% transportoirul melcat
Tq1=5;
Ktm=0.12+Tap*n;
Tb=60;
% transportorul cu cupe
K=0.9;
T=10;
Omegatc=1;
tm=60/Omegatc;
% doza gravimetrica
Kg= 0.16;
Tg=2;
% ventilul pneumatic
Kpg= Tap;
KceKu=0.025;
Tv=4;
% cuptorul
Kc=200+2*n;
Kphit=0.8;
Kphiz=0.3;
Tphiz=120;
Tc=600+5*n;
Tauc=0.15*Tc;
Tphit=(100+2*n);
TauT=0.05*T;
% traductorul de temperatura
Kphim=0.16;
Tphim=4;
Kphic=0.1;
Tphic=16;
Tm=Tm2+Ttomega;
Km=K2*Ktomega;
%% 1. Calculul regulatorului cu metoda repartitiei pol-zero( Guilleman Truxel) necorectat
close all;
estp=0;
sigma=0.15;
tr=1.2;
wb=12;
estv=0.15;
cv=1/estv;
tita=abs(log(sigma))/sqrt(log(sigma)^2+pi^2);
wn=2*tita*cv;
wb=wn*sqrt(1-2*tita^2+sqrt(2-4*tita^2+4*tita^4));
tr=4/(tita*wn);

Hf=tf(Kap*K1*Km,conv(conv([Tap 1],[Tm1 1]),[Tm 1]));
Ho2=tf(wn*wn,[1 2*tita*wn wn*wn]);
Hd=tf(wn*wn,[1 2*tita*wn 0]);
HR1=tf(wn/(2*tita),[1/(2*tita*wn) 1 0]);
HR1=minreal(HR1/Hf); %primul HR1
HOinchis1=minreal(feedback(Hd,1));
figure
step(HOinchis1)
figure
lsim(HOinchis1,0:0.1:20,0:0.1:20)

%tita si wn nou
titap=0.7051;
wnp=11.811;

HR11=tf(wnp/(2*titap),[1/(2*titap*wnp) 1 0]);
HR11= HR11/Hf; %HR1'
Ho21=tf(wnp*wnp,[1 2*titap*wnp wnp*wnp]);
Hd11=tf(wnp*wnp,[1 2*titap*wnp 0]);
HOinchis2=minreal(feedback(Hd11,1));

Hf11=tf(conv([Tm1 1],[Tm+Tap 1]),Kap*K1*Km);
HR12=tf(wnp/(2*titap),[1/(2*titap*wnp) 1 0]);
HR12=HR12*Hf11; %HR1"
Ho22=tf(wnp*wnp,[1 2*titap*wnp wnp*wnp]);
HOinchis3=minreal(feedback(Hd11,1));


figure
step(HOinchis1,'r')
hold on
step(HOinchis2,'g')
step(HOinchis3,'b')
legend("Ho1","Ho2","Ho3")
figure
lsim(HOinchis1,0:0.1:20,0:0.1:20,'r')
hold on
lsim(HOinchis2,0:0.1:20,0:0.1:20,'g')
lsim(HOinchis3,0:0.1:20,0:0.1:20,'b')

legend("Ho1","Ho2","Ho3")
%% 1.3 Calculul regulatorului cu metoda repartitiei pol-zero( Guilleman Truxel) corectat
close all;
estp=0;
sigma=0.05;
tr=1;
wb=12;
estv=0.05;
cv=1/estv;
tita=abs(log(sigma))/sqrt(log(sigma)^2+pi^2);
wn=2*tita*cv;
wb=wn*sqrt(1-2*tita^2+sqrt(2-4*tita^2+4*tita^4));
tr=4/(tita*wn);
tita2=0.6901;
wn2=wb/(sqrt(1-2*tita2^2+sqrt(2-4*tita2^2+4*tita2^4)));

Hf=tf(Kap*K1*Km,conv(conv([Tap 1],[Tm1 1]),[Tm 1]));
Ho2=tf(wn*wn,[1 2*tita*wn wn*wn]);

pc=4.923;
zc=4.91;

Hoc=minreal(tf([wn*wn*pc wn*wn*zc*pc],conv([1 2*tita*wn wn^2],[zc zc*pc])));
% figure
% step(Hoc)

HR2=zpk(minreal(Hoc/(1-Hoc)*(1/Hf)));
HR21=tf(0.02456*0.0262*conv([1 100 ],conv([1 16.66],conv([1 4.91],[1 1.785]))),[1 4.862 0]);
HR22=tf(9.3989*conv([0.5702 1],[0.2637 1]),[1 4.862 0]);
HR23=tf(0.2456*1.4132*conv([1 1.754],[1 3.792]),[1/4.862 0]);

Hd11=Hf*HR21;
Hd12=Hf*HR22;
Hd13=Hf*HR23;
Ho21=minreal(feedback(Hd11,1));
Ho22=minreal(feedback(Hd12,1));
Ho23=minreal(feedback(Hd13,1));

figure
step(Hoc,Ho21,Ho22,Ho23)
legend("Ho2","Ho21","Ho22","Ho23")
figure
lsim(Hoc,Ho21,Ho22,Ho23,0:0.1:20,0:0.1:20)
legend("Ho2","Ho21","Ho22","Ho23")


%%  Metode frecventiale
%P
close all
Tf=Tap +Tm1;
Kf=Kap*K1;
T=Tm/Km;
Hf=tf(Kf,[Tf*T T 0]);
bode(Hf)
sigma=0.05;
wf=1/Tf;
wt=5.81;

tita=abs(log(sigma))/sqrt(log(sigma)^2+pi^2);
A=1/(4*tita^2);
A=db(A);
F=-10.2;
FN=abs(F)-abs(A);
kp=db2mag(FN);
wn=2*tita*wt;
tr=4/(tita*wn);
Hfo=minreal(feedback(Hf*kp,1))
step(Hfo)

figure
lsim(Hfo,0:0.1:20,0:0.1:20)

%bode(Hf*kp)

%%
%PI
close all
Tf=Tap+Tm1;
Kf=Kap*K1;
T=Tm/Km;
Hf=tf(Kf,[Tf*T T 0]);
bode(Hf)
sigma=0.075;
wf=1/Tf;
wt=6.03;
cv1=db2mag(16);
cv2=db2mag(10);

tita=abs(log(sigma))/sqrt(log(sigma)^2+pi^2);
A=1/(4*tita^2);
A=db(A);
F=-10.2;
FN=abs(F)-abs(A);
kpi=db2mag(FN)*cv2/cv1;
wn=2*tita*wt;
tr=4/(tita*wn);

wz=0.1*wt;
wp=(cv1*wz)/10;
Tz=1/wz;
Tp=1/wp;
Hd=tf(kp*kpi*[Tz 1],[Tp 1]);
Hfo=minreal(feedback(Hf*Hd,1));

figure
step(Hfo)
figure
lsim(Hfo,0:0.1:20,0:0.1:20)

%%
%PD
close all
Tf=Tap +Tm1;
Kf=Kap*K1;
T=Tm/Km;
Hf=tf(Kf,[Tf*T T 0]);
bode(Hf)
sigma=0.01;
wf=1/Tf;
wt1=5.8;
tita=abs(log(sigma))/sqrt(log(sigma)^2+pi^2);
tr=2/(tita^2*wt1);
wt2=wt1*tr/0.3;

A=1/(4*tita^2);
A=db(A);
F=-10.2;
FN=abs(F)-abs(A);
kpd=db2mag(FN*wt2/wt1);
wn=2*tita*wt1;
Td=0.07;
TN=Tf*0.3/tr;
Hd=tf(kpd*kp*[Td 1],[TN 1]);
Hfo=minreal(feedback(Hf*Hd,1));

figure
step(Hfo)
figure
lsim(Hfo,0:0.1:20,0:0.1:20)

%%
%PD (PID)
close all
Tf=Tap +Tm1;
Kf=Kap*K1;
T=Tm/Km;
Hf=tf(Kf,[Tf*T T 0]);
wf=1/Tf;
wt1=5.8;
tita=0.63615;
tr=0.5;
wt2=wt1*tr/0.3;

A=1/(4*tita^2);
A=db(A);
F=-10.2;
FN=abs(F)-abs(A);
kpd=db2mag(FN*wt2/wt1);
Td=0.07;
TN=Tf*0.3/tr;
Hpd=tf(kpd*[Td 1],[TN 1]);

bode(Hpd)

%%
%PID
close all
Tf=Tap+Tm1;
Kf=Kap*K1;
T=Tm/Km;
Hf=tf(Kf,[Tf*T T 0]);
bode(Hf)
sigma=0.075;
tita=abs(log(sigma))/sqrt(log(sigma)^2+pi^2);
A=1/(4*tita^2);
tr=0.5;
wt1=5.8;
wt2=wt1*tr/0.3;
cv=14.6;
cvi=12;
wz=0.1*wt2;
wp=cv/cvi*wz;
A=db(A);
F=-10.2;
FN=abs(F)-abs(A);
kpd=db2mag(FN*wt2/wt1);
Tz=1/wz;
Tp=1/wp;
tauD=Tf;
TN=tauD*wt2/wt1;
kpi=cv/cvi;
HPID=tf(kpi*[Tz 1],[Tp 1]);

% HPID=tf(kpd*kpi*conv([tauD 1],[Tz 1]),conv([TN 1],[Tp 1]));
 Hfo=minreal(feedback(Hf*HPID*Hpd,1));

figure
step(Hfo)
figure
lsim(Hfo,0:0.1:20,0:0.1:20)

%% Metode frecventiale faza impusa
%PI
close all
Hf=tf(KceKu*Kc*Kphim,conv([Tv*Tc Tv+Tc 1],[Tphic 1]),'iodelay',Tauc);
bode(Hf)
wprim=0.00606;
phase=-180+15+50;
hf=db2mag(-13.4);
Kpi=1/hf;
taui=4/wprim;
HPI=tf(Kpi*[taui 1],[taui 0])
Hd=Hf*HPI

figure
bode(Hd)
%%
%PD
close all
Hf=tf(KceKu*Kc*Kphim,conv([Tv*Tc Tv+Tc 1],[Tphic 1]),'iodelay',Tauc);
bode(Hf)
w0=0.0148;
beta=0.12;
hf=db2mag(21.1);
Kpd=sqrt(beta)/hf;
tauD=1/w0*1/sqrt(beta);
tauN=beta*tauD;
HPD=tf(Kpd*[tauD 1],[tauN 1]);
Hd=Hf*HPD;
figure
bode(Hd)
%%
%PID
close all
Hf=tf(KceKu*Kc*Kphim,conv([Tv*Tc Tv+Tc 1],[Tphic 1]),'iodelay',Tauc);
bode(Hf)
beta=0.12;
hf1=db2mag(-13.4);
hf2=db2mag(21.1);
Kpid=0.228/hf2;
w0=0.0148;
T0=2*pi/w0;
tauI=1.2*T0;
tauD=0.5*T0;
HPID=tf(Kpid*conv([tauD 1],[tauI 1]),conv([beta*tauD 1],tauI));
Hd=Hf*HPID;
figure
bode(Hd)
