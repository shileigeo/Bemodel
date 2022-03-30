clear all;
%% definition of the paramaters 
% GlobalRBe:volume-weighted average value of 10Be/9Be in global ocean basins
% qdel: the fraction of Be that survives the coastal scavenging
% Rkdtok: the ratio of Kd to k 
% Driv***: The modern denudation rate 
% Ariv***:the drainage area of ocean basins
% Aoc***:the surface area of ocean basins
% SBe***:The Be concentration in the surface ocean water
% SRBe***:The 10Be/9Be in the surface ocean water
% F10Be***:10Be flux per unit area
% Alpha***:the ratio of physical erosion to denudation rate
% Cl***:coastalline length measured with a 100km-long ruler
% Here *** could be one of the follwing: all,Arc, Mete, Natl, Satl,Spaci. Npaci, Barent,Bering
% Fram or Antar
% They represent the global oceans, the Arctic, the Mediterranean, the North Atlantic, 
% the South Atlantic, the North Pacific, the South Pacific,Barents Sea, 
% Bering Strait, Fram Strait,and the Antarctic Ocean respectively
% NA:Avogadro constant
% CBecrust:the concentration of Be in the continental crust
% fmobile: the fraction of mobile (dissolved+ active) Be in the total 
% riverine Be flux
% A:the Mass number of Be
% Nsec:numbers of seconds per year
%% the value of the constant
NA=6.022*10^23;CBecrust=2.5*10^(-6);fmobile=0.2;A=9;
Nsec=3600*24*365;
Arivall=8.43;ArivArc=1.44;ArivMete=0.394;ArivNatl=2.81;ArivSatl=1;ArivSpaci=0.140;
ArivNpaci=1.26;
Aocall=36;AocArc=1.62;AocMete=0.25;AocNatl=4.88;AocSatl=4.1;AocSpaci=9.4;
AocNpaci=8.3;
F10Beall=9.67*10^15;F10BeArc=4*10^15;F10BeMete=10*10^15;F10BeNatl=12*10^15;
F10BeSatl=8*10^15;F10BeSpaci=10*10^15;F10BeNpaci=12*10^15;
Alphaall=145/182;AlphaArc=16/37;AlphaMete=138/160;AlphaNatl=104/141;AlphaSatl=32/49;
AlphaSpaci=329/400;AlphaNpaci=264/316;
Clall=255100;ClArc=30700;ClMete=14300;ClNatl=67400;ClSatl=11600;ClSpaci=18800;
ClNpaci=50400;
%% calibration of model
Drivall=182+182*0.2*(2*rand(1,1000000)-1);
GlobalRBe=10^-7*(1.035+2*0.06*(2*rand(1,1000000)-1));
 qdel=(Aocall/Arivall*F10Beall/NA)./(GlobalRBe.*Drivall.*CBecrust*10^6*fmobile/A...
     -F10Beall/NA);
 % calculation of Kd/k
 Rkdtok=((1./qdel)-1)./Drivall/Alphaall/Arivall/10^7*Clall;

%% Be concentration and isotopic ratios in surface oceans
SRBeNatl=4.58+2*1.09*(2*rand(1,1000000)-1);SBeNatl=28.9+6.9*(2*rand(1,1000000)-1);
SRBeSatl=12.7+2*2.1*(2*rand(1,1000000)-1);SBeSatl=17+9.3*(2*rand(1,1000000)-1);
SRBeNpaci=13.8+2*3.8*(2*rand(1,1000000)-1);SBeNpaci=9.6+5.3*(2*rand(1,1000000)-1);
SRBeSpaci=26.2+2*5.2*(2*rand(1,1000000)-1);SBeSpaci=7.6+2.2*(2*rand(1,1000000)-1);
SRBeMet=1.06+2*0.24*(2*rand(1,1000000)-1);SBeMet=61.7+6.4*(2*rand(1,1000000)-1);
SRBeArc=4.72+2*3.12*(2*rand(1,1000000)-1);SBeArc=29.4+20.1*(2*rand(1,1000000)-1);
SRBeBarent=8.05+2*2.03*(2*rand(1,1000000)-1);SBeBarent=12.1+3*(2*rand(1,1000000)-1);
SRBeBering=11.3+2*3.5*(2*rand(1,1000000)-1);SBeBering=18.7+5*(2*rand(1,1000000)-1);
SRBeFram=2.81+2*0.1*(2*rand(1,1000000)-1);SBeFram=23.2+2.1*(2*rand(1,1000000)-1);
SRBeAntar=17; SBeAntar=12.9; 

%% Mete;advection from North Atlantic
advectionMeteFNalt=0.78+0.46*(2*rand(1,1000000)-1);
 DrivMete=160+160*0.2*(2*rand(1,1000000)-1);
 qdel1=1./(1+Rkdtok.*DrivMete*AlphaMete.*ArivMete*10^7/ClMete);
RBeMete=(F10BeMete/NA*(AocMete/ArivMete)+F10BeMete/NA.*qdel1)./(DrivMete*CBecrust*fmobile*1000000/A.*qdel1);
%correct the RBe for advection
favectedMete1=advectionMeteFNalt.*SBeNatl*Nsec*A/10^9./(advectionMeteFNalt.*SBeNatl*Nsec*A/10^9+...
    (DrivMete*CBecrust*fmobile.*qdel1)*ArivMete*10^7);
CRbeMete=RBeMete.*(1-favectedMete1)+10^-8*SRBeNatl.*favectedMete1;
%calculate the mean and error
mRbeMete=mean(RBeMete);
eRbeMete=std(RBeMete);
mCRbeMete=mean(CRbeMete);
eCRbeMete=std(CRbeMete);

%% Arctic;advection from Berning strait, Barent sea opening and Fram Strait
advectionArcFBeringStait=0.8+0.2*(2*rand(1,1000000)-1);
advectionArcFBerentsSeaOpening=2;
advectionArcFFramStrait=1.8+0.1*(2*rand(1,1000000)-1);
DrivArc=37+37*0.2*(2*rand(1,1000000)-1);
qdel2=1./(1+Rkdtok.*DrivArc*AlphaArc.*ArivArc*10^7/ClArc);
RBeArc=(F10BeArc/NA*(AocArc/ArivArc)+F10BeArc/NA.*qdel2)./(DrivArc*CBecrust*fmobile*1000000/A.*qdel2);
%correct the RBe for advection
favectedArc1=(advectionArcFBeringStait.*SBeBering*Nsec*A/10^9)./(advectionArcFBeringStait.*SBeBering*Nsec*A/10^9+...
    advectionArcFBerentsSeaOpening*SBeBarent*Nsec*A/10^9+...
    advectionArcFFramStrait.*SBeFram*Nsec*A/10^9+...
     DrivArc*CBecrust*fmobile.*qdel2*ArivArc*10^7);
 favectedArc2=(advectionArcFBerentsSeaOpening*SBeBarent*Nsec*A/10^9)./(advectionArcFBeringStait.*SBeBering*Nsec*A/10^9+...
    advectionArcFBerentsSeaOpening*SBeBarent*Nsec*A/10^9+...
    advectionArcFFramStrait.*SBeFram*Nsec*A/10^9+...
     DrivArc*CBecrust*fmobile.*qdel2*ArivArc*10^7);
  favectedArc3=(advectionArcFFramStrait.*SBeFram*Nsec*A/10^9)./(advectionArcFBeringStait.*SBeBering*Nsec*A/10^9+...
    advectionArcFBerentsSeaOpening*SBeBarent*Nsec*A/10^9+...
    advectionArcFFramStrait.*SBeFram*Nsec*A/10^9+...
     DrivArc*CBecrust*fmobile.*qdel2*ArivArc*10^7);
 CRBeArc=RBeArc.*(1-favectedArc1-favectedArc2-favectedArc3)+10^-8*favectedArc1.*SRBeBering+10^-8*favectedArc2.*SRBeBarent+...
 10^-8*favectedArc3.*SRBeFram;
%calculate the mean and error
mRbeArc=mean(RBeArc);
eRbeArc=std(RBeArc);
mCRbeArc=mean(CRBeArc);
eCRbeArc=std(CRBeArc);


%% Npaci;advection from Spaci
advectionNpaciFSpaci=2.25+2.25*(2*rand(1,1000000)-1);
DrivNpaci=316+316*0.2*(2*rand(1,1000000)-1);
qdel3=1./(1+Rkdtok.*DrivNpaci*AlphaNpaci.*ArivNpaci*10^7/ClNpaci);
RBeNpaci=(F10BeNpaci/NA*(AocNpaci/ArivNpaci)+F10BeNpaci/NA.*qdel3)./(DrivNpaci*CBecrust*fmobile*1000000/A.*qdel3);
%correct the RBe for advection
favectedNpaci=(advectionNpaciFSpaci.*SBeSpaci*Nsec*A/10^9)./(advectionNpaciFSpaci.*SBeSpaci*Nsec*A/10^9+...
DrivNpaci*CBecrust*fmobile.*qdel3*ArivNpaci*10^7);
CRBeNpaci=RBeNpaci.*(1-favectedNpaci)+10^-8*SRBeSpaci.*favectedNpaci;
%calculate the mean and error
mRbeNpaci=mean(RBeNpaci);
eRbeNpaci=std(RBeNpaci);
mCRbeNpaci=mean(CRBeNpaci);
eCRbeNpaci=std(CRBeNpaci);

%% Natl£» advection from Salt
advectionNatlFSalt=16+2*(2*rand(1,1000000)-1);
DrivNatl=141+141*0.2*(2*rand(1,1000000)-1);
qdel4=1./(1+Rkdtok.*DrivNatl*AlphaNatl.*ArivNatl*10^7/ClNatl);
RBeNatl=(F10BeNatl/NA*(AocNatl/ArivNatl)+F10BeNatl/NA.*qdel4)./(DrivNatl*CBecrust*fmobile*1000000/A.*qdel4);
%correct the RBe for advection
favectedNatl=advectionNatlFSalt.*SBeSatl*Nsec*A/10^9./(advectionNatlFSalt.*SBeSatl*Nsec*A/10^9+...
DrivNatl*CBecrust*fmobile.*qdel4*ArivNatl*10^7);
CRBeNatl=RBeNatl.*(1-favectedNatl)+10^-8*SRBeSatl.*favectedNatl;
%calculate the mean and error
mRbeNatl=mean(RBeNatl);
eRbeNatl=std(RBeNatl);
mCRbeNatl=mean(CRBeNatl);
eCRbeNatl=std(CRBeNatl);
%% Salt
advectionSatlFAntar=16+3*(2*rand(1,1000000)-1);
DrivSatl=49+49*0.2*(2*rand(1,1000000)-1);
qdel5=1./(1+Rkdtok.*DrivSatl*AlphaSatl.*ArivSatl*10^7/ClSatl);
RBeSatl=(F10BeSatl/NA*(AocSatl/ArivSatl)+F10BeSatl/NA.*qdel5)./(DrivSatl*CBecrust*fmobile*1000000/A.*qdel5);
%correct the RBe for advection
favectedSaltl=advectionSatlFAntar.*SBeAntar*Nsec*A/10^9./(advectionSatlFAntar.*SBeAntar*Nsec*A/10^9+...
    DrivSatl*CBecrust*fmobile.*qdel5*ArivSatl*10^7);
CRBeSatl=RBeSatl.*(1-favectedSaltl)+10^-8*SRBeAntar.*favectedSaltl;
%calculate the mean and error
mRbeSatl=mean(RBeSatl);
eRbeSatl=std(RBeSatl);
mCRbeSatl=mean(CRBeSatl);
eCRbeSatl=std(CRBeSatl);
%% Spaci
advectionSpaciFAntar=19+5*(2*rand(1,1000000)-1);

DrivSpaci=400+400*0.2*(2*rand(1,1000000)-1);
qdel6=1./(1+Rkdtok.*DrivSpaci*AlphaSpaci*ArivSpaci*10^7/ClSpaci);
RBeSpaci=(F10BeSpaci/NA*(AocSpaci/ArivSpaci)+F10BeSpaci/NA.*qdel6)./(DrivSpaci*CBecrust*fmobile*1000000/A.*qdel6);
%correct the RBe for advection
favectedSpaci=advectionSpaciFAntar.*SBeAntar*Nsec*A/10^9./(advectionSpaciFAntar.*SBeAntar*Nsec*A/10^9+...
    DrivSpaci*CBecrust*fmobile.*qdel6*ArivSpaci*10^7);
CRBeSpaci=RBeSpaci.*(1-favectedSpaci)+10^-8*SRBeAntar.*favectedSpaci;
%calculate the mean and error
mRbeSpaci=mean(RBeSpaci);
eRbeSpaci=std(RBeSpaci);
mCRbeSpaci=mean(CRBeSpaci);
eCRbeSpaci=std(CRBeSpaci);

