clear all;
%% definition of the paramaters
% GlobalRBe:volume-weighted average value of 10Be/9Be in global ocean basins
% qdel: the fraction of Be that survives the coastal scavenging
% Rkdtok: the ratio of Kd to k 
% Driv: The modern denudation rate (constant)
% NA:Avogadro constant
% CBecrust:the concentration of Be in the continental crust
% fmobile: the fraction of mobile (dissolved+ active) Be in the total 
% riverine Be flux
% A:the Mass number of Be
% Nsec:numbers of seconds per year
% Cl:coastalline length measured with a 100km-long ruler
% Drivall:The gloabl denudation rate (a variable)
%% the value of the constant 
NA=6.022*10^23;CBecrust=2.5*10^(-6);fmobile=0.2;A=9;
Nsec=3600*24*365;
Ariv=8.43;Aoc=36;F10Be=9.67*10^15;Alpha=145/182;Cl=255100;

%% calibration of model
Driv=182+182*0.2*(2*rand(1,1000000)-1);
GlobalRBe=10^-7*(1.035+2*0.06*(2*rand(1,1000000)-1));
 qdel=(Aoc/Ariv*F10Be/NA)./(GlobalRBe.*Driv.*CBecrust*10^6*fmobile/A-F10Be/NA);
 % calculation of Kd/k
 Rkdtok=((1./qdel)-1)./Driv/145*182/Ariv/10^7*Cl;
meanRkdtoK=mean(Rkdtok);
Drivall=1:0.1:10^4;
MRBeallocean1=(F10Be/NA*(Aoc/Ariv)+F10Be/NA./(1+meanRkdtoK.*Drivall*Alpha ...
*Ariv*10^7/Cl))./(Drivall*CBecrust*fmobile*10^6/A./(1+meanRkdtoK.*Drivall*Alpha ...
*Ariv*10^7/Cl));
WBeocean1=Drivall*CBecrust*fmobile*10^3;
Qdel1=1./(1+meanRkdtoK.*Drivall*Alpha*Ariv*10^7/Cl);
MRBeallocean2=(F10Be/NA*(Aoc/Ariv)+F10Be/NA./(1+meanRkdtoK.*Drivall.*(0.142.*log(Drivall)+0.0283)...
*Ariv*10^7/Cl))./(Drivall*CBecrust*fmobile*10^6/A./(1+meanRkdtoK.*Drivall.*(0.142.*log(Drivall)+0.0283)...
*Ariv*10^7/Cl));
Qdel2=1./(1+meanRkdtoK.*Drivall.*(0.142.*log(Drivall)+0.0283)*Ariv*10^7/Cl);

A(1)=182;B(1)=0.0000001035;
upli=10^-7*1.21*ones(1,99991);lowli=10^-7*0.63*ones(1,99991);
subplot(3,1,1);plot(Drivall,MRBeallocean1,'r'); hold on;plot(Drivall,MRBeallocean2);plot(A(1),B(1),'o');plot(Drivall,upli);plot(Drivall,lowli);hold off;
set(gca,'XScale','log');set(gca,'YScale','log');ylim([10^-8 10^-5]);xlabel('Denudation rate (t/km2/yr)');ylabel({'Seawater','10Be/9Be',});
subplot(3,1,2); plot(Drivall,WBeocean1);set(gca,'YScale','log');set(gca,'XScale','log');xlabel('Denudation rate (t/km2/yr)');ylabel({'Be weathering rate','kg/km2/yr',});ylim([0.0001 100]);
subplot(3,1,3); plot(Drivall,Qdel1,'r');hold on;plot(Drivall,Qdel2);hold off;set(gca,'YScale','log');set(gca,'XScale','log');xlabel('Denudation rate (t/km2/yr)');ylabel('qdel');

