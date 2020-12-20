clear all;
%% definition of the paramaters
% GlobalRBe:volume-weighted average value of 10Be/9Be in global ocean basins
% qdel: the fraction of Be that survives the coastal scavenging
% Rkdtok: the ratio of Kd to k 
% Driv: The modern denudation rate 
% NA:Avogadro constant
% CBecrust:the concentration of Be in the continental crust
% fmobile: the fraction of mobile (dissolved+ active) Be in the total 
% riverine Be flux
% A:the Mass number of Be
% Nsec:numbers of seconds per year
% Cl:coastalline length measured with a 100km-long ruler
%% the value of the constant
NA=6.022*10^23;CBecrust=2.5*10^(-6);fmobile=0.2;A=9;
Nsec=3600*24*365;
Ariv=8.43;Aoc=36;F10Be=9.67*10^15;Alpha=145/182;Cl=255100;

%% calibration of model
Driv=182+182*0.2*(2*rand(1,1000000)-1);
GlobalRBe=10^-7*(1.035+2*0.06*(2*rand(1,1000000)-1));
 qdel=(Aoc/Ariv*F10Be/NA)./(GlobalRBe.*Driv.*CBecrust*10^6*fmobile/A...
     -F10Be/NA);
 % calculation of Kd/k
 Rkdtok=((1./qdel)-1)./Driv/Alpha/Ariv/10^7*Cl;
MRkdtok=(max(Rkdtok)+min(Rkdtok))/2;
WRkdtok=max(Rkdtok)-min(Rkdtok);
load('RBenormalized.mat');

%% Model run under 'present-day' Scenario
SmoothedRBeN = smooth(RBeN,10);
GlobalRBemax=10^-7*(1.035+0.06);
GlobalRBemin=10^-7*(1.035-0.06);
RkdtoK=mean(Rkdtok); 
for i=1:138
   RBenormalizedmax=SmoothedRBeN(i)+0.30;
   RBenormalizedmin=SmoothedRBeN(i)-0.30;
 t=1;solution=[];
while (size(solution))<100000

Denudationrate=750.5+749.5*(2*rand-1);
 ModeledRBe=(F10Be/(NA)*(Aoc/Ariv)+F10Be/(NA)/(1+RkdtoK*Denudationrate*(0.142*log(Denudationrate)+0.0283)...
*Ariv*10^7/Cl))/(Denudationrate*CBecrust*fmobile*1000000/A./(1+RkdtoK*Denudationrate*(0.142*log(Denudationrate)+0.0283)...
*Ariv*10^7/Cl));
if  ModeledRBe<RBenormalizedmax*GlobalRBemax && ModeledRBe>RBenormalizedmin*GlobalRBemin
    
    solution(t)=Denudationrate;
    t=t+1;
    i
    t
end
end
sumsolution(i,:)=solution;

end


for i=1:138
maxsolution(i)=max(sumsolution(i,:));
minsolution(i)=min(sumsolution(i,:));


end

% plot(Age,maxsolution,'r',Age,minsolution,'b')
Age1=Age;
maxsolution1=maxsolution;
for i=139:276
    Age1(i)=Age1(277-i);
    maxsolution1(i)=minsolution(277-i);
end
FSmoothedRBeN=SmoothedRBeN-0.30;

for i=139:276
    FSmoothedRBeN(i)=SmoothedRBeN(277-i)+0.30;
end

%% Model run under 'reworked-clay' Scenario
for i=1:138
   RBenormalizedmax2=SmoothedRBeN(i)+0.30;
   RBenormalizedmin2=SmoothedRBeN(i)-0.30;
 t=1;solution2=[];
while (size(solution2))<100000

Denudationrate2=750.5+749.5*(2*rand-1);
 ModeledRBe2=(F10Be/(NA)*(Aoc/Ariv)+F10Be/(NA)/(1+RkdtoK*Denudationrate2*Alpha ...
*Ariv*10^7/Cl))/(Denudationrate2*CBecrust*fmobile*1000000/A/(1+RkdtoK*Denudationrate2*Alpha ...
*Ariv*10^7/Cl));
if  ModeledRBe2<RBenormalizedmax2*GlobalRBemax && ModeledRBe2>RBenormalizedmin2*GlobalRBemin
    
    solution2(t)=Denudationrate2;
    t=t+1;
    i
    t
end
end
sumsolution2(i,:)=solution2;

end


for i=1:138
maxsolution2(i)=max(sumsolution2(i,:));
minsolution2(i)=min(sumsolution2(i,:));


end

% plot(Age,maxsolution,'r',Age,minsolution,'b')
Age2=Age;
maxsolution22=maxsolution2;
for i=139:276
    Age2(i)=Age2(277-i);
    maxsolution22(i)=minsolution2(277-i);
end
%%  plot the results

ConstantW=182*2.5*10^-3*0.2*ones(1,138);
ConstantD=182*ones(1,138);
subplot(3,2,1);plot(Age,ConstantD);hold on; fill(Age1,maxsolution1,'r');xlim([0 12.33]);ylim([0 1500]);hold off; xlabel('Age(Ma)');ylabel({'Continental','denudaton rate','t/km2/yr',});
subplot(3,2,2);plot(Age,ConstantW);hold on;fill(Age1,maxsolution1*2.5*10^-3*0.2,'r');xlim([0 12.33]);ylim([0 0.75]);xlabel('Age(Ma)');ylabel({'Be weathering rate','kg/km2/yr',});
subplot(3,2,3);plot(Age,ConstantD);hold on; fill(Age2,maxsolution22,'r');xlim([0 12.33]);ylim([0 1500]);hold off;xlabel('Age(Ma)');ylabel({'Continental','denudaton rate','t/km2/yr',});
subplot(3,2,4);plot(Age,ConstantW);hold on;fill(Age2,maxsolution22*2.5*10^-3*0.2,'r');xlim([0 12.33]);ylim([0 0.75]);xlabel('Age(Ma)');ylabel({'Be weathering rate','kg/km2/yr',});
subplot(3,2,5);errorbar(Age,RBeN,eRBeN,'o'); hold on;fill(Age1,FSmoothedRBeN,'r'); plot(Age,SmoothedRBeN);
plot(Age,SmoothedRBeN-0.30);plot(Age,SmoothedRBeN+0.30);xlim([0 12.33]);xlabel('Age(Ma)');ylabel({'Normalized','10Be/9Be','in seawater'});




