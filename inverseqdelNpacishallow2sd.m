clear all;
%% Npaci;advection from Spaci
%% definition of the paramaters
% F10Be:10Be flux per unit area NA:Avogadro constant
% Ariv:the drainage area of ocean basin
% Aoc:the surface area of ocean basin
% Driv: the denudation rate 
% advectionNpaciFSpaci: the advection flux from South Pacific ocean
% to North Pacific Ocean per second
% qdel:the fraction of Be that survives the coastal scavenging
% CBecrust:the concentration of Be in the continental crust
% fmobile: the fraction of mobile (dissolved+ active) Be in the total 
% riverine Be flux
% A:the Mass number of Be
% SBeSpaci: the Be concentration of the surface water in the South Pacific ocean
% SRBeSpaci:the 10Be/9Be in the surface water in the South Pacific  ocean
% Nsec:numbers of seconds per year
% CRBeNpaci:caculated 10Be/9Be in the surface water in the North Pacifc ocean
% QdelNpaci:inversely caculated qdel 
%% model
% the value of the constant 
F10Be=12*10^15;NA=6.022*10^23;Ariv=1.26;Aoc=8.3;
CBecrust=2.5*10^(-6);fmobile=0.2;A=9;
Nsec=3600*24*365;
%define the range of the paramaters
%         Driv ; advectionNpaciFSpaci ;   qdel      SBeSpaci   SRBeSpaci 
x1(1:5) = [316         2.25                0.5       7.6         26.2    ]';
s1(1:5) = [316*0.2     2.25                0.5       2.2         2*5.2   ]';
x1=x1';s1=s1';
t=1;
for i=1:1:2000000
    
% radom samplings within the assumed range of each paramaters   
k=rand(5,1)*2-1;kk=diag(k);ss1=kk*s1;X=x1+ss1;
%calcalate the uncorrected Beryllium istope ratios in seawater 
RBeNpaci=(F10Be/NA*(Aoc/Ariv)+F10Be/NA*X(3))./...
    (X(1)*CBecrust*fmobile*1000000/A*X(3));
%correct the RBe for advection
favectedNpaci=(X(2)*X(4)*Nsec*A/10^9)./(X(2)*X(4)*Nsec*A/10^9+...
X(1)*CBecrust*fmobile*X(3)*Ariv*10^7);
CRBeNpaci=(10^8*RBeNpaci*(1-favectedNpaci)+X(5)*favectedNpaci);
% Save the results that can solve the isotopic mass balance equations
if CRBeNpaci<(13.8+2*3.8) && CRBeNpaci>(13.8-2*3.8)
 XXNpaci(t,:)=X;
 QdelNpaci(t)=X(3);
 t=t+1;
end
end
save('inverseqdelNpacishallow2sd.mat', 'QdelNpaci');
