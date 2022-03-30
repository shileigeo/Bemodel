clear all;

%% Spaci
%% definition of the paramaters
% F10Be:10Be flux per unit area NA:Avogadro constant
% Ariv:the drainage area of ocean basin
% Aoc:the surface area of ocean basin
% Driv: the denudation rate 
% advectionSpaciFAnta: the advection flux from Antarctic ocean
% to South Pacific per second
% qdel:the fraction of Be that survives the coastal scavenging
% CBecrust:the concentration of Be in the continental crust
% fmobile: the fraction of mobile (dissolved+ active) Be in the total 
% riverine Be flux
% A:the Mass number of Be
% SBeAntar: the Be concentration of the surface water in the Antarctic ocean
% SRBeAntar:the 10Be/9Be in the surface water in the Antarctic ocean
% Nsec:numbers of seconds per year
% CRBeSpaci:caculated 10Be/9Be in the surface water in the South Pacific ocean
% QdelSpaci:inversely caculated qdel 
%% model
% the value of the constant 
F10Be=10*10^15;NA=6.022*10^23;Ariv=0.14;Aoc=9.4;
CBecrust=2.5*10^(-6);fmobile=0.2;A=9;SBeAntar=12.9;  
SRBeAntar=1.7*10^(-7);Nsec=3600*24*365;
%define the range of the paramaters
%         Driv ;advectionSpaciFAntar;        qdel 
x1(1:3) = [400       19                   0.5       ]';
s1(1:3) = [400*0.2     5                    0.5       ]';
x1=x1';s1=s1';
t=1;
for i=1:1:2000000
    
% radom samplings within the assumed range of each paramaters   
k=rand(3,1)*2-1;kk=diag(k);ss1=kk*s1;X=x1+ss1;
%calcalate the uncorrected Beryllium istope ratios in seawater 
RBeSpaci=(F10Be/NA*Aoc/Ariv+F10Be/NA*X(3))./...
    (X(1)*CBecrust*fmobile*1000000/A*X(3));
%correct the RBe for advection
favectedSpacil=X(2)*SBeAntar*Nsec*A/10^9./(X(2)*SBeAntar*Nsec*A/10^9+...
    X(1)*CBecrust*fmobile*X(3)*Ariv*10^7);
CRBeSpaci=(10^8*RBeSpaci*(1-favectedSpacil)+10^8*SRBeAntar*favectedSpacil);
if CRBeSpaci<(26.2+5.2*2) && CRBeSpaci>(26.2-5.2*2)
 XXSpaci(t,:)=X;
 QdelSpaci(t)=X(3);
 t=t+1;
end
end
