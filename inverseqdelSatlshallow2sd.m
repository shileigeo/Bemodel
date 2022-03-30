clear all;

%% Salt
%% definition of the paramaters
% F10Be:10Be flux per unit area NA:Avogadro constant
% Ariv:the drainage area of ocean basin
% Aoc:the surface area of ocean basin
% Driv: the denudation rate 
% advectionSatlFAntar: the advection flux from Antarctic ocean
% to south Atlantic per second
% qdel:the fraction of Be that survives the coastal scavenging
% CBecrust:the concentration of Be in the continental crust
% fmobile: the fraction of mobile (dissolved+ active) Be in the total 
% riverine Be flux
% A:the Mass number of Be
% SBeAntar: the Be concentration of the surface water in the Antarctic ocean
% SRBeAntar:the 10Be/9Be in the surface water in the Antarctic ocean
% Nsec:numbers of seconds per year
% CRBeSalt:caculated 10Be/9Be in the surface water in the South Atlantic ocean
% QdelSatl:inversely caculated qdel 
%% model
% the value of the constant 
F10Be=8*10^15;NA=6.022*10^23;Ariv=1;Aoc=4.1;
CBecrust=2.5*10^(-6);fmobile=0.2;A=9;SBeAntar=12.9;  
SRBeAntar=1.7*10^(-7);Nsec=3600*24*365;
%define the range of the paramaters
%         Driv ;advectionSatlFAntar;       qdel 
x1(1:3) = [49        16                    0.5       ]';
s1(1:3) = [49*0.2     3                    0.5       ]';
x1=x1';s1=s1';
t=1;
for i=1:1:2000000
    
% radom samplings within the assumed range of each paramaters   
k=rand(3,1)*2-1;kk=diag(k);ss1=kk*s1;X=x1+ss1;
%calcalate the uncorrected Beryllium istope ratios in seawater 
RBeSalt=(F10Be/NA*(Aoc/Ariv)+F10Be/NA*X(3))./...
    (X(1)*CBecrust*fmobile*1000000/A*X(3));
%correct the RBe for advection
favectedSaltl=X(2)*SBeAntar*Nsec*A/10^9./(X(2)*SBeAntar*Nsec*A/10^9+...
    X(1)*CBecrust*fmobile*X(3)*Ariv*10^7);
CRBeSalt=(10^8*RBeSalt*(1-favectedSaltl)+10^8*SRBeAntar*favectedSaltl);
if CRBeSalt<(12.7+2*2.1) && CRBeSalt>(12.7-2*2.1)
 XXSatl(t,:)=X;QdelSatl(t)=X(3);t=t+1;
end

end
save('inverseqdelSaltshallow2sd.mat', 'QdelSatl');