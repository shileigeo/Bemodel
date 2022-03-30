clear all;
%% Mete;advection from North Atlantic
%% definition of the paramaters
% F10Be:10Be flux per unit area NA:Avogadro constant
% Ariv:the drainage area of ocean basin
% Aoc:the surface area of ocean basin
% Driv: the denudation rate 
% advectionNatlFSalt: the advection flux from Nouth Atlantic ocean
% to North Atlantic Ocean per second
% qdel:the fraction of Be that survives the coastal scavenging
% CBecrust:the concentration of Be in the continental crust
% fmobile: the fraction of mobile (dissolved+ active) Be in the total 
% riverine Be flux
% A:the Mass number of Be
% SBeNatl: the Be concentration of the surface water in the North Atlantic ocean
% SRBeNatl:the 10Be/9Be in the surface water in the North Atlantic ocean
% Nsec:numbers of seconds per year
% CRBeMete:caculated 10Be/9Be in the surface water in the Mediterranean ocean
% Qdelmete:inversely caculated qdel 
%% model
% the value of the constant 
F10Be=10*10^15;NA=6.022*10^23;Ariv=0.394;Aoc=0.25;
CBecrust=2.5*10^(-6);fmobile=0.2;A=9;
Nsec=3600*24*365;
%define the range of the paramaters
%         Driv ; advectionMetelFNalt ;qdelmete   SRBeNatl   SBeNalt
x1(1:5) = [160          0.78            0.5       4.58            28.9]';
s1(1:5) = [160*0.2      0.46            0.5       2*1.09           6.9 ]';
x1=x1';s1=s1';t=1;
for i=1:1:2000000
 % radom samplings within the assumed range of each paramaters   
k=rand(5,1)*2-1;kk=diag(k);ss1=kk*s1;X=x1+ss1;
%calcalate the uncorrected Beryllium istope ratios in seawater 
RBeMete=(F10Be/NA*(Aoc/Ariv)+F10Be/NA*X(3))./...
       (X(1)*CBecrust*fmobile*1000000/A*X(3));
%correct the RBe for advection
favectedMete1=X(2)*X(5)*Nsec*A/10^9/(X(2)*X(5)*Nsec*A/10^9+...
    (X(1)*CBecrust*fmobile*X(3))*Ariv*10^7);
CRbeMete=(10^8*RBeMete*(1-favectedMete1)+X(4)*favectedMete1);
CRbe(i)=CRbeMete;
% Save the results that can solve the isotopic mass balance equations
if CRbeMete<(1.06+2*0.24) && CRbeMete>(1.06-2*0.24) 
 XXMete(t,:)=X; Qdelmete(t)=X(3);t=t+1;
end
end

save('inverseqdelMeteshallow2sd.mat', 'Qdelmete');
