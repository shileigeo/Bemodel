clear all;

%% Natl£» advection from Salt
%% definition of the paramaters
% F10Be:10Be flux per unit area NA:Avogadro constant
% Ariv:the drainage area of ocean basin
% Aoc:the surface area of ocean basin
% Driv: the denudation rate 
% advectionNatlFSalt: the advection flux from South Atlantic ocean
% to North Atlantic Ocean per second
% qdel:the fraction of Be that survives the coastal scavenging
% CBecrust:the concentration of Be in the continental crust
% fmobile: the fraction of mobile (dissolved+ active) Be in the total 
% riverine Be flux
% A:the Mass number of Be
% SBeSatl: the Be concentration of the surface water in the South Atlantic ocean
% SRBeSatl:the 10Be/9Be in the surface water in the South Atlantic ocean
% Nsec:numbers of seconds per year
% CRBeNatl:caculated 10Be/9Be in the surface water in the North atlantic ocean
% QdelNatl:inversely caculated qdel 
%% model
% the value of the constant 
F10Be=12*10^15;NA=6.022*10^23;Ariv=2.81;Aoc=4.88;
CBecrust=2.5*10^(-6);fmobile=0.2;A=9;
Nsec=3600*24*365;
%define the range of the paramaters x:mean, s:error
%         Driv ;advectionNatlFSalt;        qdel       SRBeSatl  SBeSatl 
x1(1:5) = [141         16                   0.5         12.7      17]';
s1(1:5) = [141*0.2     2                   0.5          2*2.1     9.3]';
x1=x1';s1=s1';t=1;
for i=1:1:2000000
    
% radom samplings within the assumed range of each paramaters   
k=rand(5,1)*2-1;kk=diag(k);ss1=kk*s1;X=x1+ss1;
%calcalate the uncorrected Beryllium istope ratios in seawater 
RBeNalt=(F10Be/NA*(Aoc/Ariv)+F10Be/NA*X(3))./...
       (X(1)*CBecrust*fmobile*1000000/A*X(3));
%correct the RBe for advection
favectedNatl=X(2)*X(5)*Nsec*A/10^9./(X(2)*X(5)*Nsec*A/10^9+...
X(1)*CBecrust*fmobile*X(3)*Ariv*10^7);
CRBeNat=(10^8*RBeNalt*(1-favectedNatl)+X(4)*favectedNatl);
if CRBeNat<(4.58+2*1.09) && CRBeNat>(4.58-2*1.09)
 XXNatl(t,:)=X;
 QdelNatl(t)=X(3);
 t=t+1;
end
end
save('inverseqdelNatlshallow2sd.mat', 'QdelNatl');
