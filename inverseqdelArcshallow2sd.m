clear all;
 
%% Arctic;advection from Berning strait, Barent sea opening and Fram Strait
%% Mete;advection from North Atlantic
%% definition of the paramaters
% F10Be:10Be flux per unit area NA:Avogadro constant
% Ariv:the drainage area of ocean basin
% Aoc:the surface area of ocean basin
% Driv: the denudation rate 
% AdFBering:the advection flux from Bering Strait
% AdFBerents:the advection flux from Barents Sea
% AdFFram:the advection flux from Fram Strait
% qdel:the fraction of Be that survives the coastal scavenging
% CBecrust:the concentration of Be in the continental crust
% fmobile: the fraction of mobile (dissolved+ reactive) Be in the total 
% riverine Be flux
% A:the Mass number of Be
% SBeBarent: the Be concentration of the surface water in the Barents Sea 
% SRBeBarent:the 10Be/9Be in the surface water in the Barents sea 
% SBeBering: the Be concentration of the surface water in the Bering Strait 
% SRBeBering:the 10Be/9Be in the surface water in the Bering Strait
% SBeFram:the Be concentration of the surface water in the Fram Strait
% SRBeFram:the 10Be/9Be in the surface water in the Fram Strait
% Nsec:numbers of seconds per year
% CRBeArc:caculated 10Be/9Be in the surface water in the Arctic Ocean
% QdelArc:inversely caculated qdel 
%% model
% the value of the constant
F10Be=4*10^15;NA=6.022*10^23;Ariv=1.44;Aoc=1.62;
CBecrust=2.5*10^(-6);fmobile=0.2;A=9;
Nsec=3600*24*365;
%define the range of the paramaters
%         ;Driv        AdFBering;      AdFBerents;       AdFFram;   qdel   SRBeBarents  SBeBarents   SRBeBering   SBeBering     SRBeFram  SBeFram        
x1(1:11) = [37          0.8              2              1.8       0.5       8.05        12.1         11.3        18.7            2.81     23.2 ]';
s1(1:11) = [37*0.2      0.2           0.00001           0.1       0.5      2*2.03         3          2*3.5         5             2*0.1     2.1]';
x1=x1';s1=s1';t=1;
for i=1:1:2000000
% radom samplings within the assumed range of each paramaters   
k=rand(11,1)*2-1;kk=diag(k);ss1=kk*s1;X=x1+ss1;
%calcalate the uncorrected Beryllium istope ratios in seawater 
RBeArc=(F10Be/NA*(Aoc/Ariv)+F10Be/NA*X(5))./...
       (X(1)*CBecrust*fmobile*1000000/A*X(5));
%correct the RBe for advection
favectedArc1=(X(2)*X(9)*Nsec*A/10^9)/(X(2)*X(9)*Nsec*A/10^9+...
    X(3)*X(7)*Nsec*A/10^9+...
       X(4)*X(11)*Nsec*A/10^9+...
    X(1)*CBecrust*fmobile*X(5)*Ariv*10^7);
 favectedArc2=(X(3)*X(7)*Nsec*A/10^9)/(X(2)*X(9)*Nsec*A/10^9+...
    X(3)*X(7)*Nsec*A/10^9+...
       X(4)*X(11)*Nsec*A/10^9+...
    X(1)*CBecrust*fmobile*X(5)*Ariv*10^7);
  favectedArc3=(X(4)*X(11)*Nsec*A/10^9)/(X(2)*X(9)*Nsec*A/10^9+...
    X(3)*X(7)*Nsec*A/10^9+...
    X(4)*X(11)*Nsec*A/10^9+...
    X(1)*CBecrust*fmobile*X(5)*Ariv*10^7);
 CRBeArc=(10^8*RBeArc*(1-favectedArc1-favectedArc2-favectedArc3)+favectedArc1*X(8)+favectedArc2*X(6)+...
 favectedArc3*X(10));
% Save the results that can solve the isotopic mass balance equations
if CRBeArc<(4.72+2*3.12) && CRBeArc>(4.72-2*3.12)
 XXArc(t,:)=X;QdelArc(t)=X(5);t=t+1;
end
i
t
end

