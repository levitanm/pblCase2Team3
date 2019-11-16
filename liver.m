function cvectorout=liver(G,W,A,cvectorin)
%This function will deliver output volumetric flow rates for the concentration of
%the 8 components (mol/mLmin) out of the liver
%Input: G=Gender (1 if female, 0 if male)
%W=weight of person (kg); A=age of person (years)
%cvector=[cE,cO2,cCO2,cHCO3,cGlucose,cNa,cCa,cFe)
%Outputs are all in concentration flow rates (mol/mLmin) out

%Initializes cvectorout
cvectorout=[0,0,0,0,0,0,0,0];

%Wet weight of liver (g) due to person's demographics
%LW=Liver weight (g)
%make this into a seperate function if can although i really believe liver
%needs demographic info in

if (18<=A) && (A<50)
   LW=452+16.34*W+11.85*A-166*G;
else
   LW=1390+15.94*W-12.86*A;
end

%Volumetric flow rate in and out of liver (mL/min)
%Vi=Vj=V=100mL/min per 100g liver
%Blood flow distribution will change during iron deficiency anemia but we
%don't have a function modeling this yet
V=LW;

%Calculates concentration flow rate of red blood cells out 
%p=blood density in amount per mL
%Vbody=Volume of body in mL. Calculating assuming 7% of total body weight
%is blood, then divided by density of blood to find volume of blood
%RBCj=raw numerical flow rate of red blood cells out/min
%RBCi=raw numerical flow rate of RBC in/min
%nRBCj=molar flow rate of RBC out/min
%RBCcons=consumption rate of RBC/min 
%No generation of RBC in liver
if G==1
    p=5.4e9;
else
    p=4.8e9;
end
    Vbody=0.07*W*1000/1.056;
    RBCi=p*100*LW/100;
    RBCcons=5.55555e-6*p*Vbody;
    RBCj=RBCi-RBCcons;
    nRBCj=RBCj/(6.022e23);
    cvectorout(1)=nRBCj/V;

%Calculates concentration flow rate of O2 out 
%No generation of O2
%nO2i=molar flow rate of O2/min in
%nO2j=molar flow rate of O2/min out
%nO2cons=molar consumption rate of O2/min
nO2i=0.0006286842*LW/100;
nO2cons=0.0002389*LW/100;
nO2j=nO2i-nO2cons;
cvectorout(2)=nO2j/V;

%Calculates concentration flow rate of CO2 out
%nCO2i=molar flow rate of CO2/min in
%nCO2j=molar flow rate of CO2/min out
%nCO2gen=molar flow rate of CO2/min generated
%RQ=Respiratory Quotient=mol CO2 produced/mol O2 consumed
RQ=0.825;
nCO2gen=RQ*nO2cons;
nCO2i=cvectorin(3)*V;
nCO2j=nCO2i+nCO2gen;
cvectorout(3)=nCO2j/V;

%Calculates the required intake of protein in grams/min given a
%person's demographic information
%dprotein=daily protein intake in g/day
%protein=protein intake in g/min
dprotein=0.8*W;
protein=dprotein/1440;

%Calculates the concentration flow rate of Bicarbonate out 
%Assume no generation of bicarbonate
%%Also liver just takes in bicarbonate from the heart and stomach
%1mole bicarbonate ions consumed by liver per 100g protein/day=0.06944g/min
nHCO3i=cvectorin(4)*V;
nHCO3j=nHCO3i-protein/0.06944;
cvectorout(4)=nHCO3j/V;

%Calculates the concentration flow rate of Na2+ out 
%Assume no generation of Na2+, only in, out, and consumption
%nNai=molar flow rate of sodium in 
%nNacons=molar flow rate of sodium consumed given that 10% of bile salts
%are consumed, approx 600mL are produced per day, and there are 135mM Na+
%ions in bile.
%nNaj=molar flow rate of sodium out
nNacons=9.375e-9;
nNai=cvectorin(6)*V;
nNaj=nNai-nNacons;
cvectorout(6)=nNaj/V;

%Calculates the concentration flow rate of Ca2+ out
%Assume no generation nor consumption.
%cCai=concentration flow rate of Ca2+ in mol/Lmin in
%cCaj=concentration flow rate of Ca2+ in mol/Lmin out
cvectorout(7)=cvectorin(7);

%Calculates the concentration flow rate of Fe out
%Assume no generation since

%Still need to do Iron and Glucose in Liver




