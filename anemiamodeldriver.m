function anemiamodeldriver()

age = 30; %Age in years
gender = 0; %Gender, 0=male, 1=female
mass = 70; %Weight in kg
anemia = 0; %Is the patient anemic? 1 if yes, 0 if no

RQ=0.825;%reaction quotients for the general body and 
RQheart=0.7;

    if (40>=age)&&(age>0) %if/else statement that sets heart rate based on age, in bpm
        heartrate = 60;  
    else
        heartrate=70;
    end
    
bloodweight=0.07*mass; %sets mass of blood in the body based on percentage composition of blood and bodyweight


SV=mass/10; %stroke volume, blood pumped out of heart per beat in mL

bloodflow0 = SV*heartrate; %blood flow out of heart per minute

% cvector = concentration vector, [cE, cNa, cCa, cIron, cGlucose, cO2,
% cCO2, cHCO3], where c = "concentration of", and the letters are our
% components


    if gender == 0 %if loop sets volume percentage of red blood cells based on gender
   cE = 0.4345;
    elseif gender == 1
   cE = 0.402;  
    end
 
%set concentration of O2 and CO2 in mol/mL at BTP in venous blood
cO2 = 0.00000569745;
cCO2 = 0.0000235;
%concentration of bicarbonate in venous blood, in mol/L
cHCO3 = 0.00002033;
%concnetartion of glucose in mol/mL
cGlucose = 0.0000055;
%Na concentration in mol/mL
cNa = 0.000137;
%Concentration of Ca in mol/mL
cCa = 0.00000118;
%Concentration of iron in mol/mL
cIron = 0.00012592;

cvector0 = [cE cO2 cCO2 cHCO3 cGlucose cNa cCa cIron]; %concentration in the blood, in moles? mL/mL of blood

cEtrack=[cvector0(1)];
cO2track=[cvector0(2)];
cCO2track=[cvector0(3)];
cHCO3track=[cvector0(4)];
cGlucosetrack=[cvector0(5)];
cNatrack=[cvector0(6)];
cCatrack=[cvector0(7)];
cIrontrack=[cvector0(8)];

for loop=1:1
%Run initial venous blood through the heart
[bloodflow, cvector] = lungs(bloodflow0, cvector0);
%Run blood from lungs to the heart
[bloodflow, cvector] = heart(bloodflow, cvector);
%These three values are the blood flows that go from the heart to each
%organ
BFbraini=0.15*bloodflow;
BFliveri=0.25*bloodflow;
BFotherbloodi=bloodflow-BFbraini-BFliveri;
%Process each of these three blood flows in their respective organs
[BFbrainj, cvectorbrainj] = brain(BFbraini, cvector);
[BFliverj, cvectorliverj] = liver(BFliveri, cvector);
[BFotherbloodj, cvectorotherbloodj] = otherblood(BFotherbloodi, cvector); 
%Redirect blood from liver to other blood and mix the two
Mvectorotherbloodj=(cvectorotherbloodj*BFotherbloodj)+(cvectorliverj*BFliverj);
%Add the blood flow from other blood and liver to get new blood volume
BFotherbloodj=BFotherbloodj+BFliverj;
%Do a weighted average to get a new cvector
cvectorotherbloodj=Mvectorotherbloodj/BFotherbloodj;
%Send 25% of the blood flow to the kidneys
BFkidneyi=0.25*(BFbrainj+BFotherbloodj);
%Subtract blood sent to kidneys from other blood 
BFotherbloodj=BFotherbloodj-BFkidneyi;
%Have the kidneys process the blood they receive
[BFkidneyj, cvectorkidneyj] = kidney(BFkidneyi, cvectorotherbloodj);
%Send blood processed in kidneys back to other blood, pool with brain, create Mvector
Mvectorotherblood=(BFotherbloodj*cvectorotherblood)+(BFkidneyj*cvectorkidneyj)+(BFbrainj*cvectorbrainj);
%Reset the bloodflow back into the lungs to the blood flow we just computed
bloodflow0=BFotherbloodj+BFkidneyj+BFbrainj;
%Recompute cvector0 by dividing Mvectorotherblood by bloodflow0
cvector0=Mvectorotherblood/bloodflow0;

cEtrack=[cEtrack cvector0(1)];
cO2track=[cO2track cvector0(2)];
cCO2track=[cCO2track cvector0(3)];
cHCO3track=[cHCO3track cvector0(4)];
cGlucosetrack=[cGlucosetrack cvector0(5)];
cNatrack=[cNatrack cvector0(6)];
cCatrack=[cCatrack cvector0(7)];
cIrontrack=[cIrontrack cvector0(8)];

end

figure
plot(0:loop(end),cEtrack)
title('Erythrocyte Levels Over Time')
xlabel('Time in Minutes')
ylabel('Erythrocyte Concentration %Volume')

figure 
plot(0:loop(end),cO2track)
title('O2 Levels Over Time')
xlabel('Time in Minutes')
ylabel('O2 Concentration in mol/mL')

figure 
plot(0:loop(end),cCO2track)
title('CO2 Levels Over Time')
xlabel('Time in Minutes')
ylabel('CO2 Concentration in mol/mL')

figure 
plot(0:loop(end),cHCO3track)
title('HCO3 Levels Over Time')
xlabel('Time in Minutes')
ylabel('HCO3 Concentration in mol/mL')

figure 
plot(0:loop(end),cGlucosetrack)
title('Glucose Levels Over Time')
xlabel('Time in Minutes')
ylabel('Glucose Concentration in mol/mL')

figure
plot(0:loop(end),cNatrack)
title('Na Levels Over Time')
xlabel('Time in Minutes')
ylabel('Na Concentration in mol/mL')

figure
plot(0:loop(end),cCatrack)
title('Ca Levels Over Time')
xlabel('Time in Minutes')
ylabel('Ca Concentration in mol/mL')

figure
plot(0:loop(end),cIrontrack)
title('Iron Levels Over Time')
xlabel('Time in Minutes')
ylabel('Iron Concentration in mol/mL')





end

function [bloodflowj, cvectorj] = lungs(bloodflowi, cvectori) 

end

function [bloodflowj, cvectorj] = brain(bloodflowi, cvectori, mass)
%Blood volume conserved
bloodflowj=bloodflowi;

brainmass=0.02*mass*1000; %gives brain mass in grams as 2% body mass
Mvector=bloodflowi*cvectori;
%Create a null vector for cvectorj
cvectorj=zeros(1,8);
%Sets inflows equal to outflows for erythrocytes, Na, Ca, Iron
cvectorj(1)=cvectori(1);
cvectorj(6)=cvectori(6);
cvectorj(7)=cvectori(7);
cvectorj(8)=cvectori(8);

%Consume oxygen at 3.5 mL (in moles) per gram per minute
MO2j=Mvector(2)-0.0001375246*brainmass;
%Consume glucose at a rate of 120g/day
MGlucosej=Mvector(5)-(0.6660746/86400);%moles of glucose in - glucose requirements per day in moles/minutes in the day
%Calculate CO2 produced based on glucose consumed via respiration equation
MCO2j=Mvector(3)+6*(0.6660746/86400);
%Calculate HCO3 made when the CO2 is produced
rHCO3CO2=19.3/21.5; %ratio of bicarbonate to carbon dioxide in blood leaving lungs
HCO3j=Mvector(4)+rHCO3CO2*6*(0.6660746/86400);

%Calculate values above as concentrations and add them to cvectorj
cvectorj(2)=MO2j/bloodflowj;
cvectorj(5)=MGlucosej/bloodflowj;
cvectorj(3)=MCO2j/bloodflowj;
cvectorj(4)=HCO3j/bloodflowj;

end

function [bloodflowj, cvectorj] = heart(bloodflowi, cvectori)

end

function [cvectorj]=liver(G,W,A,cvectori)
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
end


function [bloodflowj, cvectorj] = kidneys(bloodflowi, cvectori, RQ)
T=86400;%Multiplier that scales up the time period of interest to one day (required for glucose equation)
Kidneymass=300; %mass of the kidneys combined in grams

%Initialize all intakes
Mvector=bloodflowi*cvectori;
bloodflowj=bloodflowi;

%Consume relevant compounds to create outlet values

MEj=Mvector(1);
MIronj=Mvector(8);

MO2j=0.85*MVector(2);
MNaj=(MNain-(Kidneymass/100)*((10*(25450*(Mvector(2)-MO2j)))-5))/1000;
MCaj=0.98*Mvector(7);
MGlucosej=Mvector(5)-(Mvector(5)*(0.226/((Mvector(5))*T)));
MCO2j=Mvector(3)+((Mvector(2)-MO2j)*RQ);
MHCO3j=Mvector(4)-0.15*(Mvector(4))+(.004/T);

%Recompute all concentrations using original mass of blood and new values

cIronj=MIronj/bloodflowj;
cNaj=MNaj/bloodflowj;
cO2j=MO2j/bloodflowj;
cCO2j=MCO2j/bloodflowj;
cEj=MEj/bloodflowj;
cGlucosej=MGlucosej/bloodflowj;
cHCO3j=MHCO3j/bloodflowj;
cCaj=MCaj/bloodflowj;

cvectorj=[cEj cO2j cCO2j cHCO3j cGlucosej cNaj cCaj cIronj];

end

function [bloodflowj, cvectorj] = otherblood(bloodflowi, cvectori)

end

