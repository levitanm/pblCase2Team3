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

bloodflow = SV*heartrate; %blood flow out of heart per minute

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

cvector = [cE cO2 cCO2 cHCO3 cGlucose cNa cCa cIron]; %concentration in the blood, in moles? mL/mL of blood

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

function [bloodflowj, cvectorj] = liver(bloodflowi, cvectori)

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

