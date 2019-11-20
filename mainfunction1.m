function mainfunction1()

age = 30; %Age in years
gender = 0; %Gender, 0=male, 1=female
mass = 70; %Weight in kg
anemia = 0; %Is the patient anemic? 1 if yes, 0 if no

RQ=0.825;%reaction quotients for the general body and 
%RQheart=0.7;

    if (40>=age)&&(age>0) %if/else statement that sets heart rate based on age, in bpm
        heartrate = 60;  
    else
        heartrate=70;
    end

    
%Intake of various nutrients
carbs=130;
calciumintake=1000;
sodiumintake=500;

if gender==1
    ironintake=18;
else
    ironintake=8;
end

bloodweight=0.07*mass; %sets mass of blood in the body based on percentage composition of blood and bodyweight


%SV=mass; %stroke volume, blood pumped out of heart per beat in mL

bloodflow0 = 1000*bloodweight/1.06; %blood flow out of heart per minute in mL, divided by density of 1.06 g/mL (steady state value)

% cvector = concentration vector, [cE, cNa, cCa, cIron, cGlucose, cO2,
% cCO2, cHCO3], where c = "concentration of", and the letters are our
% components


if gender == 0 %if loop sets volume percentage of red blood cells based on gender
cE = 0.4345;
elseif gender == 1
cE = 0.402;  
end
 
%set concentration of O2 and CO2 in mol/mL at BTP in venous blood
cO2 = 0.00000785855;
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

for loop=1:1440
%Run initial venous blood through the heart
[bloodflow, cvector] = lungs(bloodflow0, cvector0);
%Run blood from lungs to the heart
[bloodflow, cvector1] = heart(bloodflow, cvector, mass);
%These three values are the blood flows that go from the heart to each
%organ





if (18<=age) && (age<50)
   LW=452+16.34*mass+11.85*age-166*gender;
else
   LW=1390+15.94*mass-12.86*age;
end

%Volumetric flow rate in and out of liver (mL/min)
%Vi=Vj=V=100mL/min per 100g liver
%Blood flow distribution will change during iron deficiency anemia but we
%don't have a function modeling this yet
V=LW;






BFbraini=0.15*bloodflow;
BFliveri=V;
BFotherbloodi=bloodflow-BFbraini-0.3*BFliveri;
%Process each of these three blood flows in their respective organs
[BFbrainj, cvectorbrainj] = brain(BFbraini, cvector1, mass);
[BFotherbloodj, cvectorotherbloodj] = otherblood(BFotherbloodi, cvector1, carbs, calciumintake, sodiumintake, ironintake); 

Mvectorotherbloodliver=0.7*V*cvectorotherbloodj;
Mvectorheart=0.3*V*cvector1;
cvectorliverin=(Mvectorotherbloodliver+Mvectorheart)/V;

[BFliverj, cvectorliverj] = liver(BFliveri,cvectorliverin,gender,mass,LW);
%Redirect blood from liver to other blood and mix the two
Mvectorotherbloodj=(cvectorotherbloodj*(BFotherbloodj-0.7*V))+(cvectorliverj*BFliverj);
%Add the blood flow from other blood and liver to get new blood volume
BFotherbloodj=BFotherbloodj-0.7*V+BFliverj;
%Do a weighted average to get a new cvector
cvectorotherbloodj=Mvectorotherbloodj/BFotherbloodj;
%Send 25% of the blood flow to the kidneys
BFkidneyi=0.25*(BFbrainj+BFotherbloodj);
%Subtract blood sent to kidneys from other blood 
BFotherbloodj=BFotherbloodj-BFkidneyi;
%Have the kidneys process the blood they receive
[BFkidneyj, cvectorkidneyj] = kidney(BFkidneyi, cvectorotherbloodj,RQ,mass);
%Send blood processed in kidneys back to other blood, pool with brain, create Mvector
Mvectorotherblood=(BFotherbloodj*cvectorotherbloodj)+(BFkidneyj*cvectorkidneyj)+(BFbrainj*cvectorbrainj);
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

A = [cvector; cvector1; cvectorbrainj; cvectorotherbloodj; cvectorliverj; cvectorkidneyj; cvector0];
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




function [bloodout, Cout] = lungs(vblood, Cvector)
    % This function finds the output terms for each component coming out of
    % the lungs.
    % Input: W = weight;
    %        Cvector - vector containing initial concentrations of each component
    %        entering heart
    % Output: Cout - vector containing final concentrations of each
    %         component leaving heart
    
    % Finding volumetric flow rate of blood (volumetric flow rate of blood
    % in = out)
%     mblood = 0.07*W; %in kg
%     pblood = 1.06; %in kg/L %this value is dependent on other factors - can we use it? should we use 1?
%     vblood = mblood/pblood; %in L/min
    
    Cout = [];
    % Finding volumetric flow rate out of oxygen
    CiO2 = 0.000001964637; %5mL/100mL becomes .000001964637 mol
    %CdeoxygenatedO2 = 0.16; %from graph and partial pressure of oxygen in entering deoxygenated blood being...
                            %40 mmHg (this might also depend on hemoglobin)
    nO2i = vblood*(Cvector(2)+CiO2); %O2 in, mL/min
    nO2cons = 0.05*nO2i;%.0002082515; %5.3 mL/min becomes .0002082515 mol/min
    nO2j = nO2i - nO2cons; %O2 out, mol/min
    Cout(2) = nO2j/vblood; %mol/mL
    
    % Finding volumetric flow rate out of carbon dioxide
    nCO2i = vblood*Cvector(3);
    nCO2cons = vblood*(.0000015716531); %vblood(mL/min) and concentration of blood consumed (mol/mL)
    %CjCO2 = 0.48; %48 mL/100 mL, should be dependent on hemoglobin, oxygen, things like that
    nCO2j = nCO2i - nCO2cons; %CO2 out, L/min
    Cout(3) = nCO2j/vblood;
    
    % Finding volumetric flow rate out of bicarbonate
    rHCO3CO2 = Cvector(4)/Cvector(3); %ratio of bicarbonate to carbon dioxide in blood leaving lungs
    nHCO3j = rHCO3CO2*nCO2j; %HCO3 out, mol/min
    Cout(4) = nHCO3j/vblood;
    
    % Finding volumetric flow rate out of calcium
    Cout(7) = Cvector(7); %calcium out, L/min
    
    % Finding volumetric flow rate out of iron
%      Cerythrocytes = 0.45; %45 mL/100 mL, this concentration changes depending on hemoglobin, but we need...
% %                           %to figure out this relationship
%      Chemoglobin = 0.335; %g/mL, this concentration also changes, but maybe only depending on demographics
% %                          %and anemia
%     Mhemoglobin = 65000; %g/mol
%     MFe = 55.845; %g/mol
%     CFe = (MFe*Chemoglobin*Cerythrocytes)/Mhemoglobin;
%     vFej = vblood*CFe; %iron out, L/min

    Cout(8) = Cvector(8); %iron concentration doesn't change in lungs?
    % Finding volumetric flow rate out of sodium
    Cout(6) = Cvector(6); %sodium out, L/min
    
    % Finding volumetric flow rate out of erythrocytes
    %vErythrocytesj = vblood*Cerythrocytes; %erythrocytes out, L/min, once again this concentration changes
    Cout(1) = Cvector(1);
    
    % Finding volumetric flow rate out of glucose
%     Chemoglobinblood = Chemoglobin*Cerythrocytes; %concentration of hemoglobin in blood, g/mL
%     Chemoglobinbloodmol = Chemoglobinblood/Mhemoglobin; %concentration of hemoglobin in blood, mol/mL
%     CO2bloodmol = 4*Chemoglobinbloodmol; %concentration of oxygen in blood dependent on this value of...
%                                          %hemoglobin, mol/mL
    %nO2cons = CO2bloodmol*vblood*1000; %molar flow rate of oxygen consumed, mol/min - can we somehow...
                                       %convert this into mL/min to get a
                                       %more accurate value of oxygen
                                       %consumed in the oxygen accounting
                                       %equation?
    nglucosecons = 6*nO2cons; %molar flow rate of glucose consumed, mol/min
%     Mglucose = 180.18; %g/mol
%     mglucosecons = nglucosecons*Mglucose; %g/min
%     pglucose = 1560; %g/L, is this what other people are using?
%     vglucosecons = mglucosecons/pglucose; %L/min
    nGlucosei = Cvector(5)*vblood;
    nGlucosej = nGlucosei - nglucosecons; %molar flow rate of glucose out, mol/min
    Cout(5) = nGlucosej/vblood;
    bloodout = vblood;
end
    
    %*thoughts I had while coding that we should consider: putting in
    %checks to make sure outflows that go directly to the next organ are
    %the same as those inflows - if they're not the same, find the problem
    %and/or change it to be the same
    % make sure everyone codes volumetric flow rate of blood and density of
    % blood in the same way





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

%Consume oxygen at 3.5 mL (in moles) per 100 gram per minute
MO2j=Mvector(2)-0.0001375246*(brainmass/100);
%Consume glucose at a rate of 120g/day
MGlucosej=Mvector(5)-(0.6660746/1440);%moles of glucose in - glucose requirements per day in moles/minutes in the day
%Calculate CO2 produced based on glucose consumed via respiration equation
MCO2j=Mvector(3)+6*(0.6660746/1440);
%Calculate HCO3 made when the CO2 is produced
rHCO3CO2=19.3/21.5; %ratio of bicarbonate to carbon dioxide in blood leaving lungs
HCO3j=Mvector(4)+rHCO3CO2*6*(0.6660746/1440);

%Calculate values above as concentrations and add them to cvectorj
cvectorj(2)=MO2j/bloodflowj;
cvectorj(5)=MGlucosej/bloodflowj;
cvectorj(3)=MCO2j/bloodflowj;
cvectorj(4)=HCO3j/bloodflowj;

end






%I don't know about modeling different conditions based on flow at a
%certain time, would it be better if we used t (for time, i.e. what step we
%are at)
%Also if carbon dioxide concentration changes, then I am pretty sure
%bicarbonate concentration will also change (by how much I'm not sure, I'm
%tired)
function [outflow, Cout] = heart(flow,Cvector,weight)
    % This function finds the output terms for each component coming out of
    % the heart.
    % Input: Cvector - vector containing initial concentrations of each component
    %        entering heart
    %        Flow - number that represents volumetric? flow of first inlet
    %        stream from lungs to heart
    % Output: Cout - vector containing final concentrations of each
    %         component leaving heart
    
    %how do I account for total blood/blood flow
    %so the Heart has multiple flows to it, but the heart will actually consume only once -- the first time because we are assuming that everything is happening in one heartbeat? 
Cout = [];
    %if flow == Flow1LungsToHeart
    
        % components that don't change:
        
        % Finding concentration out of glucose
Cout(5)=Cvector(5);
% Finding concentration out of erythrocytes
Cout(1)=Cvector(1);
% Finding concentration out of sodium
Cout(6)=Cvector(6);
% Finding concentration out of iron
Cout(8)=Cvector(8);
% Finding concentration out of calcium
Cout(7)=Cvector(7);
% Finding concentration out of bicarbonate, however I think
% bicarbonate will change because the concentration of carbon
% dioxide changes
%Cout(4)=Cvector(4);

heartmass = weight*.003;
nO2cons = 32.5*(heartmass/300)*0.0000392927; %mL/min*mol/mL=mol/min
%vblood = (.07*weight)/1.06; %L/min, if the total volumetric blood flow leaving the lungs is the same as volumetric blood flow entering the heart
Cout(2) = (Cvector(2)*flow - nO2cons)/flow;
CO2gen = .7*nO2cons; %using respiratory quotient in the heart
Cout(3) = (Cvector(3)*flow + CO2gen)/flow;
rHCO3CO2 = Cvector(4)/Cvector(3); %ratio of bicarbonate to carbon dioxide in blood leaving lungs
nHCO3j = rHCO3CO2*CO2gen; %HCO3 out, mol/min
Cout(4) = (Cvector(4)*flow + nHCO3j)/flow;

    %components that do change
    % Finding volumetric flow rate out of oxygen
    % Finding volumetric flow rate out of carbon dioxide
    %heartRQ=0.7; %we are assuming that heart muscle cell metabolism comes soley from lipids
    %heartRQ=CO2produced/O2consumed;
    % O2 consumption: 30 - 35 mL/min per 300 g
    %based off of RQ find the CO2produced
    %vO2j=vO2-O2consumed but do units match
    %vCO2j=vCO2 + CO2produced
       
%     else 
%         % no components will change?
%         %Concentration2HeartToOther
%         %Concentration2HeartToBrain
%         %Concentration2HeartToLiver
%         %Concentration3BrainToHeart
%         %Concentration4LiverToHeart
%         %Concentration7OtherToHeart
%         %Concentration8HeartToLungs
%         
%         Cout(5)=Cvector(5);
%         Cout(1)=Cvector(1);
%         Cout(6)=Cvector(6);
%         Cout(8)=Cvector(8);
%         Cout(7)=Cvector(7);
%         Cout(4)=Cvector(4);
%         Cout(3)=Cvector(3);
%         Cout(2)=Cvector(2);
%     end
outflow = flow;
end




function [bloodflowj, cvectorout]=liver(V,cvectori,G,mass,LW)
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
    p=1.2e10; %number of red blood cells per mL of PURE red blood cells
else
    p=1.2e10;
end
    %Calculate volume of the body
    Vblood=0.07*mass*1000/1.056;
    %Calculate the number of RBCs entering liver based on concentration,
    %flow rate, and pure RBC density
    RBCi=p*V*cvectori(1);
    %Consume some based on body size
    RBCcons=5.55555e-6*cvectori(1)*Vblood;%this had vickis original p here, confirm
    %Find pure volume of RBCs we now have based on number and density 
    VRBCj=(RBCi-RBCcons)/p;
    %Calculate volume of plasma that entered
    Vplasma=V*(1-cvectori(1));
    %Calculate new concentrations based on volumeRBCs/total volume
    cvectorout(1)=VRBCj/(VRBCj+Vplasma);
    %Bloodflow out = blood cells out + plasma volume
    bloodflowj=Vplasma+VRBCj;

%Calculates concentration flow rate of O2 out 
%No generation of O2
%nO2i=molar flow rate of O2/min in
%nO2j=molar flow rate of O2/min out
%nO2cons=molar consumption rate of O2/min
nO2i= cvectori(2)*V; %0.0006286842*(LW/100);
nO2cons=0.0002389*(LW/100);
nO2j=nO2i-nO2cons;
cvectorout(2)=nO2j/V;

%Calculates concentration flow rate of CO2 out
%nCO2i=molar flow rate of CO2/min in
%nCO2j=molar flow rate of CO2/min out
%nCO2gen=molar flow rate of CO2/min generated
%RQ=Respiratory Quotient=mol CO2 produced/mol O2 consumed
RQ=0.825;
nCO2gen=RQ*nO2cons;
nCO2i=cvectori(3)*V;
nCO2j=nCO2i+nCO2gen;
cvectorout(3)=nCO2j/V;

%Calculates the required intake of protein in grams/min given a
%person's demographic information
%dprotein=daily protein intake in g/day
%protein=protein intake in g/min
dprotein=0.8*mass;
protein=dprotein/1440;

%Calculates the concentration flow rate of Bicarbonate out 
%Assume no generation of bicarbonate
%%Also liver just takes in bicarbonate from the heart and stomach
%1mole bicarbonate ions consumed by liver per 100g protein/day
rHCO3CO2 = cvectori(4)/cvectori(3); %ratio of bicarbonate to carbon dioxide in blood leaving lungs
nHCO3gen = rHCO3CO2*nCO2gen; %HCO3 out, mol/min
nHCO3i=cvectori(4)*V;
nHCO3j=nHCO3i-(protein/100)+nHCO3gen;
cvectorout(4)=nHCO3j/V;

%Calculates the concentration flow rate of Na2+ out 
%Assume no generation of Na2+, only in, out, and consumption
%nNai=molar flow rate of sodium in 
%nNacons=molar flow rate of sodium consumed given that 10% of bile salts
%are consumed, approx 600mL are produced per day, and there are 135mM Na+
%ions in bile.
%nNaj=molar flow rate of sodium out
nNacons=9.375e-9;
nNai=cvectori(6)*V;
nNaj=nNai-nNacons;
cvectorout(6)=nNaj/V;

%Calculates the concentration flow rate of Ca2+ out
%Assume no generation nor consumption.
%cCai=concentration flow rate of Ca2+ in mol/Lmin in
%cCaj=concentration flow rate of Ca2+ in mol/Lmin out
cvectorout(7)=cvectori(7);
%Calculates the concentration flow rate of Fe out
%Assume no generation since

%Still need to do Glucose in Liver
%Placeholder that sets in=out
cvectorout(5)=cvectori(5);
cvectorout(8)=cvectori(8);
end


function [bloodflowj, cvectorj] = kidney(bloodflowi, cvectori, RQ, mass)
T=1440;%Multiplier that scales up the time period of interest to one day (required for glucose equation)
Kidneymass=(300/70000)*mass; %mass of the kidneys combined in grams

%Initialize all intakes
Mvector=bloodflowi*cvectori;
bloodflowj=bloodflowi;

%Consume relevant compounds to create outlet values

MEj=Mvector(1);
MIronj=Mvector(8);

MO2j=0.85*Mvector(2);
% Find mol/min O2 consumed
%MO2cons = Mvector(2) - MO2j;
% Convert to mL/min per 100g
%MO2cons = (MO2cons/0.0000392927)/(Kidneymass/100);
% Find mmol/min per 100g Na consumed
%MNacons = (MO2cons - 0.5)/0.1;
% Convert to mol/min
%MNacons = (MNacons/1000)*(Kidneymass/100);
MNacons = (100/1440)/1000;
MNaj = Mvector(6) - MNacons;
%MNaj=(Mvector(6)-(Kidneymass/100)*((10*(25450*(Mvector(2)-MO2j)))-5))/1000;
MCaj=0.98*Mvector(7);
MGlucosej=Mvector(5)-(Mvector(5)*(0.226/((Mvector(5))*T)));
MCO2j=Mvector(3)+((Mvector(2)-MO2j)*RQ);
rHCO3O2=Mvector(4)/Mvector(3);
MHCO3gen=rHCO3O2*((Mvector(2)-MO2j)*RQ);
MHCO3j=Mvector(4)-0.15*(Mvector(4))+(.004/T)+MHCO3gen;

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

function [outflow, Cout] = otherblood(flow, Cin, carbs, calciumintake, sodiumintake, ironintake)
Cout = [];
Cout(1) = Cin(1);
%CO2cons = 
Cout(2) = Cin(2);%-CO2cons;
Cout(3) = Cin(3);
Cout(4) = Cin(4);
Cout(5) = (flow*Cin(5) + (.8*carbs)/(180.156*1440))/flow; %130g of carbs, 80 percent is directly translated to glucose, 20% is fructose and galactose which goes to the liver and may be converted to glucose at a later stage so another term is incoming
Cout(6) = (flow*Cin(6) + (sodiumintake/(1000*22.99*1440)))/flow; %500 mg of sodium needed for vital functions so start with sodium intake being 500 mg
Cout(7) = (flow*Cin(7) + (calciumintake*.26)/(1000*40.08*1440))/flow; %calcium intake should be around 1000 mg
Cout(8) = (flow*Cin(8) + (ironintake*.18)/(1000*55.845*1440)-1/(1440*1000*55.845))/flow;%iron intake will be approximately 8 mg for males and 18 mg for females
outflow = flow;
end