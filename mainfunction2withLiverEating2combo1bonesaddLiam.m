function mainfunction2withLiverEating2combo1bonesaddLiam()

age = 50;%30; %Age in years
gender = 0; %Gender, 0=male, 1=female
mass = 70; %Weight in kg
anemia = 1; %Is the patient anemic? 1 if yes, 0 if no
bleed = 0.01; %bleeding rate in mL/min, if anemic
FerritinStores = 0;%0.0179; %Stored moles of iron in liver as ferritin

RQ=0.825;%reaction quotients for the general body and
%RQheart=0.7;

if (40>=age) %if/else statement that sets heart rate based on age, in bpm
    baseheartrate = 60;
else 
    baseheartrate= 70;
end

%Intake of various nutrients
carbs=130;
calciumintake=1000;
sodiumintake=500;

if gender==1
    ironintake=8;
else
    ironintake=18;
end

basebloodweight=0.07*mass; %sets mass of blood in the body based on percentage composition of blood and bodyweight
%SV=mass; %stroke volume, blood pumped out of heart per beat in mL
bloodweight=basebloodweight;
bloodflow0 = 1000*bloodweight/1.06; %blood flow out of heart per minute in mL, divided by density of 1.06 g/mL (steady state value)



% cvector = concentration vector, [cE, cNa, cCa, cIron, cGlucose, cO2,
% cCO2, cHCO3], where c = "concentration of", and the letters are our
% components

if gender == 0 %if loop sets volume percentage of red blood cells based on gender
    cE0 = 0.4345;
elseif gender == 1
    cE0 = 0.402;
end

basehemoglobin=34.52243959*basebloodweight/10*cE0; %hemoglobin, in grams, that should be present in the body
hemoglobin=basehemoglobin;
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
cCa = 0.00000248;
%Concentration of iron in mol/mL
cIron = 7.2231e-06;

cvector0 = [cE0 cO2 cCO2 cHCO3 cGlucose cNa cCa cIron]; %concentration in the blood, in moles? mL/mL of blood

cEtrack=[cvector0(1)];
cO2track=[cvector0(2)];
cCO2track=[cvector0(3)];
cHCO3track=[cvector0(4)];
cGlucosetrack=[cvector0(5)];
cNatrack=[cvector0(6)];
cCatrack=[cvector0(7)];
cIrontrack=[cvector0(8)];
heartratetrack = [baseheartrate];
bloodflowvec=[bloodflow0];
iter = 0;
%while cvector0(1) >=.375
for loop=1:1440
    heartrate=(baseheartrate*cE0*basebloodweight)/(cvector0(1)*bloodweight);
    bloodflow0 = (heartrate/baseheartrate)*bloodflow0;
    % all blood per min for 60 beats per min -> all blood per 60 beats
    %(heartrate/60)*bloodflow0
    
    
    %Run initial venous blood through the heart
    [bloodflow, cvector, Ci] = lungs(bloodflow0, cvector0, anemia, basehemoglobin, hemoglobin);
    %Run blood from lungs to the heart
    [bloodflow, cvector1] = heart(bloodflow, cvector, mass, Ci, basehemoglobin, hemoglobin);
    %These three values are the blood flows that go from the heart to each
    %organ
    
    
    
    %Volumetric flow rate in and out of liver (mL/min)
    %Vi=Vj=V=100mL/min per 100g liver
    %Blood flow distribution will change during iron deficiency anemia but we
    %don't have a function modeling this yet
    
    if (18<=age) && (age<50)
        V=452+16.34*mass+11.85*age-166*gender;
    else
        V=1390+15.94*mass-12.86*age;
    end
    
    BFbraini=0.15*bloodflow;
    BFliveri=V;
    BFotherbloodi=bloodflow-BFbraini-0.3*BFliveri;
    %Process each of these three blood flows in their respective organs
    [BFbrainj, cvectorbrainj] = brain(BFbraini, cvector1, mass, Ci, basehemoglobin, hemoglobin);
    
    [BFotherbloodj, cvectorotherbloodj, hemoglobinout] = otherblood(BFotherbloodi, cvector1, carbs, calciumintake, sodiumintake, ironintake, Ci, RQ, anemia, bleed, basehemoglobin, hemoglobin);
    
    Mvectorotherbloodliver=0.7*V*cvectorotherbloodj;
    Mvectorheart=0.3*V*cvector1;
    cvectorliverin=(Mvectorotherbloodliver+Mvectorheart)/V;
    
    [BFliverj, cvectorliverj, FerritinStores1] = liver(BFliveri,cvectorliverin,gender,mass,Ci,FerritinStores, basehemoglobin, hemoglobin);
    %Redirect blood from liver to other blood and mix the two
    FerritinStores=FerritinStores1;
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
    
    [BFkidneyj, cvectorkidneyj] = kidney(BFkidneyi, cvectorotherbloodj,RQ,mass, Ci, V, basehemoglobin, hemoglobin, anemia);
    %Send blood processed in kidneys back to other blood, pool with brain, create Mvector
    Mvectorotherblood=(BFotherbloodj*cvectorotherbloodj)+(BFkidneyj*cvectorkidneyj)+(BFbrainj*cvectorbrainj);
    %Reset the bloodflow back into the lungs to the blood flow we just computed
    bloodflow0=BFotherbloodj+BFkidneyj+BFbrainj;
    bloodflowvec=[bloodflowvec bloodflow0];
    bloodflow0 = 1000*bloodweight/1.06;
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
    heartratetrack = [heartratetrack heartrate];
    A = [cvector; cvector1; cvectorbrainj; cvectorotherbloodj; cvectorliverj; cvectorkidneyj; cvector0];
    iter = iter + 1;
    BP1=(38.178*log(bloodweight/basebloodweight)+100)*(1.8886*(cvector0(1)/cE0)-0.8886);
    hemoglobin=hemoglobin-hemoglobinout;
end
% figure
% plot(0:loop(end),heartratetrack)

figure
plot(0:loop(end),cEtrack,'k','LineWidth',2, 'DisplayName', 'Erythrocytes Concentration')
title('Erythrocyte Levels Over Time')
xlabel('Time in Minutes')
ylabel('Erythrocyte Concentration %Volume')
v1 = [0 .43; loop(end) .43; loop(end) .57; 0 .57;];
f1 = [1 2 3 4];
patch('Faces',f1,'Vertices',v1,'FaceColor','green','FaceAlpha',.1, 'LineStyle',':', 'DisplayName','Normal Physiological Values Range');
legend


figure
plot(0:loop(end),cO2track,'k','LineWidth',2, 'DisplayName', 'Oxygen Concentration')
title('O2 Levels Over Time')
xlabel('Time in Minutes')
ylabel('O2 Concentration in mol/mL')
v1 = [0 7.12e-6; loop(end) 7.12e-6; loop(end) 9.24e-6; 0 9.24e-6;];
f1 = [1 2 3 4];
patch('Faces',f1,'Vertices',v1,'FaceColor','green','FaceAlpha',.1,'LineStyle',':', 'DisplayName','Normal Physiological Values Range');
legend

figure
plot(0:loop(end),cCO2track,'k','LineWidth',2, 'DisplayName', 'Carbon Dioxide Concentration')
title('CO2 Levels Over Time')
xlabel('Time in Minutes')
ylabel('CO2 Concentration in mol/mL')
v1 = [0 2.15e-5; loop(end) 2.15e-5; loop(end) 2.33e-5; 0 2.33e-5;];
f1 = [1 2 3 4];
patch('Faces',f1,'Vertices',v1,'FaceColor','green','FaceAlpha',.1, 'LineStyle',':', 'DisplayName','Normal Physiological Values Range');
legend

figure
plot(0:loop(end),cHCO3track,'k','LineWidth',2, 'DisplayName', 'Bicarbonate Concentration')
title('HCO3 Levels Over Time')
xlabel('Time in Minutes')
ylabel('HCO3 Concentration in mol/mL')
v1 = [0 1.93e-5; loop(end) 1.93e-5; loop(end) 2.03e-5; 0 2.03e-5;];
f1 = [1 2 3 4];
patch('Faces',f1,'Vertices',v1,'FaceColor','green','FaceAlpha',.1, 'LineStyle',':', 'DisplayName','Normal Physiological Values Range');
legend

figure
plot(0:loop(end),cGlucosetrack,'k','LineWidth',2, 'DisplayName', 'Glucose Concentration')
title('Glucose Levels Over Time')
xlabel('Time in Minutes')
ylabel('Glucose Concentration in mol/mL')
v1 = [0 3.8855e-6; loop(end) 3.8855e-6; loop(end) 9.991e-6; 0 9.991e-6;];
f1 = [1 2 3 4];
patch('Faces',f1,'Vertices',v1,'FaceColor','green','FaceAlpha',.1, 'LineStyle',':', 'DisplayName','Normal Physiological Values Range');
legend

figure
plot(0:loop(end),cNatrack,'k','LineWidth',2, 'DisplayName', 'Sodium Concentration')
title('Na Levels Over Time')
xlabel('Time in Minutes')
ylabel('Na Concentration in mol/mL')
v1 = [0 1.35e-4; loop(end) 1.35e-4; loop(end) 1.45e-4; 0 1.45e-4;];
f1 = [1 2 3 4];
patch('Faces',f1,'Vertices',v1,'FaceColor','green','FaceAlpha',.1, 'LineStyle',':', 'DisplayName','Normal Physiological Values Range');
legend

figure
plot(0:loop(end),cCatrack,'k','LineWidth',2, 'DisplayName', 'Calcium Concentration')
title('Ca Levels Over Time')
xlabel('Time in Minutes')
ylabel('Ca Concentration in mol/mL')
v1 = [0 2.2e-6; loop(end) 2.2e-6; loop(end) 2.7e-6; 0 2.7e-6;];
f1 = [1 2 3 4];
patch('Faces',f1,'Vertices',v1,'FaceColor','green','FaceAlpha',.1, 'LineStyle',':', 'DisplayName','Normal Physiological Values Range');
legend

figure
plot(0:loop(end),cIrontrack,'k','LineWidth',2, 'DisplayName', 'Iron Concentration')
title('Iron Levels Over Time')
xlabel('Time in Minutes')
ylabel('Iron Concentration in mol/mL')
end

function [bloodout, Cout, Ci] = lungs(vblood, Cvector, anemia, basehemoglobin, hemoglobin)










Cout = [];
% anemia needs to be an input for the lungs function now
Co = 0.21e-4; %concentration of oxygen in alveoli (21% b/c air - it's 0.21e-4 b/c otherwise the order of
              %magnitude was way off (not sure if this is valid but it works))
P = 7.39e-3*60; %cm/min, membrane permeability constant for alveoli
if anemia == 1
    A = 1000000; %cm^2, surface area of alveoli under extreme conditions
else
    A = 700000; %cm^2, surface area of alveoli under normal conditions
end
Ci = (Cvector(2)*0.08206*310.15)/(0.0526*1000); %mL/mL, concentration of oxygen in blood, convert mol to mL
                                                %with partial pressure of
                                                %oxygen going into lungs
                                                %(40 mmHg)
ViO2 = P*A*(Co - Ci); %mL/min
ViO2 = ViO2*(hemoglobin/basehemoglobin);
%convert ViO2 to mol/min
niO2 = (ViO2*1)/(0.08206*310.15); %mol/min; convert mL to mol with pressure being atmospheric pressure
niO2 = niO2/15; %divide molar flow rate by respiratory rate
C = niO2/(vblood*1.5); %mol/mL, divide by conversion factor of 1.5

Ci = C*vblood*(hemoglobin/basehemoglobin); %concentration of oxygen from which other organs consume

nO2i = vblood*Cvector(2) + Ci;
nO2cons = 0.05*Ci/(hemoglobin/basehemoglobin);%we changed this when we added mayas O2 code 
%nO2cons = 0.05*C*vblood; %we removed this when we added mayas O2 code
nO2j = nO2i - nO2cons; %oxygen out, mol/min
Cout(2) = nO2j/vblood; %mol/mL













% Cout = [];
% % Finding volumetric flow rate out of oxygen
% CiO2 = 0.000001964637; %5mL/100mL becomes .000001964637 mol
% %CdeoxygenatedO2 = 0.16; %from graph and partial pressure of oxygen in entering deoxygenated blood being...
% Ci=CiO2*vblood;                       %40 mmHg (this might also depend on hemoglobin)
% nO2i = vblood*(Cvector(2)+CiO2); %O2 in, mL/min
% nO2cons = 0.05*CiO2*vblood;%.0002082515; %5.3 mL/min becomes .0002082515 mol/min
% nO2j = nO2i - nO2cons; %O2 out, mol/min
% Cout(2) = nO2j/vblood; %mol/mL













% Finding volumetric flow rate out of carbon dioxide
nCO2i = vblood*Cvector(3);
nCO2cons = nCO2i*(1-(48/52));%vblood*(.0000015716531); %vblood(mL/min) and concentration of blood consumed (mol/mL)
%CjCO2 = 0.48; %48 mL/100 mL, should be dependent on hemoglobin, oxygen, things like that
nCO2j = nCO2i - nCO2cons; %CO2 out, L/min
Cout(3) = nCO2j/vblood;


% Finding volumetric flow rate out of bicarbonate
rHCO3CO2 = Cvector(4)/Cvector(3); %ratio of bicarbonate to carbon dioxide in blood leaving lungs
nHCO3j = rHCO3CO2*nCO2j; %HCO3 out, mol/min
Cout(4) = nHCO3j/vblood;

% Finding volumetric flow rate out of calcium
Cout(7) = Cvector(7); %calcium out, L/min

Cout(8) = Cvector(8); %iron concentration doesn't change in lungs
% Finding volumetric flow rate out of sodium
Cout(6) = Cvector(6); %sodium out, L/min

% Finding volumetric flow rate out of erythrocytes
Cout(1) = Cvector(1);

% Finding volumetric flow rate out of glucose

nglucosecons = nO2cons/6; %molar flow rate of glucose consumed, mol/min, in terms of O2 consumed
nGlucosei = Cvector(5)*vblood;
nGlucosej = nGlucosei - nglucosecons; %molar flow rate of glucose out, mol/min
Cout(5) = nGlucosej/vblood;
bloodout = vblood;
end

function [bloodflowj, cvectorj] = brain(bloodflowi, cvectori, mass, Ci, basehemoglobin, hemoglobin)
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

%Consume oxygen at 18.4% of intake rate
MO2j=Mvector(2)-0.184*Ci;
%Consume glucose at a rate of 5.6 mg per 100g brain mass per minute
MGlucosej=Mvector(5)-(0.00003108*(brainmass/100));%moles of glucose in - glucose requirements per mon in moles/minute
%Calculate CO2 produced based on glucose consumed via respiration equation
MCO2j=Mvector(3)+6*(0.6660746/1440);
%Calculate HCO3 made when the CO2 is produced
rHCO3CO2=Mvector(4)/Mvector(3); %ratio of bicarbonate to carbon dioxide in blood leaving lungs
HCO3j=Mvector(4)+rHCO3CO2*6*(0.6660746/1440);

%Calculate values above as concentrations and add them to cvectorj
cvectorj(2)=MO2j/bloodflowj;
cvectorj(5)=MGlucosej/bloodflowj;
cvectorj(3)=MCO2j/bloodflowj;
cvectorj(4)=HCO3j/bloodflowj;

end

function [outflow, Cout] = heart(flow,Cvector,weight, Ci, basehemoglobin, hemoglobin)

Cout = [];

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

heartmass = weight*.003;
nO2cons = 0.116*Ci/(hemoglobin/basehemoglobin); %mL/min*mol/mL=mol/min
Cout(2) = (Cvector(2)*flow - nO2cons)/flow;
CO2gen = .7*nO2cons; %using respiratory quotient in the heart
Cout(3) = (Cvector(3)*flow + CO2gen)/flow;
rHCO3CO2 = Cvector(4)/Cvector(3); %ratio of bicarbonate to carbon dioxide in blood leaving lungs
nHCO3j = rHCO3CO2*CO2gen; %HCO3 out, mol/min
Cout(4) = (Cvector(4)*flow + nHCO3j)/flow;

outflow = flow;
end




function [bloodflowj, cvectorout, FerritinStores1]=liver(V,cvectori,G,mass,Ci,FerritinStores, basehemoglobin, hemoglobin)
%This function will deliver output volumetric flow rates for the concentration of
%the 8 components (mol/mLmin) out of the liver
%Input: G=Gender (1 if female, 0 if male)
%W=weight of person (kg); A=age of person (years)
%cvector=[cE,cO2,cCO2,cHCO3,cGlucose,cNa,cCa,cFe)
%Outputs are all in concentration flow rates (mol/mLmin) out

%Initializes cvectorout
cvectorout=[0,0,0,0,0,0,0,0];



%Calculates concentration flow rate of red blood cells out
%p=blood density in amount per mL
%Vbody=Volume of body in mL. Calculating assuming 7% of total body weight
%is blood, then divided by density of blood to find volume of blood
%RBCj=raw numerical flow rate of red blood cells out/min
%RBCi=raw numerical flow rate of RBC in/min
%nRBCj=molar flow rate of RBC out/min
%RBCcons=consumption rate of RBC/min
%No generation of RBC in liver

p=1.2e10; %number of red blood cells per mL of PURE red blood cells

%Calculate volume of the body
Vblood=0.07*mass*1000/1.056;
%Calculate the number of RBCs entering liver based on concentration,
%flow rate, and pure RBC density
RBCi=p*V*cvectori(1);

%Consume some based on body size
RBCcons=5.55555e-6*cvectori(1)*Vblood;
%Find pure volume of RBCs we now have based on number and density
VRBCj=(RBCi-RBCcons)/p;
%Calculate volume of plasma that entered
Vplasma=V*(1-cvectori(1));
%Calculate new concentrations based on volumeRBCs/total volume
cvectorout(1)=VRBCj/(VRBCj+Vplasma);
%Bloodflow out = blood cells out + plasma volume
bloodflowj=V;%Vplasma+VRBCj;

%Calculates concentration flow rate of O2 out
%No generation of O2
%nO2i=molar flow rate of O2/min in
%nO2j=molar flow rate of O2/min out
%nO2cons=molar consumption rate of O2/min
nO2i= cvectori(2)*V;
nO2cons=0.204*Ci/(hemoglobin/basehemoglobin);
nO2j=nO2i-nO2cons;
cvectorout(2)=nO2j/V;

%Calculates the required intake of protein in grams/min given a
%person's demographic information
%dprotein=daily protein intake in g/day
%protein=protein intake in g/min
dprotein=0.8*mass;
protein=dprotein/1440;

%Calculates concentration flow rate of CO2 out
%nCO2i=molar flow rate of CO2/min in
%nCO2j=molar flow rate of CO2/min out
%nCO2gen=molar flow rate of CO2/min generated
%RQ=Respiratory Quotient=mol CO2 produced/mol O2 consumed
RQ=0.825;
rHCO3CO2 = cvectori(4)/cvectori(3);
nCO2gen=RQ*nO2cons;
nCO2i=cvectori(3)*V;
nCO2j=nCO2i+nCO2gen-((protein/100)/rHCO3CO2);
cvectorout(3)=nCO2j/V;


%Calculates the concentration flow rate of Bicarbonate out
%Assume no generation of bicarbonate
%%Also liver just takes in bicarbonate from the heart and stomach
%1mole bicarbonate ions consumed by liver per 100g protein/day
%rHCO3CO2 = cvectori(4)/cvectori(3); %ratio of bicarbonate to carbon dioxide in blood leaving lungs
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
%cCai=concentration flow rate of Ca2+ in mol/mLmin in
%cCaj=concentration flow rate of Ca2+ in mol/mLmin out
cvectorout(7)=cvectori(7);
%Calculates the concentration flow rate of Fe out
%Assume no generation since

% Glucose needs a big ol generation term
consumeorgenerate=0.8144*(cvectori(5)*10^6)^3-9.622*(cvectori(5)*10^6)^2+108.74*(cvectori(5)*10^6)-480.69;
scalingterm = 300/(26*mass*10^-6);
if consumeorgenerate < 0
    nGlucosegen=abs(consumeorgenerate)/scalingterm;
    nGlucosecons=0;
elseif consumeorgenerate > 0
    nGlucosecons=abs(consumeorgenerate)/scalingterm;
    nGlucosegen=0;
else
    nGlucosegen=0;
    nGlucosecons=0;
end
nGlucosebasalcons = nO2cons/6;
nGlucosei = cvectori(5)*V;
nGlucosej = nGlucosei + nGlucosegen - nGlucosecons - nGlucosebasalcons;
cvectorout(5) = nGlucosej/V;

% Iron gets packaged into ferritin in the liver, to be stored and used
% later if we run low


%determine molar flow rate of iron in
Mironin=V*cvectori(8);
Mirondif=(V*7.2231e-06)-Mironin;%this number is desired steady state iron concentration in the blood

if Mirondif < 0
    FerritinStores1=FerritinStores-Mirondif;
    cvectorout(8)=7.2231e-06;
elseif Mirondif == 0
    cvectorout(8)=cvectori(8);
    FerritinStores1=FerritinStores;
    
elseif FerritinStores > Mirondif && Mirondif>0
    Mironadd=Mirondif;
    FerritinStores1=FerritinStores-Mironadd;
    cvectorout(8)=(Mironadd+Mironin)/V;

elseif FerritinStores < Mirondif && Mirondif>0
    Mironadd=FerritinStores;
    cvectorout(8)=(Mironadd+Mironin)/V;
    FerritinStores1=0;
    
elseif FerritinStores == Mirondif && Mirondif>0
    FerritinStores1=0;
    Mironadd=FerritinStores;
    cvectorout(8)=(Mironadd+Mironin)/V;

elseif FerritinStores == 0 && Mirondif>0
    FerritinStores1=0;
    cvectorout(8)=cvectori(8);
end

end


function [bloodflowj, cvectorj] = kidney(bloodflowi, cvectori, RQ, mass, Ci, V, basehemoglobin, hemoglobin, anemia)
T=1440;%Multiplier that scales up the time period of interest to one day (required for glucose equation)
Kidneymass=(300/70000)*mass; %mass of the kidneys combined in grams

%Initialize all intakes
Mvector=bloodflowi*cvectori;
bloodflowj=bloodflowi;

%Consume relevant compounds to create outlet values

MEj=Mvector(1);
MIronj=Mvector(8);

MO2j=Mvector(2)-(0.072*Ci/(hemoglobin/basehemoglobin));
% Find mol/min O2 consumed
%MO2cons = Mvector(2) - MO2j;
% Convert to mL/min per 100g
%MO2cons = (MO2cons/0.0000392927)/(Kidneymass/100);
% Find mmol/min per 100g Na consumed
%MNacons = (MO2cons - 0.5)/0.1;
% Convert to mol/min
%MNacons = (MNacons/1000)*(Kidneymass/100);

%this is new sodium stuff

PNaconc=cvectori(6)*1e6; %plasma sodium concentration, converted to mmol/L for the sake of the equation
%Naremoved=0.0043*(PNaconc)^3 - 1.6942*(PNaconc)^2 + 223.58*PNaconc - 9925; %sodium removed, micromol per minute
Naremoved=0.0003*exp(0.0817*PNaconc);%sodium removed, micromol per minute
MNaremoved=Naremoved/1e6;%sodium removed in moles
MNaj = Mvector(6) - MNaremoved;

%this is old sodium stuff
% MNacons = (100/1440)/1000;
% MNaj = Mvector(6) - MNacons;




x = (1000*.26)/(1000*40.08*1440);
y = (17/5)*(bloodflowi-((3/34)*V));
b = (((x/y)*(y-(.7*V))+(((.7*x)/y)*V))*(y+(.3*V)-bloodflowi))/((y+(.3*V)));
c = (1-(bloodflowi/(y+(.3*V)-bloodflowi)));
MCaj= Mvector(7) - ((((x/y)*(y-(.7*V))+(((.7*x)/y)*V))*(y+(.3*V)))/((y+(.3*V))));
MGlucosej=Mvector(5)-(Mvector(5)*(0.226/((Mvector(5))*T)));
rHCO3O2=Mvector(4)/Mvector(3);
if anemia == 1
    MHCO3cons = 2.9e-06*bloodflowi;
    MCO2cons = (MHCO3cons/rHCO3O2);
else
    MHCO3cons = 0;
    MCO2cons = 0;
end
MCO2j=Mvector(3)+((Mvector(2)-MO2j)*RQ)-MCO2cons;
MHCO3gen=rHCO3O2*((Mvector(2)-MO2j)*RQ);
MHCO3j=Mvector(4)+MHCO3gen-MHCO3cons;%-0.15*(Mvector(4))+(.004/T)
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

function [outflow, Cout, hemoglobinout] = otherblood(flow, Cin, carbs, calciumintake, sodiumintake, ironintake, Ci, RQ, anemia, bleed, basehemoglobin, hemoglobin)
Cout = [];
outflow=flow;

%Accounts for the bleeding of a certain volume at a steady rate when the
%patient is anemic    
if anemia == 0
    Cout(1) = Cin(1);
    ironout=0;
    hemoglobinout=0;
else
    RBCin = flow*Cin(1)*1.2e10;
    vRBCin = flow*Cin(1);
    RBClost = 521400000;%5214000; %bleed*Cin(1)*1.2e10;
    RBCout = RBCin - RBClost;
    vRBCout = RBCout/(1.2e10);
    ironout=1.197929*RBClost/1.2e10*1.79067e-5;
    hemoglobinout=ironout*1.538461538e-5;
Cout(1) = vRBCout/outflow; 
end

Cout(2) = (Cin(2)*flow-(0.374*Ci/(hemoglobin/basehemoglobin)))/flow;%-CO2cons;
nCO2gen = 0.374*Ci*RQ;
Cout(3) = (Cin(3)*flow + nCO2gen)/flow;
rHCO3CO2 = Cin(4)/Cin(3);
Cout(4) = (Cin(4)*flow + nCO2gen*rHCO3CO2)/flow;
Cout(5) = (flow*Cin(5) + (.8*carbs)/(180.156*1440))/flow; %130g of carbs, 80 percent is directly translated to glucose, 20% is fructose and galactose which goes to the liver and may be converted to glucose at a later stage so another term is incoming
Cout(6) = (flow*Cin(6) + (sodiumintake/(700*22.99*1440)))/flow; %500 mg of sodium needed for vital functions so start with sodium intake being 500 mg
if anemia == 0
    Cout(7) = (flow*Cin(7) + (calciumintake*.26)/(1000*40.08*1440))/flow; %calcium intake should be around 1000 mg
elseif anemia == 1
    Cainput = (flow*Cin(7) + (calciumintake*.15)/(1000*40.08*1440))/flow; %assume, under anemia, calcium retention rate becomes 15%
    Cagen = ((calciumintake*.11)/(1000*40.08*1440))/flow; %intake from bones because, when we retain less calcuim from diet, the blood maintains steady concentration of calcium by taking calcium from bones
    Cout(7) = Cainput + Cagen;
end
Cout(8) = (flow*Cin(8)+ (ironintake*.18)/(1000*55.845*1440)-1/(1440*1000*55.845)-(ironout))/flow;%iron intake will be approximately 8 mg for males and 18 mg for females
end