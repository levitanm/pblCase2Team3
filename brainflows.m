 function [vO2j, vCO2j, vHCO3j, vCaj, vFej, vNaj, vErythrocytesj, vGlucosej] = lungs(W, vCai, vNai, vGlucosei)
    % This function finds the output amount of each component in the blood from the brain
    % Input: W = weight of human; vCai = volumetric flow rate of calcium going in;
    % vNai = volumetric flow rate of sodium going in; vGlucosei =
    % volumetric flow rate of glucose going in
    % Output: vO2j = volumetric flow rate of oxygen going out; vCO2 =
    % volumetric flow rate of carbon dioxide going out; vHCO3j = volumetric
    % flow rate of bicarbonate going out; vCaj = volumetric flow rate of
    % calcium going out; vFej = volumetric flow rate of iron going out; 
    % vNaj = volumetric flow rate of sodium going out; vErythrocytesj = 
    % volumetric flow rate of erythrocytes going out; vGlucosej = 
    % volumetric flow rate of glucose going out
    
    % Finding volumetric flow rate of blood (volumetric flow rate of blood
    % in = out)
    mblood = 0.07*W; %in kg
    pblood = 1.06; %perhaps we state this as one of our assumptions. If this is the average density of blood then it
    %seems plausible that we could use this for steady-state even if we do
    %have to change it for the anemic model
    vblood = mblood/pblood; %in L/min
    
    % Finding volumetric flow rate out of oxygen
    mBrain=.02*W %assuming that brain is 2% of human body weight
    CiO2 = 0.05; %5mL/100mL %this should be dependent on hemoglobin - find more on this later
    vO2i = vblood*(CiO2)
    vO2cons = mBrain*.035; %Brain consumes 3.5 ml of oxygen per gram of brain 
    %tissue per minute
    vO2j = vO2i - vO2cons; %O2 out, L/min
    
    % Finding volumetric flow rate out of carbon dioxide
    CO2con=.0814 %this is in grams per minute and takes the 80 grams per day consumption
    %of glucose and then calculates using the cellular respiration equation
    vCO2i=(CiCO2*vblood)
    vCO2j-vCO2i-CO2con %this is not going to work unit wise but the thought is there
    
    % Finding volumetric flow rate out of bicarbonate
    rHCO3CO2 = 19.3/21.5; %ratio of bicarbonate to carbon dioxide in blood leaving lungs
    vHCO3j = rHCO3CO2*vCO2j; %HCO3 out, L/min
    
    % Finding volumetric flow rate out of calcium
    vCaj = vCai; %calcium out, L/min

    % Finding volumetric flow rate out of sodium
    vNaj = vNai; %sodium out, L/min
    
    % Finding volumetric flow rate out of erythrocytes
    vErythrocytesj = vErythrocytesi %I am pretty sure that the concentration
    %of erythrocytes does not change as blood passes through the brain and
    %thus inflow should equal outflow
    
    % Finding volumetric flow rate out of glucose
    Gcon=80/(24x60) %this shows the amount used in grams so I am not sure how 
    %to get this component in terms of volume
    Glucosei=1*vBlood %this is assuming standard glucose level of 100 mg/dL
    Glucosej=Glucosei-Gcon