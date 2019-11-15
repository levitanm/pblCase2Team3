function [vO2j, vCO2j, vHCO3j, vCaj, vFej, vNaj, vErythrocytesj, vGlucosej] = lungs(W, vCai, vNai, vGlucosei)
    % This function finds the output terms for each component coming out of
    % the lungs.
    % Input: W = weight; vCai = volumetric flow rate of calcium going in;
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
    pblood = 1.06; %in kg/L %this value is dependent on other factors - can we use it? should we use 1?
    vblood = mblood/pblood; %in L/min
    
    % Finding volumetric flow rate out of oxygen
    CiO2 = 0.05; %5mL/100mL %this should be dependent on hemoglobin - find more on this later
    CdeoxygenatedO2 = 0.16; %from graph and partial pressure of oxygen in entering deoxygenated blood being...
                            %40 mmHg (this might also depend on hemoglobin)
    vO2i = vblood*(CiO2+CdeoxygenatedO2); %O2 in, L/min
    vO2cons = 0.0053; %O2 consumed (should depend on hemoglobin, other things), L/min
    vO2j = vO2i - vO2cons; %O2 out, L/min
    
    % Finding volumetric flow rate out of carbon dioxide
    CjCO2 = 0.48; %48 mL/100 mL, should be dependent on hemoglobin, oxygen, things like that
    vCO2j = vblood*CjCO2; %CO2 out, L/min
    
    % Finding volumetric flow rate out of bicarbonate
    rHCO3CO2 = 19.3/21.5; %ratio of bicarbonate to carbon dioxide in blood leaving lungs
    vHCO3j = rHCO3CO2*vCO2j; %HCO3 out, L/min
    
    % Finding volumetric flow rate out of calcium
    vCaj = vCai; %calcium out, L/min
    
    % Finding volumetric flow rate out of iron
    Cerythrocytes = 0.45; %45 mL/100 mL, this concentration changes depending on hemoglobin, but we need...
                          %to figure out this relationship
    Chemoglobin = 0.335; %g/mL, this concentration also changes, but maybe only depending on demographics
                         %and anemia
    Mhemoglobin = 65000; %g/mol
    MFe = 55.845; %g/mol
    CFe = (MFe*Chemoglobin*Cerythrocytes)/Mhemoglobin;
    vFej = vblood*CFe; %iron out, L/min
    
    % Finding volumetric flow rate out of sodium
    vNaj = vNai; %sodium out, L/min
    
    % Finding volumetric flow rate out of erythrocytes
    vErythrocytesj = vblood*Cerythrocytes; %erythrocytes out, L/min, once again this concentration changes
    
    % Finding volumetric flow rate out of glucose
    Chemoglobinblood = Chemoglobin*Cerythrocytes; %concentration of hemoglobin in blood, g/mL
    Chemoglobinbloodmol = Chemoglobinblood/Mhemoglobin; %concentration of hemoglobin in blood, mol/mL
    CO2bloodmol = 4*Chemoglobinbloodmol; %concentration of oxygen in blood dependent on this value of...
                                         %hemoglobin, mol/mL
    nO2cons = CO2bloodmol*vblood*1000; %molar flow rate of oxygen consumed, mol/min - can we somehow...
                                       %convert this into mL/min to get a
                                       %more accurate value of oxygen
                                       %consumed in the oxygen accounting
                                       %equation?
    nglucosecons = 6*nO2cons; %molar flow rate of glucose consumed, mol/min
    Mglucose = 180.18; %g/mol
    mglucosecons = nglucosecons*Mglucose; %g/min
    pglucose = 1560; %g/L, is this what other people are using?
    vglucosecons = mglucosecons/pglucose; %L/min
    vGlucosej = vGlucosei - vglucosecons; %volumetric flow rate of glucose out, L/min
    
    %*thoughts I had while coding that we should consider: putting in
    %checks to make sure outflows that go directly to the next organ are
    %the same as those inflows - if they're not the same, find the problem
    %and/or change it to be the same
    % make sure everyone codes volumetric flow rate of blood and density of
    % blood in the same way