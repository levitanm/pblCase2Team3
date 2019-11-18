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
    nO2cons = .0002082515; %5.3 mL/min becomes .0002082515 mol/min
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
     Cerythrocytes = 0.45; %45 mL/100 mL, this concentration changes depending on hemoglobin, but we need...
%                           %to figure out this relationship
     Chemoglobin = 0.335; %g/mL, this concentration also changes, but maybe only depending on demographics
%                          %and anemia
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
    Chemoglobinblood = Chemoglobin*Cerythrocytes; %concentration of hemoglobin in blood, g/mL
    Chemoglobinbloodmol = Chemoglobinblood/Mhemoglobin; %concentration of hemoglobin in blood, mol/mL
    CO2bloodmol = 4*Chemoglobinbloodmol; %concentration of oxygen in blood dependent on this value of...
                                         %hemoglobin, mol/mL
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