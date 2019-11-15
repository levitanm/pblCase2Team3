%im dumb and i dont know how to code
%this broken heart is gonna kill someone

function [vO2j, vCO2j, vHCO3j, vCaj, vFej, vNaj, vErythrocytesj, vGlucosej] = heart(flow, vO2i, vCO2i, vHCO3i, vCai, vFei, vNai, vErythrocytesi, vGlucosei)
    % This function finds the output terms for each component coming out of
    % the heart.
    % Input: 
    % Output: vO2j = volumetric flow rate of oxygen going out; vCO2 =
    % volumetric flow rate of carbon dioxide going out; vHCO3j = volumetric
    % flow rate of bicarbonate going out; vCaj = volumetric flow rate of
    % calcium going out; vFej = volumetric flow rate of iron going out; 
    % vNaj = volumetric flow rate of sodium going out; vErythrocytesj = 
    % volumetric flow rate of erythrocytes going out; vGlucosej = 
    % volumetric flow rate of glucose going out
    
    %how do I account for total blood/blood flow
    %so the Heart has multiple flows to it, but the heart will actually consume only once -- the first time because we are assuming that everything is happening in one heartbeat? 
    if flow == Flow1LungsToHeart
    
        % components that don't change:
        
        % Finding volumetric flow rate out of glucose
        vGlucosej=vGlucosei;
        % Finding volumetric flow rate out of erythrocytes
        vErythrocytesj=vErythrocytesi;
        % Finding volumetric flow rate out of sodium
        vNaj=vNai;
        % Finding volumetric flow rate out of iron
        vFej=vFei;
        % Finding volumetric flow rate out of calcium
        vCaj=vCai;
        % Finding volumetric flow rate out of bicarbonate
        vHCO3j=vHCO3i;

        %components that do change
        % Finding volumetric flow rate out of oxygen
        % Finding volumetric flow rate out of carbon dioxide
        heartRQ=0.7; %we are assuming that heart muscle cell metabolism comes soley from lipids
        heartRQ=CO2produced/O2consumed;
        % O2 consumption: 30 - 35 mL/min per 300 g
        %based off of RQ find the CO2produced
        %vO2j=vO2-O2consumed but do units match
        %vCO2j=vCO2 + CO2produced
       
    else 
        % no components will change?
        %Flow2HeartToOther
        %Flow2HeartToBrain
        %Flow2HeartToLiver
        %Flow3BrainToHeart
        %Flow4LiverToHeart
        %Flow7OtherToHeart
        %Flow8HeartToLungs
        
        vGlucosej=vGlucosei;
        vErythrocytesj=vErythrocytesi;
        vNaj=vNai;
        vFej=vFei;
        vCaj=vCai;
        vHCO3j=vHCO3i;
        vCO2j=vCO2i;
        vO2j=vO2i;
    end

end