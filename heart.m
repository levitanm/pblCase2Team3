%I don't know about modeling different conditions based on flow at a
%certain time, would it be better if we used t (for time, i.e. what step we
%are at)
%Also if carbon dioxide concentration changes, then I am pretty sure
%bicarbonate concentration will also change (by how much I'm not sure, I'm
%tired)
function [outflow, Cout] = heart(flow, Cvector,weight)
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
Cout(4)=Cvector(4);

heartmass = weight*.003;
O2cons = 32.5/(heartmass/300); %mL/min
%vblood = (.07*weight)/1.06; %L/min, if the total volumetric blood flow leaving the lungs is the same as volumetric blood flow entering the heart
Cout(2) = Cvector(2)*flow - O2cons/1000;
CO2gen = .7*O2cons; %using respiratory quotient in the heart
Cout(3) = Cvector(3)*flow + CO2gen/1000;

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