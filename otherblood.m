function [t, vO2j, vCO2j, vCaj, vFej, vHCO3j, vGlucosej, vNaj, vErythrocytesj] = otherblood(t, vO2i, vCO2i, vCai, vGlucosei, vFei, vHCO3i, vNai, vErythrocytesi)
%for all of the components, adjustments should be made based on what is
%being inputted to other blood and what the output values should be to the
%other organs. For example, when t = 2, blood is entering from the heart
%and then this blood is entering the liver so, based on the output value of
%the components from the heart and the input values needed for the liver,
%adjustments will be made to each component in other blood. 
%Another direction we could take is trying to find specific values for adjustments
%of each component at each time (heart to liver, liver to kidneys, kidneys
%to heart) but this will take a lot of research to determine. 
if t == 2
    %This will be the blood traveling to the liver
    %vblood = ?, do I just assume it is the vblood from the heart or should
    %I adjust it based on the volume of blood in 'otherblood'?
    % CGlucose = x; %this will be the concentration of glucose entering
    % from the digestive tract, will also depend on how much nutrients the
    % patient is eating
    % vGlucosegen = CGlucose * vblood
    % CCa = y %concentration of calcium entering from bones
    % vCagen = CCa * vblood
    vErythrocytesj = vErythrocytesi * 2000000 * 60; %number of RBC entering through the blood marrow per min (though in research it does mention that the same number is being recycled/cleared so this probably needs to be changed)
    vO2j = vO2i;
    vCO2j = vCO2i;
    vCaj = vCai * vCagen;
    vNaj = vNai;
    vGlucosej = vGlucosei + VGlucosegen; %VGlucosegen will be addition of glucose from digestive tract, still need to find the value
    vFej = vFei;
    vHCO3j = vHCO3i;
    t = t + 1;
elseif t == 4
    %blood is entering from liver and leaving to kidneys
    vErythrocytesj = vErythrocytesi;
    vO2j = vO2i;
    vCO2j = vCO2i;
    vCaj = vCai;
    vNaj = vNai;
    vGlucosej = vGlucosei;
    vFej = vFei;
    vHCO3j = vHCO3i;
    t = t + 1;
elseif t == 6
    %blood is entering from kidneys and going to heart
    vErythrocytesj = vErythrocytesi;
    vO2j = vO2i;
    vCO2j = vCO2i;
    vCaj = vCai;
    vNaj = vNai;
    vGlucosej = vGlucosei;
    vFej = vFei;
    vHCO3j = vHCO3i;
    t = t + 1;
end
end