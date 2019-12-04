function oxygeninlungs()
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
%convert ViO2 to mol/min
niO2 = (ViO2*1)/(0.08206*310.15); %mol/min; convert mL to mol with pressure being atmospheric pressure
niO2 = niO2/15; %divide molar flow rate by respiratory rate
C = niO2/(vblood*1.5); %mol/mL, divide by conversion factor of 1.5

Ci = C*vblood; %concentration of oxygen from which other organs consume

nO2i = vblood*(Cvector(2) + C);
nO2cons = 0.05*C*vblood;
nO2j = nO2i - nO2cons; %oxygen out, mol/min
Cout(2) = nO2j/vblood; %mol/mL