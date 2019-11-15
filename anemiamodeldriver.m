function anemiamodeldriver()

age = x;
gender = y;
weight = z;

digestivesystem(age, gender, weight) = (ca, fe, na, gluc);

bones() = [];

lungs(o2, co2, bloodflow) = [o2, co2, bloodflow];

heart(o2, co2) = [o2, co2];

brain(gluc, o2, bloodflow) = [gluc, o2, bloodflow];

liver(bicarb, co2, o2) = [bicarb, co2, o2];

kidneys(na, ca, gluc, o2, co2) = [na, ca, gluc, o2, co2];

end

function lungs(o2, co2, bloodflow) = [o2, co2, bloodflow]

end

function heart(o2, co2) = [o2, co2]

end

function brain(gluc, o2, bloodflow) = [gluc, o2, bloodflow]

end

function liver(bicarb, co2, o2) = [bicarb, co2, o2]

end

function kidneys(na, ca, gluc, o2, co2) = [na, ca, gluc, o2, co2]

end

function digestivesystem(age, gender, weight) = (ca, fe, na, gluc)

end

function excretions(age, gender, weight) = ()

end

function bones() = []

end
