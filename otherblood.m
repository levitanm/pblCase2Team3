function [Cout] = otherblood(flow, Cin, carbs, calciumintake, sodiumintake, ironintake)
Cout = [];
Cout(1) = Cin(1);
Cout(2) = Cin(2);
Cout(3) = Cin(3);
Cout(4) = Cin(4);
Cout(5) = (flow*Cin(5) + (.8*carbs)/180.156)/flow; %80 percent is directly translated to glucose, 20% is fructose and galactose which goes to the liver and may be converted to glucose at a later stage so another term is incoming
Cout(6) = (flow*Cin(6) + (sodiumintake/(1000*22.99)))/flow; %500 mg of sodium needed for vital functions so start with sodium intake being 500 mg
Cout(7) = (flow*Cin(7) + (calciumintake*.26)/(1000*40.08))/flow; %calcium intake should be around 1000 mg
Cout(8) = (flow*Cin(8) + (ironintake*.18)/(1000*55.845))/flow;%iron intake will be approximately 8 mg for males and 18 mg for females
end