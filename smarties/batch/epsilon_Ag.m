function e = epsilon_Ag(lambda)
% 10.1103/PhysRevB.91.235137

tmp = load('epsAg2015PRBCCdata.mat');
epsload = tmp.epsAg2015PRBCCdata;

e1 = interp1(epsload(:,1),epsload(:,2),lambda,"linear");
e2 = interp1(epsload(:,1),epsload(:,3),lambda,"linear");
e=e1+1i*e2;

end

