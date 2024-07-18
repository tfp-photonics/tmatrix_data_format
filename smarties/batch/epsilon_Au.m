function e = epsilon_Au(lambda)
% 10.1103/PhysRevB.86.235147

tmp = load('epsAu2012SCdata.mat');
epsload = tmp.epsAu2012SCdata(1:164,:); % 1:164 to avoid repeated entries

e1 = interp1(epsload(:,1),epsload(:,2),lambda,'linear');
e2 = interp1(epsload(:,1),epsload(:,3),lambda,'linear');
e=e1+1i*e2;

end

