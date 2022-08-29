function f = objFunLabourMkt20220623(t_actual,unknowns,data,options,yinit,CasesByMonth,prior)

params.rhoE = unknowns(1);
params.rhoU = unknowns(2);
params.lambdaE = unknowns(3);
params.lambdaU = unknowns(4);

tspan = [t_actual(1),t_actual(end)];
Infections = @(t)interp1(t_actual,CasesByMonth,t);
[t,y] = ode45(@(t,y)LabourMkt20220623(t,y, params,Infections),tspan,yinit,options);
Em = interp1(t,y(:,1),t_actual);
Un = interp1(t,y(:,2),t_actual);
% aux = (norm(data(1,:))+norm(data(2,:)));
% f = [(data(1,:)-Em)/aux,(data(2,:)-Un)/aux,1E-5*(unknowns)/norm(unknowns)];
f = [data(1,:)-Em,data(2,:)-Un,1E-2*(unknowns-prior)];