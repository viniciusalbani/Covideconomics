function f = ObjFun_beta6(t_actual,params,unknowns,data,options,priors,beta,yinit)
Number = params.NumberOfAgeClasses;
NP = params.NumberOfPlaces;
tspan = [t_actual(1),t_actual(end)];
beta = @(t)betaA([beta;unknowns(1:NP)],t_actual,t,NP);

N = params.N;

sigma = params.sigma;

[~,y] = ode45(@(t,y)seir_death_age_beta4(t,y, params,beta),tspan,yinit,options);

f1 =zeros(1,NP);
for jj = 1:NP
aux1 = zeros;
for ii = 1:Number
aux1 = aux1+sigma*y(end,4*(NP*Number)+ii+(jj-1)*Number)*N+sigma*y(end,5*(NP*Number)+ii+(jj-1)*Number)*N;
end
NewInfections = aux1;
aux1 = data(end,jj);
% Stirling = 0.5*log(2*pi*aux) + aux.*log(aux) - aux;
% f(jj) = aux.*log(NewInfections) - NewInfections - Stirling;
% f(jj) = (aux-NewInfections)/pesos(jj);
f1(jj) = (log(abs(1+aux1))-log(abs(1+NewInfections)));
end

f = [f1,1E-5*(unknowns-priors)];
% f=[f1,zeros*(unknowns)];
f(isnan(f))=zeros;