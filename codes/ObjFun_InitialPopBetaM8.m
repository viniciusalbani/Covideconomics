function [f] = ObjFun_InitialPopBetaM8(t_actual,params,data,options,priors,...
                            yinit,Proportion,PropInfections,unknowns,pesos)
NP = params.NumberOfPlaces;
I_M = unknowns(1:NP);
Number = params.NumberOfAgeClasses;
a = unknowns(NP+1:2*NP);
params.a = a;
b = params.b;
c = params.c;

dist = params.dist;
aux0 = unknowns(2*NP+1:end);
aux = zeros(NP);
aux(dist==0)=ones;
aux((dist>=1 & dist<1000))=aux0(1)*ones;
aux((dist>=1000 & dist<2000))=aux0(2)*ones;
aux((dist>=2000 & dist<=3000))=aux0(3)*ones;
aux((dist>=3000 & dist<=4000))=aux0(4)*ones;
aux(dist>=4000)= aux0(5)*ones;
distance = aux;

aux = zeros(NP);
AUX = zeros(NP*Number);
for ii=1:NP
aux(ii,ii:end) = 0.5*(a(ii)*ones+a(ii:end)).*distance(ii,ii:end);
aux(ii,ii) = 0.5*aux(ii,ii);
aux2 = reshape(ones(3,1)*aux(ii,ii:end),1,3*length(aux(ii,ii:end)));
for jj = 1:Number
AUX((ii-1)*Number+jj,(ii-1)*Number+jj:end) = aux2(jj:end);
end
end

beta_M = AUX + AUX';

params.beta_M = beta_M;      % In Mild conditions
params.beta_H = b*beta_M;  % Hospitalized individuals
params.beta_I = c*beta_M; % In ICU

%%%%%

tspan = [t_actual(1),t_actual(end)];

N = params.N;
for jj = 1:NP
yinit(6*NP*Number+(jj-1)*Number+1:6*NP*Number+jj*Number) = ...
                       I_M(jj)*PropInfections((jj-1)*Number+1:jj*Number);
yinit((jj-1)*Number+1:jj*Number) = Proportion((jj-1)*Number+1:jj*Number)-...
                       I_M(jj)*PropInfections((jj-1)*Number+1:jj*Number);
end
sigma = params.sigma;

[t,y]=ode45(@(t,y)seir_death_age_beta_b4(t,y, params),tspan,yinit,options);

NewInfections = zeros(length(t_actual)-1,NP);
f =[];
for jj = 1:NP
aux = zeros;
for ii = 1:Number
aux = aux+sigma*y(:,4*(NP*Number)+ii+(jj-1)*Number)*N + sigma*y(:,5*(NP*Number)+ii+(jj-1)*Number)*N;
end
NewInfections(:,jj) = interp1(t,aux,t_actual(2:end)');

aux = data(:,jj);
% Stirling = 0.5*log(2*pi*aux) + aux.*log(aux) - aux;
% f = [f;aux.*log(NewInfections(:,jj)) - NewInfections(:,jj) - Stirling];
% f = [f;(aux-NewInfections(:,jj))/pesos(jj)];
f = [f;(log(abs(1+aux))-log(abs(1+NewInfections(:,jj))))];
end
f = [f;1E-3*(unknowns-priors)'];
% f = [f;zeros*(unknowns)'];

f(isnan(f))=zeros;