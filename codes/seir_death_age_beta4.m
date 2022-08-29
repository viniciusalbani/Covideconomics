function  yprime = seir_death_age_beta4(t,y,params,beta)

factor = params.factorWorse;
factorI = params.factorICU;
factorD = params.factorDeath;
fh = factor(t);
fI = factorI(t);
fd = factorD(t);


gamma_M = params.GetWorse_M;
gamma_H = params.GetWorse_H;
mu_I = params.Death_I;

NP = params.NumberOfPlaces;
be = beta(t);
Number = params.NumberOfAgeClasses;
b = params.b;
c = params.c;

beta_M = params.beta_M;
for jj = 1:NP
beta_M((jj-1)*Number+1:jj*Number,(jj-1)*Number+1:jj*Number) = ...
        be(jj)*beta_M((jj-1)*Number+1:jj*Number,(jj-1)*Number+1:jj*Number);
%%%%
end
beta_H = b*beta_M;
beta_I = c*beta_M;

nu_M = zeros(1,NP*Number);
% N = params.N/3;
sigma = params.sigma;
nu_H = params.nu_H*ones(size(nu_M));%max(0,ones-gamma_H);%params.Recovery_M;%params.Recovery_H;
% mu_M = max(0,fd-mu_I.*(gamma_M.*gamma_H));%params.Death_M;
% nu_M = max(0,ones-gamma_M-mu_M);%params.Recovery_M;
mu_H = params.Death_H;
nu_I = params.nu_I*ones(size(nu_M));%ones-mu_I;
%%% New Parameters (time-dependent):
LAMBDA_UN = params.lambda_Un;
LAMBDA_EM = params.lambda_Em;

lambda_Un1 = LAMBDA_UN(t);
lambda_Em1 = LAMBDA_EM(t);
lambda_Un = zeros(1,NP*Number);
lambda_Em = zeros(1,NP*Number);
NU = params.nu;
nu1 = NU(t);
nu = zeros(1,NP*Number);
% AUX = zeros(1,NP*Number);
for jj = 1:NP
lambda_Un((jj-1)*Number+1:jj*Number) = lambda_Un1(jj)*[0,1,0];
lambda_Em((jj-1)*Number+1:jj*Number) = lambda_Em1(jj)*[0,1,0];
nu((jj-1)*Number+1:jj*Number) = nu1(jj)*[0.063931519819463,0.767207875043092,0.168860605137446];
% AUX((jj-1)*Number+1:jj*Number) = sum(y((jj-1)*Number+1:jj*Number))+sum(y(Number*NP+(jj-1)*Number+1:Number*NP+jj*Number))*ones;
end
eps_Un = params.eps_Un;


%    beta_M = R_zero*gamma; 
Number = Number*NP;

Se = y(1:Number);
Su = y(Number+1:2*Number);
Ve = y(2*Number+1:3*Number);
Vu = y(3*Number+1:4*Number);
Ee = y(4*Number+1:5*Number);
Eu = y(5*Number+1:6*Number);
I_Me = y(6*Number+1:7*Number); 
I_Mu = y(7*Number+1:8*Number);   
I_He = y(8*Number+1:9*Number); 
I_Hu = y(9*Number+1:10*Number); 
I_Ie = y(10*Number+1:11*Number); 
I_Iu = y(11*Number+1:12*Number);
Re = y(12*Number+1:13*Number);
Ru = y(13*Number+1:14*Number);
%      D = y(15);

Number = Number/NP;
for jj = 1:NP
auxM = sum(I_Me((jj-1)*Number+1:jj*Number)+I_Mu((jj-1)*Number+1:jj*Number));
auxH = sum(I_He((jj-1)*Number+1:jj*Number)+I_Hu((jj-1)*Number+1:jj*Number));
gamma_M((jj-1)*Number+1:jj*Number) = min(0.5,fh(jj)*gamma_M((jj-1)*Number+1:jj*Number)/auxM); % Rate of Hospitalization
gamma_H((jj-1)*Number+1:jj*Number) = min(0.5,fI(jj)*gamma_M((jj-1)*Number+1:jj*Number)/auxH); % Rate of ICU Admission
mu_M((jj-1)*Number+1:jj*Number) = min(0.5,max(0,(fd(jj)/auxM)-mu_I((jj-1)*Number+1:jj*Number).*(gamma_M((jj-1)*Number+1:jj*Number).*gamma_H((jj-1)*Number+1:jj*Number))));%params.Death_M;
if auxM==zeros
    gamma_M((jj-1)*Number+1:jj*Number)=zeros;
    mu_M((jj-1)*Number+1:jj*Number) = zeros;
end
if auxH==zeros
    gamma_H((jj-1)*Number+1:jj*Number)=zeros;
end
% mu_M((jj-1)*Number+1:jj*Number) = max(0,fd(jj)*mu_I((jj-1)*Number+1:jj*Number)/0.4);

%nu_M((jj-1)*Number+1:jj*Number) = params.nu_M*ones;%max(0,ones-gamma_M((jj-1)*Number+1:jj*Number)-mu_M((jj-1)*Number+1:jj*Number));%params.Recovery_M;
nu_M((jj-1)*Number+1:jj*Number) = params.nu_M*ones;%max(0,ones-gamma_M((jj-1)*Number+1:jj*Number)-mu_M((jj-1)*Number+1:jj*Number));
end
Number = Number*NP;
yprime = zeros(14*Number,1);
for jj = 1:Number
AUX = Se(jj)+Su(jj);
AUX2 = ones;%Se(jj)+Ve(jj)+Ee(jj)+I_Me(jj)+Re(jj);
AUX3 = ones;%Su(jj)+Vu(jj)+Eu(jj)+I_Mu(jj)+Ru(jj);
if AUX2==0
    AUX2=ones;
end
if AUX3==0
    AUX3 = ones;
end
if nu(jj)/AUX >=1
   nu(jj) = zeros;
end
% if t>678
%     AUX3=ones;
%     AUX2=ones;
% end
Reinf = params.reinf;
% Susceptible:
yprime(jj) = -Se(jj)*(beta_M(jj,:)*(I_Me+eps_Un*I_Mu)+beta_H(jj,:)*(I_He+I_Hu)+beta_I(jj,:)*(I_Ie+I_Iu) + lambda_Un(jj)/AUX2 + nu(jj)/AUX) + lambda_Em(jj)/AUX3*Su(jj) + Reinf*(Ve(jj)+Re(jj));%+0.3*Ve(jj);
yprime(Number+jj) = -Su(jj)*(eps_Un*(beta_M(jj,:)*(I_Me+I_Mu)+beta_H(jj,:)*(I_He+I_Hu)+beta_I(jj,:)*(I_Ie+I_Iu)) + lambda_Em(jj)/AUX3 + nu(jj)/AUX) + lambda_Un(jj)/AUX2*Se(jj) + Reinf*(Vu(jj)+Ru(jj));%+0.3*Vu(jj);
% Vaccinated
yprime(2*Number+jj) = nu(jj)/AUX*Se(jj) - lambda_Un(jj)/AUX2*Ve(jj) + lambda_Em(jj)/AUX3*Vu(jj) -Reinf*Ve(jj);% - 0.3*Ve(jj);
yprime(3*Number+jj) = nu(jj)/AUX*Su(jj) + lambda_Un(jj)/AUX2*Ve(jj) - lambda_Em(jj)/AUX3*Vu(jj) -Reinf*Vu(jj);% - 0.3*Vu(jj);
% Exposed
yprime(4*Number+jj) = Se(jj)*(beta_M(jj,:)*(I_Me+eps_Un*I_Mu)+ beta_H(jj,:)*(I_He+I_Hu) + beta_I(jj,:)*(I_Ie+I_Iu))-(sigma+lambda_Un(jj)/AUX2)*Ee(jj) + lambda_Em(jj)/AUX3*Eu(jj);
yprime(5*Number+jj) = eps_Un*Su(jj)*(beta_M(jj,:)*(I_Me+I_Mu)+ beta_H(jj,:)*(I_He+I_Hu) + beta_I(jj,:)*(I_Ie+I_Iu))-(sigma+lambda_Em(jj)/AUX3)*Eu(jj) + lambda_Un(jj)/AUX2*Ee(jj);
% Mildly Infective
yprime(6*Number+jj) = sigma*Ee(jj)-(nu_M(jj)+mu_M(jj)+gamma_M(jj)+lambda_Un(jj)/AUX2)*I_Me(jj) + lambda_Em(jj)/AUX3*I_Mu(jj);
yprime(7*Number+jj) = sigma*Eu(jj)-(nu_M(jj)+mu_M(jj)+gamma_M(jj)+lambda_Em(jj)/AUX3)*I_Mu(jj) + lambda_Un(jj)/AUX2*I_Me(jj);
% Hospitalized
yprime(8*Number+jj) = gamma_M(jj)*I_Me(jj) - (nu_H(jj)+mu_H(jj)+gamma_H(jj))*I_He(jj);
yprime(9*Number+jj) = gamma_M(jj)*I_Mu(jj) - (nu_H(jj)+mu_H(jj)+gamma_H(jj))*I_Hu(jj);
% In ICU
yprime(10*Number+jj) = gamma_H(jj)*I_He(jj)-(nu_I(jj) + mu_I(jj))*I_Ie(jj);
yprime(11*Number+jj) = gamma_H(jj)*I_Hu(jj)-(nu_I(jj) + mu_I(jj))*I_Iu(jj);
% Recovered
yprime(12*Number+jj) = nu_M(jj)*I_Me(jj)+nu_H(jj)*I_He(jj)+nu_I(jj)*I_Ie(jj) - lambda_Un(jj)/AUX2*Re(jj) + lambda_Em(jj)/AUX3*Ru(jj) - Reinf*Re(jj);
yprime(13*Number+jj) = nu_M(jj)*I_Mu(jj)+nu_H(jj)*I_Hu(jj)+nu_I(jj)*I_Iu(jj) + lambda_Un(jj)/AUX2*Re(jj) - lambda_Em(jj)/AUX3*Ru(jj) - Reinf*Ru(jj);
% Deceased
yprime(14*Number+jj) = mu_M(jj)*(I_Me(jj)+I_Mu(jj))+mu_H(jj)*(I_He(jj)+I_Hu(jj))+mu_I(jj)*(I_Ie(jj)+I_Iu(jj));
end
%%%