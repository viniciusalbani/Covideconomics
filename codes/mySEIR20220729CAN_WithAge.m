clear all; clc; close all; format long e; tic;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% We estimate the parameters of a SEIR-type epidemiological model by
%%%% using a maximum a posteriori estimator. All the estimation procedures
%%%% are carried out by LSQNONLIN, although the likelihood function (or
%%%% data misfit) is log-Poisson. The model parameters are estimated from
%%%% the daily records of infections and deaths.

load data_20220729;
EmploymentRate = zeros(13,size(ParamsCA,1));
UnemploymentRate = zeros(13,size(ParamsCA,1));
for jj=1:13
EmploymentRate(jj,:) = ParamsCA(1:end,3,jj)'+ParamsCA(:,1,jj)'.*(CASES(jj,N:end));
UnemploymentRate(jj,:) = ParamsCA(1:end,4,jj)'+ParamsCA(:,2,jj)'.*(CASES(jj,N:end));
end
ParamsCA = ParamsCA(195:936,:,:);
DATA_CAN2;
data = cases;
data2 = data;
for jj=4:size(data,1)-3
for ii = 1:size(data,2)
data2(jj,ii) = mean(data(jj-3:jj+3,ii));
end
end
data2 = data2(26:end,:);
data = data(26:end,:);
t_actual = t_actual(26:end);
t_span = t_span(26:end);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% Setup for the ODE solver and the least-square function
tol = 1.e-6;  % ode solver tolerance
% set 'Stats','on' to get more info
% options = odeset('AbsTol', tol,'RelTol',tol,'MaxOrder',5,'Stats','on');

% note: set Refine switch to avoid interpolation
options = odeset('AbsTol', tol,'RelTol',tol,'MaxOrder',5,'Stats',...
                                                         'off','Refine',1);
options2 = [];%optimset('MaxFunEvals',15000,'MaxIter',10000,'TolFun',...
                                                      %  1e-30,'TolX',1e-30);

options3 = [];
% 198:939

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% Experimental Data

Vacc1 = cumsum(data2(:,53:end),1);

% DATA_CAN2;
LaborAdjusted;

%%%% Total Population of CANADA:
N = sum(demo);      
aux = sum(cases);

Number = 3;  % total number of age ranges
NP = 13; % total number of places

% %%%% Population proportion on each age range:
params.PropPopAge = PropPopAge;
for jj = 1:NP
Proportion((jj-1)*Number+1:jj*Number) = (PropPopAge(:,jj)'/sum(PropPopAge(:,jj)))*demo(jj)/N;
PropInfections((jj-1)*Number+1:jj*Number) = aux(jj)*(PropPopAge(:,jj)'/sum(PropPopAge(:,jj)))/demo(jj);
% PropInfections((jj-1)*Number+1:jj*Number) = aux(jj)*ones(1,3)*demo(jj)/N;%(PropPopAge(:,jj)'/sum(PropPopAge(:,jj)))/demo(jj);
PropHosp((jj-1)*Number+1:jj*Number) = ones(1,3);%[1.9,52.1,46]/100;
PropICU((jj-1)*Number+1:jj*Number) =  ones(1,3);%[1.2,65.3,33.5]/100;
PropDeath((jj-1)*Number+1:jj*Number) = 0.4*[0.1,15.3,84.6]/100;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Initial parameters
Su_0 = sum(Unemployed0);%0.055*N;
Ve_0 = zeros;
Vu_0 = zeros;
Ee_0 = 1;       % initial exposed
Eu_0 = 0;
I_Me0 = 1;      % Potential initial infected and mild at t=0
I_Mu0 = 0;      % 
I_He0 = 0;      % initial infected and hospitalized at t=0
I_Hu0 = 0;
I_Ie0 = 0;      % initial infected and in ICU at t=0
I_Iu0 = 0;
Re_0 = 0;       % initial recovered 
Ru_0 = 0;
D_0 = 0;       % initial deceased

%  params is a structure used to pass parameters to the
%   ODE solver

Se_0 = N-(Su_0+Ve_0+Vu_0+Ee_0+Eu_0+I_Me0+I_Mu0+I_He0+I_Hu0+I_Ie0+I_Iu0+Re_0+Ru_0+D_0);    % Suceptible pop. excluding initial infected 
params.N =   N;  % N = total population



aux = Number*NP;
yinit = zeros(1,15*aux);
yinit(1:aux) = Se_0*Proportion;
yinit(aux+1:2*aux) = reshape([zeros(1,13);Unemployed0;zeros(1,13)],1,aux);
yinit(2*aux+1:3*aux) = Ve_0*Proportion;
yinit(3*aux+1:4*aux) = Vu_0*Proportion;
yinit(4*aux+1:5*aux) = (Ee_0*Proportion);
yinit(5*aux+1:6*aux) = (Eu_0*Proportion).*reshape([zeros(1,13);ones(1,13);zeros(1,13)],1,aux);
yinit(6*aux+1:7*aux) = I_Me0*Proportion;
yinit(7*aux+1:8*aux) = (I_Mu0*Proportion).*reshape([zeros(1,13);ones(1,13);zeros(1,13)],1,aux);
yinit(8*aux+1:9*aux) = I_He0*Proportion;
yinit(9*aux+1:10*aux) = I_Hu0*Proportion;
yinit(10*aux+1:11*aux) = I_Ie0*Proportion;
yinit(11*aux+1:12*aux) = I_Iu0*Proportion;
yinit(12*aux+1:13*aux) = Re_0*Proportion;
yinit(13*aux+1:14*aux) = Ru_0*Proportion;
yinit(14*aux+1:15*aux) = D_0*Proportion;
yinit = yinit/N;

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Model Parameters

params.sigma = 1./5.1;   % inverse of mean incubation time
params.NumberOfAgeClasses = Number;
params.NumberOfPlaces = NP;
%--------------------------------------------------------------------------
% Mean time until recovery
nu_M = 1./14;
nu_H = 1./12;
nu_I = 1./9;
params.nu_M = nu_M;
params.nu_H = nu_H;
params.nu_I = nu_I;
%--------------------------------------------------------------------------
% Mean time until death
mu_M = zeros;
mu_H = zeros;
mu_I = 1./7;

%--------------------------------------------------------------------------
% Mean time until passing to another infectious Class
gamma_M = 1./1.2;
gamma_H = 1./3.5;

%--------------------------------------------------------------------------
% Proportion of individuals that will recover
p_M  = ones-PropHosp;
params.p_M = p_M; % In Mild conditions

p_H = ones-PropICU;
params.p_H = p_H; % Hospitalized individuals

p_I = ones-PropDeath;
params.p_I = p_I; % In ICU

%--------------------------------------------------------------------------
% Proportion of individuals that will die
q_M = zeros(1,Number*NP);
params.q_M = q_M;  % In Mild conditions

q_H = zeros(1,Number*NP);
params.q_H = q_H; % Hospitalized individuals

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% REINFECTION PARAMETER:
reinf = 1/360; % Eu estava usando 1/360.
params.reinf = reinf;
%--------------------------------------------------------------------------
%%% RATES
%%% Recovery Rate

Recovery_M = p_M;
Recovery_H = p_H;
Recovery_I = p_I;

params.Recovery_M = Recovery_M; % in Mild conditions
params.Recovery_H = ones-PropICU;%Recovery_H; % Hospitalized individuals
params.Recovery_I = Recovery_I; % in ICU individuals

%%% Getting Worse Rate

GetWorse_M = PropHosp;
GetWorse_H = PropICU;

params.GetWorse_M = GetWorse_M; % Mild conditions
params.GetWorse_H = GetWorse_H; % Hospitalized individuals

%%% Death Rate

Death_M = q_M;
Death_H = q_H;
Death_I = PropDeath;

params.Death_M = Death_M; % in Mild conditions
params.Death_H = Death_H; % Hospitalized individuals
params.Death_I = Death_I; % in ICU individuals

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% Time-dependent rates of hospitalization, ICU admission, and death.
Hosp = zeros(length(t_actual),NP);
Death = zeros(length(t_actual),NP);
ICU = zeros(length(t_actual),NP);
NU = zeros(length(t_actual),NP);
NU(2:end,:) = data2(:,53:end);
UnRate = zeros(length(t_actual),NP);
EmRate = zeros(length(t_actual),NP);
for jj = 1:NP
aux = 0.5*data2(2:end,jj);
Hosp(3:end,jj) = min(aux,data2(2:end,13+jj))/sum(demo);
ICU(3:end,jj) = min(aux,data2(2:end,26+jj))/sum(demo);
Death(3:end,jj) = max(0,min(aux,data2(2:end,39+jj)))/sum(demo);
NU(:,jj) = 0.7*NU(:,jj)/params.N;
UnRate(:,jj) = UnemploymentRate(jj,195:936)';
EmRate(:,jj) = EmploymentRate(jj,195:936)';
end
Hosp(isnan(Hosp)==1)=zeros;
ICU(isnan(ICU)==1)=zeros;
Death(isnan(Death)==1)=zeros;
params.factorWorse = @(t)hospA(Hosp,t_actual,t,NP);
params.factorICU = @(t)ICUA(ICU,t_actual,t,NP);
params.factorDeath = @(t)deathA(Death,t_actual,t,NP);
params.nu = @(t)hospA(NU,t_actual,t,NP); % Vaccination rate
params.eps_Un = 0.5; % Stay-at-home infection reduction
params.lambda_Un = @(t)hospA(UnRate,t_actual,t,NP); % Unemployment rate
params.lambda_Em = @(t)hospA(EmRate,t_actual,t,NP); % Employment rate

%--------------------------------------------------------------------------
% A priori value for the transmission parameter
R_zero = 1.4*2.8/0.1782;
gamma = 1/18;

%--------------------------------------------------------------------------
% Transmission parameters
beta = 2.2911*R_zero*gamma;

aux = zeros(NP);
AUX = zeros(NP*Number);
for ii=1:NP
aux(ii,ii:end) = 0.5*(PropInfections((ii-1)*Number+1)*ones+PropInfections((ii-1)*Number+1:Number:end)).*distance(ii,ii:end);
aux(ii,ii) = 0.5*aux(ii,ii);
for jj = 1:Number
    aux2 = reshape(ones(3,1)*aux(ii,ii:end),1,3*length(aux(ii,ii:end)));
AUX((ii-1)*Number+jj,(ii-1)*Number+jj:end) = aux2(jj:end);
end
end
beta_M = AUX + AUX';
params.beta_M = beta_M;      % In Mild conditions
params.beta_H = beta_M;  % Hospitalized individuals
params.beta_I = beta_M; % In ICU

params.a = ones;
params.b = 0.1;
params.c = 0.01;
params.dist = dist;
paramsOld = params;
yinitOld = yinit;
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Estimating basic parameters from day 1 until 20:
day = 20;

LB1 = [zeros(1,NP),zeros(1,NP)];
UB1 = [50*ones(1,NP)/N,100*ones(1,NP)];
unknowns10=[ones(1,NP)/N,beta*ones(1,NP)];
                          
priors1 = unknowns10;
unknowns20 = [0.5,0.25,0.125,0.625,0.03];
priors2 = unknowns20;
LB2 = zeros(size(unknowns20));
UB2 = 0.5*ones(size(unknowns20));

unknowns0 = [unknowns10,unknowns20]; % initial parameters
LB = [LB1,LB2]; % lower bounds
UB = [UB1,UB2]; % upper bounds

priors = [priors1,priors2];

pesos = ones(1,13);
for jj = 1:13
pesos(jj) = max(1E-1,norm(log(1+data2(1:day-1,jj))));
end

%%% Estimating the transmission constant parameters (M,H,I), the initial
%%% infecve population (I_M0) and the transmission matrix:

OF = @(unknowns)ObjFun_InitialPopBetaM8(t_actual(1:day),params,...
data2(1:day-1,:),options,priors,yinit,Proportion,PropInfections,unknowns,pesos);
unknowns = lsqnonlin(OF,unknowns0,LB,UB,options2);

I_Me = unknowns(1:NP);
N = params.N;
for jj = 1:NP
yinit(6*NP*Number+(jj-1)*Number+1:6*NP*Number+jj*Number) = ...
                       (I_Me(jj)*PropInfections((jj-1)*Number+1:jj*Number));
yinit((jj-1)*Number+1:jj*Number) = Proportion((jj-1)*Number+1:jj*Number)-...
                       I_Me(jj)*PropInfections((jj-1)*Number+1:jj*Number);
end
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
aux2 = reshape(ones(Number,1)*aux(ii,ii:end),1,Number*length(aux(ii,ii:end)));
for jj = 1:Number
AUX((ii-1)*Number+jj,(ii-1)*Number+jj:end) = aux2(jj:end);
end
end

beta_M = AUX + AUX';
params.beta_M = beta_M; % In Mild conditions
params.beta_H = beta_M; % Hospitalized individuals
params.beta_I = beta_M; % In ICU
params.distance = distance;
PARAMS1 = unknowns;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% Estimating a time-dependent transmission parameter (beta) from day 1
% %%%% until 19:

BETA = zeros(length(t_actual),NP);
VACC = zeros(length(t_actual),NP);
BETA(1,:) = ones(1,NP);
unknowns0 = ones(1,NP);
priors = unknowns0;
LB = 1E-3*ones(1,NP);
UB = 50*ones(1,NP);
yinit2 = yinit;
yb = zeros(length(t_actual),length(yinit));
yb(1,:) = yinit;
time = zeros;
for jj =1:day-1
t_actual2 = t_actual(jj:jj+1);
OF2 = @(unknowns2)ObjFun_beta6(t_actual2,params,unknowns2,data2(jj,:),options,priors,BETA(jj,:),yinit2);
unknowns = lsqnonlin(OF2,unknowns0,LB,UB,options3);
unknowns0 = unknowns(1:NP);
priors = unknowns0;
BETA(jj+1,:) = unknowns(1:NP);
beta2 = @(t)betaA(BETA(jj:jj+1,:),t_actual2,t,NP);

[~,y2] = ode45(@(t,y)seir_death_age_beta4(t,y,params,beta2),...
                                                 t_actual2,yinit2,options);
yinit2 = y2(end,:);
yb(jj+1,:) = yinit2;
disp(num2str(unknowns))
time = time+1;
end
auxB = 1:day-1;

dt = 20;
nmax = floor(length(t_actual)/dt);
NewCases = zeros(nmax-1,NP+1);
NewHospt = zeros(nmax-1,NP+1);
NewICUt = zeros(nmax-1,NP+1);
NewDeathst = zeros(nmax-1,NP+1);
NRES = zeros(nmax-1,1);
PERR = zeros(nmax-1,1);
[NewCases(1,:),NewHospt(1,:),NewICUt(1,:),NewDeathst(1,:),NRES(1),PERR(1)] = Backtesting_20210326(t_span,t_actual,BETA,time,Hosp,Death,ICU,UnRate,EmRate,NU,yb,params,Number,NP,options,data2,Provinces,auxB);
                     
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Estimating basic parameters from day 21 until 80:

dt = 20;
nmax = floor(length(t_actual)/dt);
PARAMS2 = zeros(nmax,13+5);
for ss = 1:nmax-1 
%%% Estimating the transmission constant parameters (M,H,I), the initial
%%% infecve population (I_M0) and the transmission matrix:
if ss == nmax-1
auxB = (day-1)+(ss-1)*dt+1:length(t_actual)-1;
aux2 = (day-1)+(ss-1)*dt+1:length(t_actual);    
else
auxB = (day-1)+(ss-1)*dt+1:(day-1)+ss*dt;
aux2 = (day-1)+(ss-1)*dt+1:(day-1)+ss*dt+1;
end

pesos = ones(1,13);

LB1 = zeros(1,NP);
UB1 = 50*ones(1,NP);
unknowns10=beta*ones(1,NP);
                          
priors1 = unknowns10;
unknowns20 = [0.5,0.25,0.125,0.625,0.03];
priors2 = unknowns20;
LB2 = zeros(size(unknowns20));
UB2 = 0.5*ones(size(unknowns20));

unknowns0 = [unknowns10,unknowns20]; % initial parameters
LB = [LB1,LB2]; % lower bounds
UB = [UB1,UB2]; % upper bounds

priors = [priors1,priors2];

yinit2 = y2(end,:);

%%% Estimating the transmission constant parameters (M,H,I), the initial
%%% infecve population (I_M0) and the transmission matrix:

OF = @(unknowns)ObjFun_BetaM8(t_actual(aux2),params,...
           data2(auxB,:),options,priors,yinit2,unknowns,PropInfections,pesos);
unknowns = lsqnonlin(OF,unknowns0,LB,UB,options2);

a = unknowns(1:NP);
params.a = a;
b = params.b;
c = params.c;

dist = params.dist;
aux0 = unknowns(NP+1:end);
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
params.beta_M = beta_M; % In Mild conditions
params.beta_H = beta_M; % Hospitalized individuals
params.beta_I = beta_M; % In ICU
params.distance = distance;

PARAMS2(ss,:) = unknowns;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% Estimating a time-dependent transmission parameter (beta) from day 21
%%%% until 80:
unknowns0 = BETA(jj+1,:);
priors = unknowns0;
LB = 1E-3*ones(1,NP);
UB = 50*ones(1,NP);
for jj = auxB
t_actual2 = t_actual(jj:jj+1);
OF2 = @(unknowns2)ObjFun_beta6(t_actual2,params,unknowns2,data2(jj,:),options,priors,BETA(jj,:),yinit2);
unknowns = lsqnonlin(OF2,unknowns0,LB,UB,options3);
unknowns0 = 0.5*ones(size(unknowns));
priors = unknowns0;
BETA(jj+1,:) = unknowns;
beta2 = @(t)betaA(BETA(jj:jj+1,:),t_actual2,t,NP);
[~,y2] = ode45(@(t,y)seir_death_age_beta4(t,y,params,beta2),...
                                                 t_actual2,yinit2,options);
yinit2 = y2(end,:);
yb(jj+1,:) = yinit2;
disp(num2str(unknowns))
time = time+1;
end
if ss< nmax-1
[NewCases(ss+1,:),NewHospt(ss+1,:),NewICUt(ss+1,:),NewDeathst(ss+1,:),NRES(ss+1),PERR(ss+1)] = ...
    Backtesting_20210326(t_span,t_actual,BETA,time,Hosp,Death,ICU,UnRate,EmRate,NU,yb,params,Number,NP,options,data2,Provinces,auxB);
end
end

for ll = 1:nmax-1
if ll ==1
time = 19;
disp('CANADA')
formatOut = 'dd-mmm';
ndays = 7;
else
auxB = (day-1)+(ll-2)*dt+1:(day-1)+(ll-1)*dt;
time = time+length(auxB);
end
period = t_actual(1:time);
t_actual2 = period(1):period(end)+ndays;
t_span2 = t_span(1) + caldays(0:length(t_actual2)-1);

np = length(period);
AUX2 = [NewCases(ll,end),NewHospt(ll,end),NewICUt(ll,end),NewDeathst(ll,end)];
AUX3 = [sum(sum(data2(time-1:length(t_actual2)-1,1:13))),sum(sum(data2(time-1:length(t_actual2)-1,14:26))),...
        sum(sum(data2(time-1:length(t_actual2)-1,27:39))),sum(sum(data2(time-1:length(t_actual2)-1,40:52)))];

disp([datestr(t_span2(np+1),formatOut),' to ',datestr(t_span2(end),formatOut),';',num2str(round(AUX2(1))),';',num2str(round(AUX2(2))),';',num2str(round(AUX2(3))),';',num2str(round(AUX2(4))),';',num2str(PERR(ll)),';',num2str(NRES(ll))])
disp([' ',';',num2str(round(AUX3(1))),';',num2str(round(AUX3(2))),';',num2str(round(AUX3(3))),';',num2str(round(AUX3(4)))])
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Estimating basic parameters from day 21 until 80:
%%% Final Number of Cases for each Age Range
%%% Total Number of Deaths for each day
sigma=params.sigma;
Susceptible = zeros(size(yb,1),NP);
NewInfections = zeros(size(yb,1),NP);
NewDeaths = zeros(size(yb,1),NP);
NewHospB = zeros(size(yb,1),NP*Number);
NewHosp = zeros(size(yb,1),NP);
NewICU =  zeros(size(yb,1),NP);
Vaccinated = zeros(size(yb,1),NP);
Unemployed = zeros(size(yb,1),NP);
Employed = zeros(size(yb,1),NP);
LabourForce = zeros(size(yb,1),NP);
UNEMPL = zeros(length(t_span),NP);
EMPL = zeros(length(t_span),NP);
LABOUR = zeros(length(t_span),NP);


close all;
dates = [datetime(2020,02,19),datetime(2020,03:27,01)];

A = caldays(between(dates(1:end-1),dates(2:end),'days'));
B = [0,cumsum(A)];
coord = [B(1:end-1)'+1,B(2:end)'];

EmplByMonth = zeros(13,length(dates));
UnemplByMonth = zeros(13,length(dates));
LabourByMonth = zeros(13,length(dates));

for ii = 1:NP
auxD = min(0.5,Death(:,ii)./sum(yb(:,6*NP*Number+(ii-1)*Number+1:6*NP*Number+ii*Number)+yb(:,7*NP*Number+(ii-1)*Number+1:7*NP*Number+ii*Number),2));
auxH = min(0.5,Hosp(:,ii)./sum(yb(:,6*NP*Number+(ii-1)*Number+1:6*NP*Number+ii*Number)+yb(:,7*NP*Number+(ii-1)*Number+1:7*NP*Number+ii*Number),2));
auxI = min(0.5,ICU(:,ii)./sum(yb(:,8*NP*Number+(ii-1)*Number+1:8*NP*Number+ii*Number)+yb(:,9*NP*Number+(ii-1)*Number+1:9*NP*Number+ii*Number),2));
auxD(isnan(auxD))=zeros;
auxD(isinf(auxD))=zeros;
auxH(isnan(auxH))=zeros;
auxH(isinf(auxH))=zeros;
auxI(isnan(auxI))=zeros;
auxI(isinf(auxI))=zeros;
for jj = 1:Number
Susceptible(:,ii) = Susceptible(:,ii)+yb(:,0*(NP*Number)+(ii-1)*Number+jj)*N + sigma*yb(:,1*(NP*Number)+(ii-1)*Number+jj)*N;
NewInfections(:,ii) = NewInfections(:,ii)+sigma*yb(:,4*(NP*Number)+(ii-1)*Number+jj)*N + sigma*yb(:,5*(NP*Number)+(ii-1)*Number+jj)*N;
NewDeaths(:,ii) = NewDeaths(:,ii) + Death_I((ii-1)*Number+jj).*yb(:,10*NP*Number+jj+(ii-1)*Number) + ...
                max(0,auxD-(Death_I((ii-1)*Number+jj)*GetWorse_M((ii-1)*Number+jj).*auxH).*(GetWorse_H((ii-1)*Number+jj)*auxI)).*yb(:,6*NP*Number+jj+(ii-1)*Number)+...
                Death_I((ii-1)*Number+jj).*yb(:,11*NP*Number+jj+(ii-1)*Number) + ...
                max(0,auxD-(Death_I((ii-1)*Number+jj)*GetWorse_M((ii-1)*Number+jj).*auxH).*(GetWorse_H((ii-1)*Number+jj)*auxI)).*yb(:,7*NP*Number+jj+(ii-1)*Number);
NewHosp(:,ii) = NewHosp(:,ii)+GetWorse_M((ii-1)*Number+jj)*(auxH.*yb(:,6*(NP*Number)+(ii-1)*Number+jj))+...
                              GetWorse_M((ii-1)*Number+jj)*(auxH.*yb(:,7*(NP*Number)+(ii-1)*Number+jj));
NewICU(:,ii) = NewICU(:,ii)+GetWorse_H((ii-1)*Number+jj)*(auxI.*yb(:,8*(NP*Number)+(ii-1)*Number+jj))+...
                            GetWorse_H((ii-1)*Number+jj)*(auxI.*yb(:,9*(NP*Number)+(ii-1)*Number+jj));
Vaccinated(:,ii) = Vaccinated(:,ii)+yb(:,2*(NP*Number)+(ii-1)*Number+jj)*N+yb(:,3*(NP*Number)+(ii-1)*Number+jj)*N;
Unemployed(:,ii) = Unemployed(:,ii)+yb(:,(NP*Number)+(ii-1)*Number+jj)*N + yb(:,3*(NP*Number)+(ii-1)*Number+jj)*N + yb(:,5*(NP*Number)+(ii-1)*Number+jj)*N+...
    yb(:,7*(NP*Number)+(ii-1)*Number+jj)*N + yb(:,9*(NP*Number)+(ii-1)*Number+jj)*N+yb(:,11*(NP*Number)+(ii-1)*Number+jj)*N+yb(:,13*(NP*Number)+(ii-1)*Number+jj)*N;
if jj==2
Employed(:,ii) = Employed(:,ii)+yb(:,(ii-1)*Number+jj)*N + yb(:,2*(NP*Number)+(ii-1)*Number+jj)*N + yb(:,4*(NP*Number)+(ii-1)*Number+jj)*N+...
    yb(:,6*(NP*Number)+(ii-1)*Number+jj)*N + yb(:,8*(NP*Number)+(ii-1)*Number+jj)*N+yb(:,10*(NP*Number)+(ii-1)*Number+jj)*N+yb(:,12*(NP*Number)+(ii-1)*Number+jj)*N;
end
end
LabourForce(:,ii) = Unemployed(:,ii)+Employed(:,ii);
UNEMPL(:,ii) = N*UnemplDay(ii,194:934)';
EMPL(:,ii) = N*EmplDay(ii,194:934)';
LABOUR(:,ii) = UNEMPL(:,ii) + EMPL(:,ii);
corr = LabourForce(1,ii)/LABOUR(1,ii);
UNEMPL(:,ii) = corr*UNEMPL(:,ii);
EMPL(:,ii) = corr*EMPL(:,ii);
LABOUR(:,ii) = corr*LABOUR(:,ii);
t_spanc = t_span;
t_spanc(end) = datetime(2022,03,01);
EmplByMonth(ii,:) = interp1(t_spanc,Employed(2:end,ii)',dates);
UnemplByMonth(ii,:) = interp1(t_spanc,Unemployed(2:end,ii)',dates);
LabourByMonth(ii,:) = interp1(t_spanc,LabourForce(2:end,ii)',dates);
end
save data_20220729_REINFECTION360;

l=0;
H = [100 100 1200 600];
formatOut = 'dd-mmm';
DT = 20;

figure
hold on
grid on
box on
title('Susceptible')
plot(t_span(1:end-l),sum(Susceptible(2:end-l,:),2),'r','LineWidth',1)
ylabel('Number of Individuals')
xtickformat('dd-MMM')
set(gca,'FontSize',16,'FontName','Arial')
set(gcf,'Position',H)
hold off
saveas(gcf,'Infections_CA.fig');
print('-dpng','Infections_CA');

figure
hold on
grid on
box on
title('Infections')
bar(t_span(1:end-l),sum(data2(1:end-l,1:NP),2),'FaceColor',[0 0.75 0.75],'EdgeColor','none','FaceAlpha',0.8)
plot(t_span(1:end-l),sum(NewInfections(2:end-l,:),2),'r','LineWidth',1)
legend('Reported','Predicted','Location','NorthWest')
ylabel('Number of Individuals')
xtickformat('dd-MMM')
set(gca,'FontSize',16,'FontName','Arial')
set(gcf,'Position',H)
hold off
saveas(gcf,'Infections_CA.fig');
print('-dpng','Infections_CA');

figure
hold on
grid on
box on
title('Deaths')
bar(t_span(1:end-l),sum(data2(1:end-l,40:52),2),'FaceColor',[0 0.75 0.75],'EdgeColor','none','FaceAlpha',0.8)
plot(t_span(1:end-l),sum(NewDeaths(2:end-l,:)*N,2),'r','LineWidth',1)
legend('Reported','Predicted','Location','NorthWest')
ylabel('Number of Individuals')
xtickformat('dd-MMM')
set(gca,'FontSize',16,'FontName','Arial')
set(gcf,'Position',H)
hold off
saveas(gcf,'Deaths_CA.fig');
print('-dpng','Deaths_CA');


figure
hold on
grid on
box on
title('Hospitalizations')
bar(t_span(1:end-l),sum(data2(1:end-l,14:26),2),'FaceColor',[0 0.75 0.75],'EdgeColor','none','FaceAlpha',0.8)
plot(t_span(1:end-l),sum(NewHosp(2:end-l,:)*N,2),'r','LineWidth',1)
legend('Reported','Predicted','Location','NorthWest')
ylabel('Number of Individuals')
xtickformat('dd-MMM')
set(gca,'FontSize',16,'FontName','Arial')
set(gcf,'Position',H)
hold off
saveas(gcf,'Hosp_CA.fig');
print('-dpng','Hosp_CA');


figure
hold on
grid on
box on
title('ICU Admissions')
bar(t_span(1:end-l),sum(data2(1:end-l,27:39),2),'FaceColor',[0 0.75 0.75],'EdgeColor','none','FaceAlpha',0.8)
plot(t_span(1:end-l),sum(NewICU(2:end-l,:)*N,2),'r','LineWidth',1)
legend('Reported','Predicted','Location','NorthWest')
ylabel('Number of Individuals')
xtickformat('dd-MMM')
set(gca,'FontSize',16,'FontName','Arial')
set(gcf,'Position',H)
hold off
saveas(gcf,'ICU_CA.fig');
print('-dpng','ICU_CA');

Vacc1 = cumsum(data2(:,53:end),1);
figure
hold on
grid on
box on
title('Vaccinated')
bar(t_span,0.7*sum(Vacc1,2),'FaceColor',[0 0.75 0.75],'EdgeColor','none','FaceAlpha',0.8)
plot(t_span,sum(Vaccinated(2:end,:),2),'r','LineWidth',1)
legend('Reported','Predicted','Location','NorthWest')
ylabel('Number of Individuals')
xtickformat('dd-MMM')
set(gca,'FontSize',16,'FontName','Arial')
set(gcf,'Position',H)
hold off
saveas(gcf,'Vacc_CA.fig');
print('-dpng','Vacc_CA');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% Labour Market Data
DATA = importdata('LabourCanada2010_2022.xlsx');
DATA = DATA.data; NP=13;

WorkForce = DATA(2:6:end,:);
Employment = DATA(3:6:end,:);
Unemployment = DATA(6:6:end,:);
PopByProvince = DATA(1:6:end,:);
Employment = Employment(:,123:148);
Unemployment = Unemployment(:,123:148);
WorkForce = WorkForce(:,123:148);
for jj=1:NP
EmplByMonth(jj,:) = (WorkForce(jj,:)/LabourByMonth(jj,1)).*EmplByMonth(jj,:);
UnemplByMonth(jj,:) = (WorkForce(jj,:)/LabourByMonth(jj,1)).*UnemplByMonth(jj,:);
end


figure
hold on
box on
grid on
title('Employed Individuals by 1000')
bar(dates,sum(Employment),'FaceColor',[0 0.75 0.75],'EdgeColor','none','FaceAlpha',0.8)
plot(dates,sum(EmplByMonth),'r','LineWidth',2)
legend('Stats Canada','Model')
xtickformat('MMM-yy')
set(gca,'FontSize',16,'FontName','Arial')
set(gcf,'Position',H)
hold off
print('-dpng','Empl_CA');

figure
hold on
box on
grid on
title('Unemployed Individuals by 1000')
bar(dates,sum(Unemployment),'FaceColor',[0 0.75 0.75],'EdgeColor','none','FaceAlpha',0.8)
plot(dates,sum(UnemplByMonth),'r','LineWidth',2)
legend('Stats Canada','Model')
xtickformat('MMM-yy')
set(gca,'FontSize',16,'FontName','Arial')
set(gcf,'Position',H)
hold off
print('-dpng','Unempl_CA');

