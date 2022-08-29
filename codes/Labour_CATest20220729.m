clear all; close all; clc;
tic;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% Demographic and Epidemiological Data
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

DaysByMonth2020 = [11,31,30,31,30,31,31,30,31,30,31];
DaysByMonth2021 = [31,28,31,30,31,30,31,31,30,31,30,31];
DaysByMonth2022 = [31,28];%,31,30,31,30,31,31,30,31,30,31];

aux = [0,DaysByMonth2020,DaysByMonth2021,DaysByMonth2022];
aux = cumsum(aux);
coord = [aux(1:end-1)'+ones,aux(2:end)'];
CasesByMonth = zeros(13,size(coord,1));
for jj = 1:13
    for ii=1:size(coord,1)
    CasesByMonth(jj,ii)=sum(data(coord(ii,1):coord(ii,2),jj));
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% Labour Market Data
DATA = importdata('LabourCanada2010_2022.xlsx');
DATA = DATA.data; NP=13;

WorkForce = DATA(2:6:end,:);
Employment = DATA(3:6:end,:);%./WorkForce;
Unemployment = DATA(6:6:end,:);%./WorkForce;
PopByProvince = DATA(1:6:end,:);
Employment = Employment./PopByProvince;
Unemployment = Unemployment./PopByProvince;
LF = zeros(NP,size(Employment,2));
for jj=1:NP
LF(jj,:) = (PropPopAge(2,jj)/sum(demo))./(Employment(1,:)+Unemployment(1,:));
end

aux = CasesByMonth;
CasesByMonth = zeros(size(WorkForce));
CasesByMonth(:,123:147)=aux;

dates = datetime(2009,12:12+150-1,01);

A = caldays(between(dates(1:end-1),dates(2:end),'days'));
B = [0,cumsum(A)];
C = 0:B(end);

EmplDay = zeros(size(Employment,1),length(C));
UnemplDay = zeros(size(Unemployment,1),length(C));
for jj = 1:size(Employment,1)
EmplDay(jj,:) = interp1(B,LF(jj,:).*Employment(jj,:),C);
UnemplDay(jj,:) = interp1(B,LF(jj,:).*Unemployment(jj,:),C);
end
datesB = datetime(2009,12,1+C);

data2 = cases;
for jj=4:size(cases,1)-3
for ii = 1:size(cases,2)
data2(jj,ii) = mean(cases(jj-3:jj+3,ii));
end
end
cases = data2;

data2 = cases;
for jj=4:size(cases,1)-3
for ii = 1:size(cases,2)
data2(jj,ii) = mean(cases(jj-3:jj+3,ii));
end
end
cases = data2;

CASES = zeros(length(B),size(cases,2));
CASES(3733:3732+size(cases,1),:)=cases;
CASES = CASES';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
EmplDay = EmplDay(:,end-1000+1:end);
UnemplDay = UnemplDay(:,end-1000+1:end);
CASES = CASES(:,end-1000+1:end);
datesB = datesB(:,end-1000+1:end);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% Model Calibration and Simulation
tol = 1.e-6;
options = odeset('AbsTol', tol,'RelTol',tol,'MaxOrder',5,'Stats',...
                                                         'off','Refine',1);
N = 4; Len = size(EmplDay,2)-N+1;
ParamsCA = zeros(Len,4,13);
EmCA = zeros(13,size(EmplDay,2));
UnCA = zeros(13,size(UnemplDay,2));
params.rho=zeros;
parfor ii = 1:13
params2 = params;
yinit = zeros(2,1); data = zeros(2,N);
unknowns0 = [0,0,1,1];
LB = -1000*[0,0,0,0];%zeros(1,4);
UB = 1000*[1,1,1,1];
Params = zeros(Len,4);
Em = zeros(1,size(EmplDay,2));
Un = zeros(1,size(UnemplDay,2));
Em(1) = EmplDay(ii,1);
Un(1) = UnemplDay(ii,1);
prior = unknowns0;

for jj = 1:Len
yinit(1) = EmplDay(ii,jj);
yinit(2) = UnemplDay(ii,jj);
data(1,:) = EmplDay(ii,jj:jj+N-1);
data(2,:) = UnemplDay(ii,jj:jj+N-1);
Cases = CASES(ii,jj:jj+N-1);
t_actual = jj-1:jj-1+size(data,2)-1;
% if jj<198
ObjFun = @(unknowns)objFunLabourMkt20220623(t_actual,unknowns,data,options,yinit,Cases,prior);
unknowns = lsqnonlin(ObjFun,unknowns0,LB,UB);
prior=unknowns;
Params(jj,:) = unknowns;
params2.rhoE = unknowns(1);
params2.rhoU = unknowns(2);
params2.lambdaE = unknowns(3);
params2.lambdaU = unknowns(4);
tspan = [t_actual(1),t_actual(end)];
Infections = @(t)interp1(t_actual,Cases,t);
[t,y] = ode45(@(t,y)LabourMkt20220623(t,y, params2,Infections),tspan,yinit,options);
Em(jj+1:jj+N-1) = interp1(t,y(:,1),t_actual(2:end));
Un(jj+1:jj+N-1) = interp1(t,y(:,2),t_actual(2:end));
end
ParamsCA(:,:,ii) = Params;
EmCA(ii,:) = Em;
UnCA(ii,:) = Un;
end

elapsedtime = toc;
disp(['Elapsed time: ',num2str(elapsed_time),' hours.'])
