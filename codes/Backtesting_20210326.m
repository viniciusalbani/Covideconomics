function [NewCasestt,NewHosptt,NewICUtt,NewDeathstt,NRES,PERR] = Backtesting_20210326(t_span,t_actual,BETA,time,Hosp,Death,ICU,UnRate,EmRate,NU,yb,params,Number,NP,options,data2,Provinces,auxB)
ndays = 7;
formatOut = 'dd-mmm';
l = 10;
%%%%
period = t_actual(1:time);
t_actual2 = period(1):period(end)+ndays;
t_span = t_span(1) + caldays(0:length(t_actual2)-1);

np = length(period);
% nt = length(t_actual2(np:end));
%%%%%%%% Case 1
BETA2 = [BETA(1:np,:);ones(ndays,1)*mean(BETA(np-round(l/3):np,:))];
beta2 = @(t)betaA(BETA2,t_actual2,t,NP);
Hosp2 = [Hosp(1:np,:);ones(ndays,1)*mean(Hosp(np-l:np,:))];
Death2 = [Death(1:np,:);ones(ndays,1)*mean(Death(np-l:np,:))];
ICU2 = [ICU(1:np,:);ones(ndays,1)*mean(ICU(np-l:np,:))];
params.factorWorse = @(t)hospA(Hosp2,t_actual2,t,NP);
params.factorICU = @(t)ICUA(ICU2,t_actual2,t,NP);
params.factorDeath = @(t)deathA(Death2,t_actual2,t,NP);
UnRate2 = [UnRate(1:np,:);ones(ndays,1)*mean(UnRate(np-l:np,:))];
EmRate2 = [EmRate(1:np,:);ones(ndays,1)*mean(EmRate(np-l:np,:))];
NU2 = [NU(1:np,:);ones(ndays,1)*mean(NU(np-l:np,:))];
params.lambda_Un = @(t)hospA(UnRate2,t_actual2,t,NP); % Unemployment rate
params.lambda_Em = @(t)hospA(EmRate2,t_actual2,t,NP); % Employment rate
params.nu = @(t)hospA(NU2,t_actual2,t,NP); % Vaccination rate



sigma = params.sigma;
N = params.N;
GetWorse_M = params.GetWorse_M;
GetWorse_H = params.GetWorse_H;
Death_I = params.Death_I;
yb2 = zeros(length(t_actual2),15*Number*NP);
yb2(1:time,:) = yb(1:time,:);
yinit2 = yb(time,:);
for jj = time:length(t_actual2)-1
t_actualB = t_actual2(jj:jj+1);
[~,y2] = ode45(@(t,y)seir_death_age_beta4(t,y,params,beta2),...
                                                 t_actualB,yinit2,options);
yinit2 = y2(end,:);
yb2(jj+1,:) = yinit2;
end

NewCases = zeros(size(yb2,1),NP);
NewDeaths = zeros(size(yb2,1),NP);
NewHosp = zeros(size(yb2,1),NP);
NewICU = zeros(size(yb2,1),NP);

disp([datestr(t_span(np+1),formatOut),' to ',datestr(t_span(end))])

% AUX2 = zeros(1,4); AUX3 = zeros(1,4);
for ii = 1:NP
auxD = min(0.5,Death2(:,ii)./sum(yb2(:,6*NP*Number+(ii-1)*Number+1:6*NP*Number+ii*Number)+yb2(:,7*NP*Number+(ii-1)*Number+1:7*NP*Number+ii*Number),2));
auxH = min(0.5,Hosp2(:,ii)./sum(yb2(:,6*NP*Number+(ii-1)*Number+1:6*NP*Number+ii*Number)+yb2(:,7*NP*Number+(ii-1)*Number+1:7*NP*Number+ii*Number),2));
auxI = min(0.5,ICU2(:,ii)./sum(yb2(:,8*NP*Number+(ii-1)*Number+1:8*NP*Number+ii*Number)+yb2(:,9*NP*Number+(ii-1)*Number+1:9*NP*Number+ii*Number),2));
auxD(isnan(auxD))=zeros;
auxD(isinf(auxD))=zeros;
auxH(isnan(auxH))=zeros;
auxH(isinf(auxH))=zeros;
auxI(isnan(auxI))=zeros;
auxI(isinf(auxI))=zeros;
for jj = 1:Number
NewCases(:,ii) = NewCases(:,ii) + sigma*yb2(:,4*NP*Number+jj+(ii-1)*Number)*N + sigma*yb2(:,5*NP*Number+jj+(ii-1)*Number)*N;
NewDeaths(:,ii) = NewDeaths(:,ii) + Death_I((ii-1)*Number+jj).*yb2(:,10*NP*Number+jj+(ii-1)*Number) + ...
                                    max(0,auxD-(Death_I((ii-1)*Number+jj)*GetWorse_M((ii-1)*Number+jj).*auxH).*(GetWorse_H((ii-1)*Number+jj)*auxI)).*yb2(:,6*NP*Number+jj+(ii-1)*Number)+...
                                    Death_I((ii-1)*Number+jj).*yb2(:,11*NP*Number+jj+(ii-1)*Number) + ...
                                    max(0,auxD-(Death_I((ii-1)*Number+jj)*GetWorse_M((ii-1)*Number+jj).*auxH).*(GetWorse_H((ii-1)*Number+jj)*auxI)).*yb2(:,7*NP*Number+jj+(ii-1)*Number);
NewHosp(:,ii) = NewHosp(:,ii) + GetWorse_M((ii-1)*Number+jj)*(auxD.*yb2(:,6*NP*Number+jj+(ii-1)*Number))+...
                                GetWorse_M((ii-1)*Number+jj)*(auxD.*yb2(:,7*NP*Number+jj+(ii-1)*Number));
NewICU(:,ii) = NewICU(:,ii) + GetWorse_H((ii-1)*Number+jj)*(auxI.*yb2(:,8*NP*Number+jj+(ii-1)*Number))+...
                              GetWorse_H((ii-1)*Number+jj)*(auxI.*yb2(:,9*NP*Number+jj+(ii-1)*Number));
end 
NewDeaths(:,ii) = NewDeaths(:,ii)*N;
NewHosp(:,ii) = NewHosp(:,ii)*N;
NewICU(:,ii) = NewICU(:,ii)*N;
end
NewCasestt = zeros(1,NP+1);
NewHosptt = zeros(1,NP+1);
NewICUtt = zeros(1,NP+1);
NewDeathstt = zeros(1,NP+1);

NewCasest = sum(NewCases(time+1:end,:))';
NewHospt = sum(NewHosp(time+1:end,:))';
NewICUt = sum(NewICU(time+1:end,:))';
NewDeathst = sum(NewDeaths(time+1:end,:))';

NewCasestt(1:NP) = NewCasest;
NewHosptt(1:NP) = NewHospt;
NewICUtt(1:NP) = NewICUt;
NewDeathstt(1:NP) = NewDeathst;

AUX2 = [NewCasest,NewHospt,NewICUt,NewDeathst];
AUX3 = [sum(data2(time-1:length(t_actual2)-1,1:13))',sum(data2(time-1:length(t_actual2)-1,14:26))',...
    sum(data2(time-1:length(t_actual2)-1,27:39))',sum(data2(time-1:length(t_actual2)-1,40:52))'];
NRES = zeros(1,NP);
PERR = zeros(1,NP);
for ll = 1:length(NewCasest)
NRES(ll)= norm(data2(auxB,ll) - NewCases(auxB+1,ll))/max(1E-1,norm(data2(auxB,ll)));
PERR(ll) = abs(AUX2(ll,1)-AUX3(ll,1))/max(1E-1,AUX3(ll,1));
% disp(Provinces(ll+1))
% disp(['predicted: ',num2str(round(AUX2(ll,:))),', CALIB ERROR: ',num2str(NRES(ll))])
% disp(['observed:  ',num2str(round(AUX3(ll,:))),', PREDIC ERROR: ',num2str(PERR(ll))])
end

NewCasest = sum(sum(NewCases(time+1:end,:)));
NewHospt = sum(sum(NewHosp(time+1:end,:)))';
NewICUt = sum(sum(NewICU(time+1:end,:)));
NewDeathst = sum(sum(NewDeaths(time+1:end,:)));

NewCasestt(end) = NewCasest;
NewHosptt(end) = NewHospt;
NewICUtt(end) = NewICUt;
NewDeathstt(end) = NewDeathst;

AUX2 = [NewCasest,NewHospt,NewICUt,NewDeathst];
AUX3 = [sum(sum(data2(time-1:length(t_actual2)-1,1:13))),sum(sum(data2(time-1:length(t_actual2)-1,14:26))),...
    sum(sum(data2(time-1:length(t_actual2)-1,27:39))),sum(sum(data2(time-1:length(t_actual2)-1,40:end)))];
disp('CANADA')
NRES = norm(data2(auxB,1:13) - NewCases(auxB+1,:))/norm(data2(auxB,:));
PERR = abs(AUX2(1) - AUX3(1))/AUX3(1);
disp(['predicted: ',num2str(round(AUX2)),', CALIB ERROR : ',num2str(NRES)])
disp(['observed:  ',num2str(round(AUX3)),' PREDIC ERROR: ',num2str(PERR(1))])
