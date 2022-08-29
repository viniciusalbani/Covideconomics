%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% Reading Daily Cases and Deaths in Canada

DATA = importdata('Provincial_Daily_Totals_20220316.txt');
data = DATA.data(1:end,:);
text = DATA.textdata;

AB = []; BC = []; MN = []; NB = []; NL = []; NT = []; NS = []; NU = [];
ON = []; PE = []; QB = []; SK = []; YK = [];

dtALB = []; dtBC = []; dtMN = []; dtNB = []; dtNL = []; dtNT = []; 
dtNS = []; dtNU = []; dtON = []; dtPEI = []; dtQB = []; dtSK = []; 
dtYK = [];

data(isnan(data)==1)=zeros; 

ii = ones;
aux = string(text(:,2));
Provinces = aux;
aux2 = datetime(text(2:end,3), 'InputFormat', 'dd/MM/yyyy')';
% aux2 = datetime(text(2:end,3), 'InputFormat', 'yyyy/MM/dd')';
t_span = min(aux2):max(aux2);
% [Cases,Hospitalizations,ICU Adm.,Deaths,Vaccinated] (daily).
H = [2,12,14,6,20];
cases = zeros(length(t_span),13*length(H));

while ii<size(data,1)
if aux(ii+1) == "AB"
AB = [AB;data(ii,H)];
elseif aux(ii+1) == "BC"
BC = [BC;data(ii,H)];
elseif aux(ii+1) == "MB"
MN = [MN;data(ii,H)];
elseif aux(ii+1) == "NB"
NB = [NB;data(ii,H)];
elseif aux(ii+1) == "NL"
NL = [NL;data(ii,H)];
elseif aux(ii+1) == "NT"
NT = [NT;data(ii,H)];
elseif aux(ii+1) == "NS"
NS = [NS;data(ii,H)];
elseif aux(ii+1) == "NU"
NU = [NU;data(ii,H)];
elseif aux(ii+1) == "ON"
ON = [ON;data(ii,H)];
elseif aux(ii+1) == "PE"
PE = [PE;data(ii,H)];
elseif aux(ii+1) == "QC"
QB = [QB;data(ii,H)];
elseif aux(ii+1) == "SK"
SK = [SK;data(ii,H)];
elseif aux(ii+1) == "YT"
YK = [YK;data(ii,H)];
end
ii=ii+1;
end
H = [1,14,27,40,53];
%%% time

cases(:,H) = AB;
cases(:,H+1) = BC;
cases(:,H+2) = MN;
cases(:,H+3) = NB;
cases(:,H+4) = NL;
cases(:,H+5) = NT;
cases(:,H+6) = NS;
cases(:,H+7) = NU;
cases(:,H+8) = ON;
cases(:,H+9) = PE;
cases(:,H+10) = QB;
cases(:,H+11) = SK;
cases(:,H+12) = YK;

%cases(2:end,:) = diff(cases);
cases = max(cases(1:end-15,:),0);
t_span = t_span(1:end-15);
t_actual = 0:length(t_span);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% Reading Distance Information
DATA = importdata('DistanceBerweenStatesCAN.xlsx');
data = DATA.data;
data(isnan(data)==1)=zeros;
data = data+data';
dist = data;
aux = data;
aux(data==0)=ones;
aux((data>=1 & data<1000))=0.5*ones;
aux((data>=1000 & data<2000))=0.25*ones;
aux((data>=2000 & data<=3000))=0.125*ones;
aux((data>=3000 & data<=4000))=0.0625*ones;
aux(data>=4000)= 0.03*ones;
distance = aux;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% Reading Demographic Information
DATA = importdata('Demographics.xlsx');
data = DATA.data;
text = DATA.textdata;
demo = data(1,:);

PropPopAge = [1080018	991471	270420	152797	96958	9138	144213	15742	3131022	25347	561852	158425	9234
2940957	2867933	881460	411086	322616	31157	578734	22738	9812596	93816	5182878	747802	26847
400901	671160	147976	115785	76510	2157	140412	873	1790396	21841	1162204	100230	3226];

%%% Source:  How to cite: Statistics Canada. Table 17-10-0005-01 Population estimates on July 1st, by age and sex
%%% DOI: https://doi.org/10.25318/1710000501-eng
