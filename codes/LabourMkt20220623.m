function yprime = LabourMkt20220623(t,y, params,Infections)

rhoE = params.rhoE;
rhoU = params.rhoU;
lambdaE = params.lambdaE;
lambdaU = params.lambdaU;


Em = y(1);
Un = y(2);
yprime = zeros(2,1);
% yprime(1) = - (lambdaU + rhoU*Infections(t))*Em + (lambdaE+rhoE*Infections(t))*Un;
% yprime(2) = - (lambdaE+rhoE*Infections(t))*Un + (lambdaU+ rhoU*Infections(t))*Em;
yprime(1) = - (lambdaU + rhoU*Infections(t))*Em + (lambdaE+rhoE*Infections(t))*Un;
yprime(2) = - (lambdaE+rhoE*Infections(t))*Un + (lambdaU+ rhoU*Infections(t))*Em;
