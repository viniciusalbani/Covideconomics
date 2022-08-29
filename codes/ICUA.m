function f = ICUA(ICU,TA,t,NP)
f = zeros(length(t),NP);
for jj = 1:NP
f(:,jj) = interp1(TA,ICU(:,jj),t);
end