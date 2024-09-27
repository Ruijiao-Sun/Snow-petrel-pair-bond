%Created by Hal Caswell, last modified by Silke van Daalen, May 2019

function out=rewards(P,R1,R2,R3,R4)
%function to compute mean, variance, and skewness of reward
% P = markov transition matrix
% R1, R2, R3, R4 = matrices of first four raw moments of reward


%size of matrices
[s,s]=size(P);
e=ones(s,1);

% Creating what we need for the equilibrium calculation
% First define Z
Z= [eye(s-1,s-1) zeros(s-1,1)];

% And I:
I= [eye(s-1,s-1)];

% and u:
u= Z*P*Z.';

%Don't forget our vector of ones:
e=ones(s,1);
% equilibria for the moments

rho1_eq = inv(I-u.')*(Z*(((P.*R1).')*e));
rho2_eq = inv(I-u.')*(Z*(((P.*R2).')*e)+2*Z*((P.*R1).')*(Z.')*rho1_eq);
rho3_eq = inv(I-u.')*(Z*(((P.*R3).')*e)+3*Z*((P.*R2).')*(Z.')*rho1_eq + ...
    3*Z*((P.*R1).')*(Z.')*rho2_eq);
rho4_eq = inv(I-u.')*(Z*(((P.*R4).')*e)+4*Z*((P.*R3).')*(Z.')*rho1_eq + ...
    6*Z*((P.*R2).')*(Z.')*rho2_eq +4*Z*((P.*R1).')*(Z.')*rho3_eq);

%calculate some statistics from the moment vectors
var_eq=rho2_eq - rho1_eq.^2;
std_eq=sqrt(var_eq);
cv_eq=std_eq./rho1_eq;
pick=find(isnan(cv_eq));
cv_eq(pick)=0;
crow_eq=var_eq./(rho1_eq.^2);
pick=find(isnan(crow_eq));
crow_eq(pick)=0;

vv=var_eq(:,end);
pick=find(vv==0);
vv(pick)=1;
vv=vv.^(-3/2);
vv(pick)=0;
dd=diag(vv);

skew_eq=dd*...
    (rho3_eq(:,end) - 3*rho2_eq(:,end).*rho1_eq(:,end) + 2*(rho1_eq(:,end).^3));

kurt_eq=(-3*rho1_eq.^4+6*rho2_eq.*rho1_eq.^2-4*rho1_eq.*rho3_eq+rho4_eq)./(var_eq.^2)-3;

% k-fold differences (lifetime)
% k=2
k2=2*normcdf(-log(2)/sqrt(2*log(cv_eq(1,:)^2+1)));
% k=10
k10=2*normcdf(-log(10)/sqrt(2*log(cv_eq(1,:)^2+1)));

% some longevity things:
N_out=inv(I-u);
ee=ones(1,s-1);
eta_out=ee*inv(I-u);
vareta_out=ee*N_out*(2*N_out-I)-(eta_out.*eta_out);

%create structure with output in it
out.rho1=rho1_eq;
out.rho2=rho2_eq;
out.rho3=rho3_eq;
out.rho4=rho4_eq;
out.var=var_eq;
out.std=std_eq;
out.cv=cv_eq;
out.crow=crow_eq;
out.skew=skew_eq;
out.kurt=kurt_eq;
out.k2=k2;
out.k10=k10;
out.eta=eta_out;
out.vareta=vareta_out; 
out.P=P;
out.R1=R1;
out.R2=R2;
out.R3=R3;






