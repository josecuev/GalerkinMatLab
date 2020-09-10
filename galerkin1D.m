%% Garlekin Method
%  1D Gerneralization - Error norms
% 
%% Problem data definition

    syms x;
    
    k=1;
    c=1;
    f=0;
    L=1;
    u=2*x+3;


%% Function definitions

% Strain energy in one dimension
function[energy]= strainenergy(u,k,c,L)
syms x
energy=0.5*int(k*diff(u,x)^2+c*u^2,x,0,L)
end