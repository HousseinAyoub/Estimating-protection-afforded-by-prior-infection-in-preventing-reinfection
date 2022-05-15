function [Lambda]=ForceofInfection(x,beta,Ncomp,VEIReinfection)

eps=1e-10;

rhoN=0;
for s=1:Ncomp %Stage
    rhoN=rhoN+x(s);
end
Lambda=beta.*(x(3)+(1-VEIReinfection)*x(6))/(rhoN+eps);


