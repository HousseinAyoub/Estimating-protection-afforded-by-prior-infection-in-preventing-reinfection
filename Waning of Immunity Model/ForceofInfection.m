function [Lambda]=ForceofInfection(x,beta,Ncomp,VEIReinfection,NR)

eps=1e-10;

rhoN=0;
for s=1:Ncomp %Stage
    rhoN=rhoN+x(s);
end
Lambda=beta.*(x(3)+(1-VEIReinfection)*sum(x(6:4:6+4*(NR-1))))/(rhoN+eps);


