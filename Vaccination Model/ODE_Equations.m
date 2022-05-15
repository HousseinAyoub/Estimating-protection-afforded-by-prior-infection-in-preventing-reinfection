function [dxL,Lambda]=ODE_Equations(t,xt,t0,a0L,Ncomp,VEIReinfection,delta,nu,alpha,mu,EffectivnessReinfection,VaccinatedRate,VES,WaningVaccinated)

%beta=R0*(delta+mu)*(nu+mu+alpha)/delta;

a1=a0L(1);
a2=a0L(2);

b1=a0L(3);
b2=b1+a0L(4);

c1=a0L(5);
c2=a0L(6);

beta=(a1/(1+exp((t-b1)/c1)))+a2/(1+exp(-c2*(t-b2)));


[Lambda]=ForceofInfection(xt,beta,Ncomp,VEIReinfection);

dx=zeros(Ncomp,1);
%%Not vaccinated
dx(1)= WaningVaccinated*xt(7+1)+mu*sum(xt(:))+alpha*(xt(3)+xt(6)+xt(7+3)+xt(7+6))-Lambda*xt(1)-mu*xt(1)-VaccinatedRate(round((t-t0)*2+1))*xt(1);             %(S)
dx(2)= Lambda*xt(1) - delta*xt(2)-mu*xt(2);       %(E)
dx(3)=delta*xt(2) - nu*xt(3)-alpha*xt(3)-mu*xt(3);      %(I)

dx(4)= WaningVaccinated*xt(7+4)+nu*xt(3)-(1-EffectivnessReinfection(1))*Lambda*xt(4)-mu*xt(4)-VaccinatedRate(round((t-t0)*2+1))*xt(4);          %(R1)
dx(5)= (1-EffectivnessReinfection(1))*Lambda*(xt(4)+xt(7)) - delta*xt(5)-mu*xt(5);       %(E1)
dx(6)=delta*xt(5) - nu*xt(6)-alpha*xt(6)-mu*xt(6);      %(I1)
dx(7)= WaningVaccinated*xt(7+7)+nu*xt(6)-mu*xt(7)-(1-EffectivnessReinfection(1))*Lambda*xt(7)-VaccinatedRate(round((t-t0)*2+1))*xt(7);       %(R1)

%%Vaccinated
dx(7+1)= VaccinatedRate(round((t-t0)*2+1))*xt(1)-(1-VES).*Lambda*xt(7+1)-mu*xt(7+1)-WaningVaccinated*xt(7+1);             %(S)
dx(7+2)= (1-VES).*Lambda*xt(7+1) - delta*xt(7+2)-mu*xt(7+2);       %(E)
dx(7+3)=delta*xt(7+2) - nu*xt(7+3)-alpha*xt(7+3)-mu*xt(7+3);      %(I)

dx(7+4)= VaccinatedRate(round((t-t0)*2+1))*xt(4)+nu*xt(7+3)-(1-VES).*(1-EffectivnessReinfection(1))*Lambda*xt(7+4)-mu*xt(7+4)-WaningVaccinated*xt(7+4);          %(R1)
dx(7+5)= (1-VES).*(1-EffectivnessReinfection(1))*Lambda*(xt(7+4)+xt(7+7)) - delta*xt(7+5)-mu*xt(7+5);       %(E1)
dx(7+6)=delta*xt(7+5) - nu*xt(7+6)-alpha*xt(7+6)-mu*xt(7+6);      %(I1)
dx(7+7)= VaccinatedRate(round((t-t0)*2+1))*xt(7)+nu*xt(7+6)-mu*xt(7+7)-(1-VES).*(1-EffectivnessReinfection(1))*Lambda*xt(7+7)-WaningVaccinated*xt(7+7);       %(R1)
 
dxL=dx;
end



