function [dxL,Lambda]=ODE_Equations(t,xt,t0,a0L,Ncomp,VEIReinfection,delta,nu,alpha,mu,EffectivnessReinfection,WanNaturalImmu,NR)

a1=a0L(1);
a2=a0L(2);

b1=a0L(3);
b2=b1+a0L(4);

c1=a0L(5);
c2=a0L(6);

beta=(a1/(1+exp((t-b1)/c1)))+a2/(1+exp(-c2*(t-b2)));


[Lambda]=ForceofInfection(xt,beta,Ncomp,VEIReinfection,NR);

dx=zeros(Ncomp,1);

dx(1)= mu*sum(xt(:))+alpha.*(xt(3)+sum(xt(6:4:6+4*(NR-1))))-Lambda*xt(1)-mu*xt(1);             %(S)
dx(2)= Lambda*xt(1) - delta*xt(2)-mu*xt(2);       %(E)
dx(3)=delta*xt(2) - nu*xt(3)-alpha*xt(3)-mu*xt(3);      %(I)

dx(4)= nu*xt(3)-(1-EffectivnessReinfection(1))*Lambda*xt(4)-mu*xt(4)-WanNaturalImmu(1)*xt(4);          %(R1)
dx(5)= (1-EffectivnessReinfection(1))*Lambda*(xt(4)+xt(7)) - delta*xt(5)-mu*xt(5)-WanNaturalImmu(1)*xt(5);       %(E1)
dx(6)=delta*xt(5) - nu*xt(6)-alpha*xt(6)-mu*xt(6)-WanNaturalImmu(1)*xt(6);      %(I1)
dx(7)= nu*sum(xt(6:4:6+4*(NR-1)))-mu*xt(7)-WanNaturalImmu(1)*xt(7)-(1-EffectivnessReinfection(1))*Lambda*xt(7);       %(R1)

for n=2:NR
    dx(4+4*(n-1))= WanNaturalImmu(n-1)*xt(4+4*(n-2))-(1-EffectivnessReinfection(n))*Lambda*xt(4+4*(n-1))-mu*xt(4+4*(n-1))-WanNaturalImmu(n)*xt(4+4*(n-1));         %(Rn)
    dx(5+4*(n-1))= WanNaturalImmu(n-1)*xt(5+4*(n-2))+(1-EffectivnessReinfection(n))*Lambda*(xt(4+4*(n-1))+xt(7+4*(n-1))) - delta*xt(5+4*(n-1))-mu*xt(5+4*(n-1))-WanNaturalImmu(n)*xt(5+4*(n-1));       %(En)
    dx(6+4*(n-1))= WanNaturalImmu(n-1)*xt(6+4*(n-2))+delta*xt(5+4*(n-1)) - nu*xt(6+4*(n-1))-alpha*xt(6+4*(n-1))-mu*xt(6+4*(n-1))-WanNaturalImmu(n)*xt(6+4*(n-1));      %(In)
    dx(7+4*(n-1))= WanNaturalImmu(n-1)*xt(7+4*(n-2))-mu*xt(7+4*(n-1))-WanNaturalImmu(n)*xt(7+4*(n-1))-(1-EffectivnessReinfection(n))*Lambda*xt(7+4*(n-1));       %(Rn)
end
%+nu*xt(6+4*(n-1))    
dxL=dx;
end



