%% Simulation the model using the estimated parameters
clc
clear all
close all

load('EstimatedParamQatarNew.mat','EstimatedParamQatarNew') %%Estimated from Qatar Model
zz=EstimatedParamQatarNew; 
a0L(1)=zz(1);
a0L(2)=zz(2);

a0L(3)=zz(3);
a0L(4)=zz(4);

a0L(5)=zz(5);
a0L(6)=zz(6);

eps=1.0e-10;

NR=24;

Ncomp=6+4*(NR-1)+1;
mu=1/(80.7*365);               % 1/life expectancy

delta=1/3.69;  % Duration of stay in E is in average 3.69 days
nu=1/3.48;    % Duration of stay in I is in average 3.48 days
alpha=1.85/10000;

t0=110;
dt=0.5; %%for the stability of the scheme
tf=t0+NR*30; %Final date
tspan=t0:dt:tf;  % Timespan to use in the ode function

%Initial conditions
N0=2.76e+06; %% Total population of Qatar
E0=1.2715e+03; %%As per Qatar Model
I0=960; %%As per Qatar Model
x0L=zeros(Ncomp,1);
x0L(2)=E0;
x0L(3)=I0;
x0L(1)=N0-x0L(3)-x0L(2);
x0=x0L;

EffectivnessRecentInfection=0.8;
for n=1:NR
EffectivnessReinfection(n)=EffectivnessRecentInfection-(n-1)*(EffectivnessRecentInfection/(NR-1));
end
plot(EffectivnessReinfection)

VEIReinfection=0;

WanNaturalImmu=(1/30).*ones(NR,1);
WanNaturalImmu(NR)=0;

%%Deterministic model
[T,x] = ode45(@(t,x)ODE_Equations(t,x,t0,a0L,Ncomp,VEIReinfection,delta,nu,alpha,mu,EffectivnessReinfection,WanNaturalImmu,NR),tspan,x0);
TT=length(T);
for i=1:TT 
     [dxL,Lambda]=ODE_Equations(T(i),x(i,:),t0,a0L,Ncomp,VEIReinfection,delta,nu,alpha,mu,EffectivnessReinfection,WanNaturalImmu,NR);
     lmbL(i)=Lambda;
end

TT=length(T);
Susceptible=x(:,1);
Latent=x(:,2);
Infected=x(:,3);

for n=1:NR
Recovered(:,n)=x(:,4+4*(n-1));
LatentR(:,n)=x(:,5+4*(n-1));
InfectedR(:,n)=x(:,6+4*(n-1));
RecoveredR(:,n)=x(:,7+4*(n-1));
end


Total=sum(x,2);
for t=1:TT
     IncidenceTotal(t)=delta.*(Latent(t)+sum(LatentR(t,:)));
    for n=1:NR        
        NumeratorRealEffectivness(t,n)=EffectivnessReinfection(n)*(Recovered(t,n)+LatentR(t,n)+InfectedR(t,n)+RecoveredR(t,n));
        DenominatorRealEffectivness(t,n)=(Recovered(t,n)+LatentR(t,n)+InfectedR(t,n)+RecoveredR(t,n));
    end  
end
for t=t0:tf-1
    IncidenceTotalT(t-t0+1)=trapz(IncidenceTotal((t-t0)*2+1:(t-t0+dt)*2+1));
    
    SusceptibleT(t-t0+1)=trapz(Susceptible((t-t0)*2+1:(t-t0+dt)*2+1));
    LatentT(t-t0+1)=trapz(Latent((t-t0)*2+1:(t-t0+dt)*2+1));
    InfectedT(t-t0+1)=trapz(Infected((t-t0)*2+1:(t-t0+dt)*2+1));
    for n=1:NR
        RecoveredT(t-t0+1,n)=trapz(Recovered((t-t0)*2+1:(t-t0+dt)*2+1,n));
        LatentRT(t-t0+1,n)=trapz(LatentR((t-t0)*2+1:(t-t0+dt)*2+1,n));
        InfectedRT(t-t0+1,n)=trapz(InfectedR((t-t0)*2+1:(t-t0+dt)*2+1,n));
        RecoveredRT(t-t0+1,n)=trapz(RecoveredR((t-t0)*2+1:(t-t0+dt)*2+1,n));
    end
    %%
    for n=1:NR
    %Numerator
    ImmunizedCasesT(t-t0+1,n)=(LatentRT(t-t0+1,n)+InfectedRT(t-t0+1,n));
    UnImmunizedControlT(t-t0+1,n)=SusceptibleT(t-t0+1);
       
    %Denominator
    ImmunizedControlT(t-t0+1,n)=(RecoveredT(t-t0+1,n)+RecoveredRT(t-t0+1,n));
    UnImmunizedCasesT(t-t0+1,n)=(LatentT(t-t0+1)+InfectedT(t-t0+1));
    end
    
    %Numerator
    TotalImmunizedCasesT(t-t0+1)=sum(LatentRT(t-t0+1,:)+InfectedRT(t-t0+1,:));
    TotalUnImmunizedControlT(t-t0+1)=SusceptibleT(t-t0+1);
       
    %Denominator
    TotalImmunizedControlT(t-t0+1)=sum(RecoveredT(t-t0+1,:)+RecoveredRT(t-t0+1,:));
    TotalUnImmunizedCasesT(t-t0+1)=(LatentT(t-t0+1)+InfectedT(t-t0+1));
   
    for n=1:NR
    ORatioNaturalInfectionT(t-t0+1,n)=(ImmunizedCasesT(t-t0+1,n)*UnImmunizedControlT(t-t0+1,n))/(ImmunizedControlT(t-t0+1,n)*UnImmunizedCasesT(t-t0+1,n)+eps);
    NaturalInfectionEffectivnessT(t-t0+1,n)=(1-ORatioNaturalInfectionT(t-t0+1,n))*100;
        
    NumeratorRealEffectivnessT(t-t0+1,n)=trapz(NumeratorRealEffectivness((t-t0)*2+1:(t-t0+dt)*2+1,n));
    DenominatorRealEffectivnessT(t-t0+1,n)=trapz(DenominatorRealEffectivness((t-t0)*2+1:(t-t0+dt)*2+1,n));
    RealNaturalInfectionEffectivnessT(t-t0+1,n)=100.*NumeratorRealEffectivnessT(t-t0+1,n)/DenominatorRealEffectivnessT(t-t0+1,n);
    end
    
     
    ORatioNaturalInfectionTAll(t-t0+1)=(TotalImmunizedCasesT(t-t0+1).*TotalUnImmunizedControlT(t-t0+1))/(TotalImmunizedControlT(t-t0+1).*TotalUnImmunizedCasesT(t-t0+1)+eps);
    NaturalInfectionEffectivnessTAll(t-t0+1)=(1-ORatioNaturalInfectionTAll(t-t0+1))*100;   
    RealNaturalInfectionEffectivnessTAll(t-t0+1)=100.*sum(NumeratorRealEffectivnessT(t-t0+1,:))/sum(DenominatorRealEffectivnessT(t-t0+1,:));

end

for n=1:NR-1
    NaturalInfectionEffectivnessTMontly(n)=NaturalInfectionEffectivnessT(600-t0+1,n);
    RealNaturalInfectionEffectivnessTMontly(n)=RealNaturalInfectionEffectivnessT(600-t0+1,n);

end
    NaturalInfectionEffectivnessTMontly(NR)=NaturalInfectionEffectivnessT(600-t0+1,NR);
    RealNaturalInfectionEffectivnessTMontly(NR)=RealNaturalInfectionEffectivnessT(600-t0+1,NR);

%%%Figure 4B
subplot(212)
hold on
bar(1.2:24.2,RealNaturalInfectionEffectivnessT(600-t0+1,:),'r')
bar(1:24,NaturalInfectionEffectivnessT(600-t0+1,:),'b')
ylabel('Estimated effectiveness of prior infection in preventing reinfection (%)')
xlabel('Months')
legend('PE_S^{true} in presence of gradual waning of protection against reinfection (using the Waning of Immunity Model)','PE_S^{test-negative} in presence of gradual waning of protection against reinfection (using the Waning of Immunity Model)')
