%% Simulation the model using the estimated parameters
clc
clear all
close all

load('EstimatedParamQatarNew.mat','EstimatedParamQatarNew') %%Estimated Parameters from the Qatar Model
zz=EstimatedParamQatarNew;
a0L(1)=zz(1);
a0L(2)=zz(2);

a0L(3)=zz(3);
a0L(4)=zz(4);

a0L(5)=zz(5);
a0L(6)=zz(6);

eps=1.0e-10;


Ncomp=7;
mu=1/(80.7*365);               % 1/life expectancy

delta=1/3.69;  % Duration of stay in E is in average 3.69 days
nu=1/3.48;    % Duration of stay in I is in average 3.48 days
alpha=1.85/10000;

t0=110;
dt=0.5; %%for the stability of the scheme
tf=t0+24*30; %Final date
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

EffectivnessReinfection=0.8; %%True PES
VEIReinfection=0;

%%Deterministic model
[T,x] = ode45(@(t,x)ODE_Equations(t,x,t0,a0L,Ncomp,VEIReinfection,delta,nu,alpha,mu,EffectivnessReinfection),tspan,x0);
TT=length(T);
for i=1:TT 
     [dxL,Lambda]=ODE_Equations(T(i),x(i,:),t0,a0L,Ncomp,VEIReinfection,delta,nu,alpha,mu,EffectivnessReinfection);
     lmbL(i)=Lambda;
end

TT=length(T);
Susceptible=x(:,1);
Latent=x(:,2);
Infected=x(:,3);


Recovered=x(:,4);
LatentR=x(:,5);
InfectedR=x(:,6);
RecoveredR=x(:,7);



Total=sum(x,2);
Bias=0.75;

for t=1:TT
     IncidenceTotal(t)=delta.*(Latent(t)+LatentR(t));

%% No Bias
%       %Numerator
        ImmunizedCases(t)=(LatentR(t)+InfectedR(t));
        UnImmunizedControl(t)=Susceptible(t);
       
        %Denominator
        ImmunizedControl(t)=(Recovered(t)+RecoveredR(t));
        UnImmunizedCases(t)=(Latent(t)+Infected(t));
        
%% Misclassification of prior exposure
%         ImmunizedCases(t)=(LatentR(t)+InfectedR(t))*(1-Bias);
%         ImmunizedControl(t)=(Recovered(t)+RecoveredR(t))*(1-Bias);
% 
%         UnImmunizedCases(t)=(Latent(t)+Infected(t))+(LatentR(t)+InfectedR(t))*(Bias);
%         UnImmunizedControl(t)=Susceptible(t)+(Recovered(t)+RecoveredR(t))*(Bias);
        
%% Misclassification of those with latent infection 
%         ImmunizedCases(t)=(LatentR(t)*(1-Bias)+InfectedR(t));
%         ImmunizedControl(t)=(Recovered(t)+RecoveredR(t))+LatentR(t)*(Bias);
% 
%         UnImmunizedCases(t)=(Latent(t)*(1-Bias)+Infected(t));
%         UnImmunizedControl(t)=Susceptible(t)+Latent(t)*(Bias);

%% Undereporting of those infected 
%         ImmunizedCases(t)=(LatentR(t)+InfectedR(t))*(1-Bias);
%         ImmunizedControl(t)=(Recovered(t)+RecoveredR(t))+(LatentR(t)+InfectedR(t))*(Bias);
% 
%         UnImmunizedCases(t)=(Latent(t)+Infected(t))*(1-Bias);
%         UnImmunizedControl(t)=Susceptible(t)+(Latent(t)+Infected(t))*(Bias);
% 
        ORatioNaturalInfection(t)= (ImmunizedCases(t)*UnImmunizedControl(t))/(ImmunizedControl(t)*UnImmunizedCases(t));
        EstimatedEfficacy(t)=(1-ORatioNaturalInfection(t))*100;
        
        NumeratorRealEfficacy(t)=EffectivnessReinfection*(Recovered(t)+LatentR(t)+InfectedR(t)+RecoveredR(t));
        DenominatorRealEfficacy(t)=(Recovered(t)+LatentR(t)+InfectedR(t)+RecoveredR(t));
        
        NumeratorAttackRate(t)=Latent(t)+Infected(t)+Recovered(t)+LatentR(t)+InfectedR(t)+RecoveredR(t);
        
        ProportionofControlmissedNumerator(t)=(Recovered(t)+RecoveredR(t))*(Bias);
        ControlDenomiator(t)=(Susceptible(t)+(Recovered(t)+RecoveredR(t)));
end  

for t=t0:tf-1
    IncidenceTotalT(t-t0+1)=trapz(IncidenceTotal((t-t0)*2+1:(t-t0+dt)*2+1));
        
        ImmunizedCasesT(t-t0+1)=trapz(ImmunizedCases((t-t0)*2+1:(t-t0+dt)*2+1));
        UnImmunizedControlT(t-t0+1)=trapz(UnImmunizedControl((t-t0)*2+1:(t-t0+dt)*2+1));
        ImmunizedControlT(t-t0+1)=trapz(ImmunizedControl((t-t0)*2+1:(t-t0+dt)*2+1));
        UnImmunizedCasesT(t-t0+1)=trapz(UnImmunizedCases((t-t0)*2+1:(t-t0+dt)*2+1));

        ORatioNaturalInfectionT(t-t0+1)=(ImmunizedCasesT(t-t0+1)*UnImmunizedControlT(t-t0+1))/(ImmunizedControlT(t-t0+1)*UnImmunizedCasesT(t-t0+1)+eps);
        EstimatedEfficacyTAll(t-t0+1)=(1-ORatioNaturalInfectionT(t-t0+1))*100;
        
        NumeratorRealEffectivnessT(t-t0+1)=trapz(NumeratorRealEfficacy((t-t0)*2+1:(t-t0+dt)*2+1));
        DenominatorRealEffectivnessT(t-t0+1)=trapz(DenominatorRealEfficacy((t-t0)*2+1:(t-t0+dt)*2+1));
        RealEfficacyTAll(t-t0+1)=100.*NumeratorRealEffectivnessT(t-t0+1)/DenominatorRealEffectivnessT(t-t0+1);
        
        AttackRate(t-t0+1)=100.*trapz(NumeratorAttackRate((t-t0)*2+1:(t-t0+dt)*2+1))/trapz(Total((t-t0)*2+1:(t-t0+dt)*2+1));
        ProportionofControlmissed(t-t0+1)=100.*trapz(ProportionofControlmissedNumerator((t-t0)*2+1:(t-t0+dt)*2+1))/trapz(ControlDenomiator((t-t0)*2+1:(t-t0+dt)*2+1));
end
% 
% 
% 
% 


%%%%Figure 3A and 3B 
figure
subplot(311)
hold on
plot(IncidenceTotalT(1:600),'b')
ylabel('Number of new infections per day')
xlabel('Day')
subplot(312)
hold on
plot(AttackRate,'b')
ylabel('Proportion of the population ever infected (%)')
xlabel('Day')
xlim([0 600])
subplot(313)

%%%Figure 4
figure
subplot(211)
hold on
plot(RealEfficacyTAll,'--r')
plot(EstimatedEfficacyTAll,'k')
ylim([0 100])
xlim([0 600])
ylabel('Effectiveness of prior infection in preventing reinfection (%)')
xlabel('Days')
legend('True effectiveness of prior infection in preventing reinfection (PE_S^{true})','PE_S^{test-negative} estimated in absence of bias (using the Baseline Model)')

%%%Figure 5
figure
subplot(411)
hold on
plot(RealEfficacyTAll,'--r')
plot(EstimatedEfficacyTAll,'b')
ylim([0 100])
xlim([0 600])
xlabel('Day')
legend('True effectiveness of prior infection in preventing reinfection (PE_S^{true})','PE_S^{test-negative} estimated in absence of bias (using the Baseline Model)')
ylabel('Estimated effectiveness of prior infection in preventing reinfection (%)')

