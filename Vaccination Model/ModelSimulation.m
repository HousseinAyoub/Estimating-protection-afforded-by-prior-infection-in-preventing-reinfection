%% Simulation the model using the estimated parameters
clc
clear all
close all

InciQatar='DataCOVIDQatarHOA.xlsx';
IncidenceQatar= xlsread(InciQatar,1);

load('EstimatedParamQatarNew.mat','EstimatedParamQatarNew') %%Estimated from Qatar Model
zz=EstimatedParamQatarNew;
a0L(1)=zz(1);
a0L(2)=zz(2);

a0L(3)=zz(3);
a0L(4)=zz(4);

a0L(5)=zz(5);
a0L(6)=zz(6);

eps=1.0e-10;


Ncomp=14;
mu=1/(80.7*365);               % 1/life expectancy

delta=1/3.69;  % Duration of stay in E is in average 3.69 days
nu=1/3.48;    % Duration of stay in I is in average 3.48 days
alpha=1.85/10000;
R0=4;

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

EffectivnessReinfection=0.8;
VEIReinfection=0;

VES=0.75;
WaningVaccinated=1/(6*30); %%6 Months based on the NEJM paper

RateofVaccination=2.6193e-05; %%Estimated from Qatar Model
VaccinationRate=zeros(length(tspan),1);
for t=t0+306:dt:length(IncidenceQatar(:,3))+t0-1+100
    VaccinationRate((t-t0)*2+1)=RateofVaccination*(t-t0-306+1);%RateofVaccination;
end
%%Deterministic model
[T,x] = ode45(@(t,x)ODE_Equations(t,x,t0,a0L,Ncomp,VEIReinfection,delta,nu,alpha,mu,EffectivnessReinfection,VaccinationRate,VES,WaningVaccinated),tspan,x0);
TT=length(T);
for i=1:TT 
     [dxL,Lambda]=ODE_Equations(T(i),x(i,:),t0,a0L,Ncomp,VEIReinfection,delta,nu,alpha,mu,EffectivnessReinfection,VaccinationRate,VES,WaningVaccinated);
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

Vaccinated=x(:,7+1);
LatentVaccinated=x(:,7+2);
InfectedVaccinated=x(:,7+3);

RecoveredVaccinated=x(:,7+4);
LatentRVaccinated=x(:,7+5);
InfectedRVaccinated=x(:,7+6);
RecoveredRVaccinated=x(:,7+7);
   
Total=sum(x,2);
BiasP=0.75;
BiasI=0.75;
Bias=0.75;
for t=1:TT
     IncidenceTotal(t)=delta.*(Latent(t)+LatentR(t)+LatentVaccinated(t)+LatentRVaccinated(t));
     DailyVaccinated(t)=VaccinationRate(t)*(Susceptible(t)+Recovered(t)+RecoveredR(t));
     CumDailyVaccinated(t)=dt*sum(DailyVaccinated(1:t));

     VaccinatedALL(t)=Vaccinated(t)+LatentVaccinated(t)+InfectedVaccinated(t)+RecoveredVaccinated(t)+...
                  LatentRVaccinated(t)+InfectedRVaccinated(t)+RecoveredRVaccinated(t); 
%         %% No Bias
        %Numerator
        ImmunizedCases(t)=(LatentR(t)+InfectedR(t)+LatentRVaccinated(t)+InfectedRVaccinated(t));
        UnImmunizedControl(t)=Susceptible(t)+Vaccinated(t);
       
        %Denominator
        ImmunizedControl(t)=(Recovered(t)+RecoveredR(t)+RecoveredVaccinated(t)+RecoveredRVaccinated(t));
        UnImmunizedCases(t)=(Latent(t)+Infected(t)+LatentVaccinated(t)+InfectedVaccinated(t));
%         
%         %%Bias only Prior infection
%         %Numerator
%         ImmunizedCases(t)=(LatentR(t)+InfectedR(t)+LatentRVaccinated(t)+InfectedRVaccinated(t))*(1-Bias);
%         UnImmunizedControl(t)=(Susceptible(t)+Vaccinated(t))+...
%                               (Recovered(t)+RecoveredR(t)+RecoveredVaccinated(t)+RecoveredRVaccinated(t))*(Bias);
%        
%         %Denominator
%         ImmunizedControl(t)=(Recovered(t)+RecoveredR(t)+RecoveredVaccinated(t)+RecoveredRVaccinated(t))*(1-Bias);
%         UnImmunizedCases(t)=(Latent(t)+Infected(t)+LatentVaccinated(t)+InfectedVaccinated(t))+...
%                             (LatentR(t)+InfectedR(t)+LatentRVaccinated(t)+InfectedRVaccinated(t))*(Bias);

%         %% Combined Biases
%         %Numerator
%         ImmunizedCases(t)=(LatentR(t)+InfectedR(t)+LatentRVaccinated(t)+InfectedRVaccinated(t))*(1-BiasP)*(1-BiasI);
%         UnImmunizedControl(t)=Susceptible(t)+Vaccinated(t)+...
%                             (Latent(t)+Infected(t)+LatentVaccinated(t)+InfectedVaccinated(t))*(BiasI)+...
%                             (Recovered(t)+RecoveredR(t)+RecoveredVaccinated(t)+RecoveredRVaccinated(t))*(BiasP)+...
%                             (LatentR(t)+InfectedR(t)+LatentRVaccinated(t)+InfectedRVaccinated(t))*(BiasP)*(BiasI);
%        
%         %Denominator
%         ImmunizedControl(t)=(Recovered(t)+RecoveredR(t)+RecoveredVaccinated(t)+RecoveredRVaccinated(t))*(1-BiasP)+...
%                             (LatentR(t)+InfectedR(t)+LatentRVaccinated(t)+InfectedRVaccinated(t))*(1-BiasP)*(BiasI);
%         
%         UnImmunizedCases(t)=(Latent(t)+Infected(t)+LatentVaccinated(t)+InfectedVaccinated(t))*(1-BiasI)+...
%                            (LatentR(t)+InfectedR(t)+LatentRVaccinated(t)+InfectedRVaccinated(t))*(BiasP)*(1-BiasI);

        

        cases(t)=ImmunizedCases(t)+UnImmunizedCases(t);
        controls(t)=UnImmunizedControl(t)+ImmunizedControl(t);
        
        ORatioNaturalInfection(t)= (ImmunizedCases(t)*UnImmunizedControl(t))/(ImmunizedControl(t)*UnImmunizedCases(t));
        NaturalInfectionEffectivness(t)=(1-ORatioNaturalInfection(t))*100;
        
        NumeratorRealEffectivness(t)=EffectivnessReinfection*(Recovered(t)+LatentR(t)+InfectedR(t)+RecoveredR(t));
        DenominatorRealEffectivness(t)=(Recovered(t)+LatentR(t)+InfectedR(t)+RecoveredR(t));
      
end
for t=t0:tf-1
    IncidenceTotalT(t-t0+1)=trapz(IncidenceTotal((t-t0)*2+1:(t-t0+dt)*2+1));
    VaccinatedALLT(t-t0+1)=trapz(VaccinatedALL((t-t0)*2+1:(t-t0+dt)*2+1));
    TotalT(t-t0+1)=trapz(Total((t-t0)*2+1:(t-t0+dt)*2+1));
    CoverageT(t-t0+1)=100.*VaccinatedALLT(t-t0+1)/TotalT(t-t0+1);
    
    CumDailyVaccinatedT(t-t0+1)=trapz(CumDailyVaccinated((t-t0)*2+1:(t-t0+dt)*2+1));
    CoverageTEver(t-t0+1)=100.*CumDailyVaccinatedT(t-t0+1)/TotalT(t-t0+1);

       ImmunizedCasesT(t-t0+1)=trapz(ImmunizedCases((t-t0)*2+1:(t-t0+dt)*2+1));
        UnImmunizedControlT(t-t0+1)=trapz(UnImmunizedControl((t-t0)*2+1:(t-t0+dt)*2+1));
        ImmunizedControlT(t-t0+1)=trapz(ImmunizedControl((t-t0)*2+1:(t-t0+dt)*2+1));
        UnImmunizedCasesT(t-t0+1)=trapz(UnImmunizedCases((t-t0)*2+1:(t-t0+dt)*2+1));
        
        NumeratorRealEffectivnessT(t-t0+1)=trapz(NumeratorRealEffectivness((t-t0)*2+1:(t-t0+dt)*2+1));
        DenominatorRealEffectivnessT(t-t0+1)=trapz(DenominatorRealEffectivness((t-t0)*2+1:(t-t0+dt)*2+1));
    
         ORatioNaturalInfectionTAll(t-t0+1)=((ImmunizedCasesT(t-t0+1).*UnImmunizedControlT(t-t0+1)))/((ImmunizedControlT(t-t0+1).*UnImmunizedCasesT(t-t0+1))+eps);
         NaturalInfectionEffectivnessTAll(t-t0+1)=(1-ORatioNaturalInfectionTAll(t-t0+1))*100;
         RealNaturalInfectionEffectivnessTAll(t-t0+1)=100.*(NumeratorRealEffectivnessT(t-t0+1))/(DenominatorRealEffectivnessT(t-t0+1));


end

figure
hold on
plot(NaturalInfectionEffectivnessTAll,'k')
plot(RealNaturalInfectionEffectivnessTAll,'r')
xlim([0 600])
ylim([0 100])
ylabel('Natural infection effectivness (%)')
xlabel('Days')
legend('Calculated overall effectivness','Real overall effectivness')



