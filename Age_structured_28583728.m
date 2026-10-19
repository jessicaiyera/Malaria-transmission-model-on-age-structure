%% Biomathematics 374
% Name: Jessica Iyera
% 
% SU #: 28583728
%% Age-structured malaria transmission model

% Purpose:
% 1. Identify differences in malaria-induced mortality in two age groups in
% children below age of 5 and older individuals aged five and above.
% 2. Identify how much each age group contributes to malaria transmission
% (mosquito infection)

% Compartments:
% children (<5): Sb, Eb, Ib, Rb
% older individuals (5+): Sa, Ea, Ia, Aa, Ra
% Mosquitoes: Sv, Ev, Iv

clear;
clc;
close all;

format long
set(0, 'defaultlinelinewidth', 2)

%% ===============================================================
% MODEL PARAMETERS
% ================================================================
% Human natural mortality
mu_h = 0.00004;

% Mosquito natural mortality
mu_v = 0.052;

% Birth and recruitment rates used in this numerical simulation
% Set equal to natural mortality to keep the background
% human and mosquito populations constant (since we're not
% focusing on population growth)
Lambda_h = 0.00004;
Lambda_v = 0.052;

% Child maturation rate
tau = 0.000183;

% Mosquito biting rate
b = 0.5;

%% Transmission

beta_b = 0.27;     % Mosquito to child

beta_a = 0.27;     % Mosquito to adult

beta_v = 0.64;     % Human to mosquito

r = 0.68;          % Relative infectiousness of asymptomatic adults, kept constant

%% Disease progression

theta = 0.637;

phi_b = 0.203;

phi_a = 0.07;

phi_v = 0.10;

%% Malaria mortality

sigma_b = 0.0005;

sigma_a = 0.0002;

%% Recovery

zeta_b = 0.0014;

zeta_s = 0.01;

zeta_a = 0.002;

%% Loss of immunity

gamma_b = 0.0027;

gamma_a = 0.0027;

%% Older individuals symptomatic/asymptomatic infection transitions

rho = 0.4;

kappa = 0.05;

%% ==============================================================
% INITIAL CONDITIONS
% ==============================================================

% Children (below 5)

Sb0 = 1950;
Eb0 = 20;
Ib0 = 20;
Rb0 = 10;

% Older individuals (5+)

Sa0 = 7850;
Ea0 = 30;
Ia0 = 30;
Aa0 = 50;
Ra0 = 40;

% Mosquito vector

Sv0 = 19500;
Ev0 = 300;
Iv0 = 200;

y0 = [Sb0; Eb0; Ib0; Rb0; Sa0; Ea0; Ia0; Aa0; Ra0; Sv0; Ev0; Iv0];

%% ==============================================================
% SIMULATION
% ==============================================================

tspan = [0 365];

[t,y] = ode45(@malariaModel,tspan,y0);

%% ==============================================================
% EXTRACT COMPARTMENTS
% ==============================================================

Sb = y(:,1);
Eb = y(:,2);
Ib = y(:,3);
Rb = y(:,4);

Sa = y(:,5);
Ea = y(:,6);
Ia = y(:,7);
Aa = y(:,8);
Ra = y(:,9);

Sv = y(:,10);
Ev = y(:,11);
Iv = y(:,12);

%% ==============================================================
% POPULATION SIZES
% ==============================================================

Nb = Sb + Eb + Ib + Rb;

Na = Sa + Ea + Ia + Aa + Ra;

Nh = Nb + Na;

%% ==============================================================
% TRANSMISSION ANALYSIS
% 
% Which age group contributes most to transmission?
% ==============================================================

childContribution = b * beta_v * (Ib./Nh);

symptomaticContribution = b * beta_v * (Ia./Nh);

asymptomaticContribution = b * beta_v * (r * Aa./Nh);

olderContribution = symptomaticContribution + asymptomaticContribution;

%% Plot

figure

plot(t , childContribution)
hold on

plot(t , symptomaticContribution)

plot(t , asymptomaticContribution)

plot(t , olderContribution)

xlabel('Time (days)')
ylabel('Contribution to mosquito infection rate')

title('Contribution to malaria transmission by age group')

legend('Children <5 (Ib)', 'Symptomatic 5+ (Ia)', 'Asymptomatic 5+ (rAa)', 'Total 5+ group', 'Location','best')

grid on;

%% Annual transmission contribution (area under the curve)

childTotal = trapz(t , childContribution);

olderTotal = trapz(t , olderContribution);

asymTotal = trapz(t , asymptomaticContribution);

overallTotal = childTotal + olderTotal;

fprintf('\n')
fprintf('TRANSMISSION RESULTS:\n')

fprintf('Children <5 contribution = %.1f%%\n', 100 * childTotal / overallTotal)

fprintf('Older group 5+ contribution = %.1f%%\n', 100 * olderTotal / overallTotal)

fprintf('Asymptomatic contribution overall = %.1f%%\n', 100 * asymTotal / overallTotal)

fprintf('Asymptomatic share within 5+ = %.1f%%\n', 100 * asymTotal / olderTotal)

%% ==============================================================
% MORTALITY ANALYSIS
% 
% Which age group experiences higher malaria mortality?
% ==============================================================

childMortality = sigma_b * Ib;

olderMortality = sigma_a * Ia;

%% Total deaths

childDeaths = trapz(t , childMortality);

olderDeaths = trapz(t , olderMortality);

%% Mortality risk

childRisk = childDeaths / Nb(1);

olderRisk = olderDeaths / Na(1);

fprintf('\n')
fprintf('MORTALITY RESULTS:\n')

fprintf('Children 5< deaths = %.2f\n', childDeaths)

fprintf('older group 5+ deaths = %.2f\n', olderDeaths)

fprintf('Children share of deaths = %.1f%%\n', 100 * childDeaths / (childDeaths + olderDeaths))

fprintf('Deaths per 100 children = %.2f\n', 100 * childRisk)

fprintf('Deaths per 100 older group = %.2f\n', 100 * olderRisk)

fprintf('Children are %.1f times more likely to die\n', childRisk / olderRisk)

%% Mortality per person

Mb = childMortality./Nb;

Ma = olderMortality./Na;

figure

plot(t , Mb)
hold on

plot(t , Ma)

xlabel('Time (days)')
ylabel('Malaria deaths per person per day')

title('Malaria-induced mortality by age group')

legend('Children <5', 'Older group 5+', 'Location','best')

grid on;

%% ==============================================================
% MODEL EQUATIONS
% ==============================================================

function dydt = malariaModel(~,y)

%% Parameters

mu_h = 0.00004;
mu_v = 0.052;

Lambda_h = 0.00004;
Lambda_v = 0.052;

tau = 0.000183;

b = 0.5;

beta_b = 0.27;
beta_a = 0.27;
beta_v = 0.64;

theta = 0.637;

phi_b = 0.203;
phi_a = 0.07;
phi_v = 0.10;

sigma_b = 0.0005;
sigma_a = 0.0002;

zeta_b = 0.0014;
zeta_s = 0.01;
zeta_a = 0.002;

gamma_b = 0.0027;
gamma_a = 0.0027;

rho = 0.4;
kappa = 0.05;

r = 0.68;

%% Compartments

Sb = y(1);
Eb = y(2);
Ib = y(3);
Rb = y(4);

Sa = y(5);
Ea = y(6);
Ia = y(7);
Aa = y(8);
Ra = y(9);

Sv = y(10);
Ev = y(11);
Iv = y(12);

%% Population size

Nh = Sb + Eb + Ib + Rb + Sa + Ea + Ia + Aa + Ra;

Nv = Sv + Ev + Iv;

%% Forces of infection

lambda_b = b*beta_b*Iv/Nh;

lambda_a = b*beta_a*Iv/Nh;

lambda_v = b*beta_v*(Ib + Ia + r*Aa)/Nh;

%% Children (<5)

dSb = Lambda_h*Nh + gamma_b*Rb - (mu_h + tau + lambda_b)*Sb;

dEb = lambda_b*Sb - (mu_h + tau + phi_b)*Eb;

dIb = phi_b*Eb - (mu_h + sigma_b + tau + zeta_b)*Ib;

dRb = zeta_b*Ib - (mu_h + tau + gamma_b)*Rb;

%% Older individuals (5+ years)

dSa = tau*Sb + gamma_a*Ra - (mu_h + lambda_a)*Sa;

dEa = lambda_a*Sa + tau*Eb - (mu_h + phi_a)*Ea;

dIa = theta*phi_a*Ea + tau*Ib - (mu_h + sigma_a + rho*zeta_s + (1-rho)*kappa)*Ia;

dAa = (1-theta)*phi_a*Ea + (1-rho)*kappa*Ia - (mu_h + zeta_a)*Aa;

dRa = zeta_a*Aa + rho*zeta_s*Ia + tau*Rb - (mu_h + gamma_a)*Ra;

%% Mosquito vector

dSv = Lambda_v*Nv - (mu_v + lambda_v)*Sv;

dEv = lambda_v*Sv - (mu_v + phi_v)*Ev;

dIv = phi_v*Ev - mu_v*Iv;

%% Return derivatives

dydt = [dSb;
        dEb;
        dIb;
        dRb;
        dSa;
        dEa;
        dIa;
        dAa;
        dRa;
        dSv;
        dEv;
        dIv];

end