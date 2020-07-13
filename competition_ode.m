function dy = competition_ode(t,y,par)
 %ba, da,mu0, mu2, bt dt mu1 mu3 aa at pt pa gb gd gm K           
%par = [32.6, 0.1, 0.067, 0.045, 11, 0.04, 0.011, 0.1, 0.2,1,.5, .5,gb, gd, gm, 60];

% Mosquito-Specific Parameters
% Albopictus Parameters
beta_a = par(1); % gross reproductive rate (eggs/female); % length of gonotrophic cycle (days)
Delta_a = par(2); % Ae. albopictus development time (days)
alpha_a = par(9); % intrinsic effect of albopictus on triseriatus
mu_1 = par(3); % Larval albopictus mortality
mu_2 = par(4); % Adult albopictus mortality

% Triseriatus Parameters
beta_t = par(5); 
Delta_t = par(6); % Ae. triseriatus development time
alpha_t = par(10); %1.0; % intrinsic effect of triseriatus on albopictus

mu_3 = par(7); % Larval triseriatus mortality
mu_4 = par(8); % Adult triseriatus mortality

%% general parameters
K = par(16); %60; % carrying capacity
pa = par(12); %0.50; % proportion of population that's female
pt = par(11);
%% Parasite Parameters

 gamma_1 = par(13); % parasite effect on fecundity
 gamma_2 = par(14); % parasite effect on development time
 gamma_3 = par(15); % parasite effect on larval mortality



%% Variables
Al = y(1);
Aa = y(2);
Tl = y(3);
Ta = y(4);

%% Equations
dAl = beta_a*Aa*pa*((K-Al-alpha_t*Tl)/K)./gamma_1 - mu_1*gamma_3*Al - Delta_a*Al./gamma_2; %larval albopictus
dAa = Delta_a*Al./gamma_2 - mu_2*Aa; %adult albopictus
dTl = beta_t*Ta*pt*((K - Tl - alpha_a*Al)/K) - mu_3*Tl - Delta_t*Tl; %larval triseriatus
dTa = Delta_t*Tl - mu_4*Ta; %adult triseriatus



%return vectors
dy = [dAl; dAa; dTl; dTa];
end


