%% Initial code was written by Emma Stump, which was changed and added to by Melody Walker

%%
% This is for Fig. 2 initial run with all values
IC = [30 0 30 0];
tspan = [0 2000];

% [beta_a, delta_a, mu_0, mu_2, beta_t, delta_t, mu_1, mu_3, alpha_a, alpha_t, p_t, p_a, gb, gd, gm]'
par = [32.6, 0.1, 0.067, 0.045, 11, 0.04, 0.011, 0.1, 0.83,0.25,.5, .5,1,1,1, 60];
 [t,y] = ode45(@(t,y)competition_ode(t,y,par ),tspan,IC);
 
 par = [32.6, 0.1, 0.067, 0.045, 11, 0.04, 0.011, 0.1, 0.42,0.73,.5, .5,1,1,1, 60];
 [t2,y2] = ode45(@(t,y)competition_ode(t,y,par ),tspan,IC);
 %%
 figure
subplot(1,2,1)
plot( t,y(:,1),'--', 'LineWidth', 2)
hold on
 
plot( t,y(:,2),':', 'LineWidth', 2)

plot( t,y(:,3), '-.','LineWidth', 2)

plot( t,y(:,4) ,'LineWidth', 2)
hold on
ylabel('Population Size')
xlabel('Time')
set(gca, 'FontSize', 18)
 axis([0,150,0, 140])
 
 
 subplot(1,2,2)
plot( t2,y2(:,1),'--', 'LineWidth', 2)
hold on
 
plot( t2,y2(:,2),':', 'LineWidth', 2)

plot( t2,y2(:,3), '-.','LineWidth', 2)

plot( t2,y2(:,4) ,'LineWidth', 2)
hold on
%ylabel('Population Size')
xlabel('Time')
set(gca, 'FontSize', 18)
 axis([0,150,0, 140])
 legend('L_a', 'A_a', 'L_t', 'A_t', 'Location', 'NorthEast')

%%
% here we set up the LHS, we use the same one for both, but set the
% parasitism to one the first time

s = RandStream('mt19937ar','Seed',0, 'NormalTransform', 'Ziggurat');

%Initial Conditions [Al Aa Tl Ta]
% Al => Albopictus Larvae
% Al => Albopictus Adult
% Tl => Triseriatus Larvae
% Ta => Triseriatus Adult

RandStream.setGlobalStream(s)
samp = 1e5;
Par = lhsdesign(15,samp);

% order of parameters in bounds
% [beta_a, delta_a, mu_0, mu_2, beta_t, delta_t, mu_1, mu_3, alpha_a, alpha_t, p_t, p_a, gb, gd, gm]'

lb = [2.5, 1/45, 0.005, 0.01, 3, 1/55, 0.002, 0.03, 0.4, 0,0.2, 0.4, 1, 1, 1 ]';

ub = [56, 1/9, 0.4, 0.065, 26, 1/13, 0.011, 0.1, 1, 0.75, 0.6, 0.55, 4, 4, 8]';

par200 = bsxfun(@plus, lb, bsxfun(@times, Par, (ub-lb))); 

IC = [30 0 30 0];
tspan = 0:.5:2000;
% options2 = odeset('NonNegative',1);
final_population = zeros(samp, 4);
final_populationg = zeros(samp, 4);
%% This runs all parameter sets without parasitism

tic
for j = 1:10
    par = [par200(1:12,(1 + 10000*(j-1)):10000*j); ones(3,10000); 60.*ones(1,10000)];
    k =  10000*(j-1);
for i = 1:10000
    % Span over which to solve ODE system
    
    [t,y] = ode45(@(t,y)competition_ode(t,y,par(:,i) ),tspan,IC);
    final_population(k + i,:) = y(end,:); 
end
%save(sprintf('lhsEm2_%d.mat',j))

end
toc
% options2 = odeset('NonNegative',1);

% All with parasitism
%
for j = 1:10
    par = [par200(:,(1 + 10000*(j-1)):10000*j); 60.*ones(1,10000)];
    k =  10000*(j-1);
parfor i = 1:10000
    % Span over which to solve ODE system
    
    [t,y] = ode45(@(t,y)competition_ode(t,y,par(:,i) ),tspan,IC);
    final_populationg(k + i,:) = y(end,:); 
end
save(sprintf('lhsEmg2_%d.mat',j))

end



%%

% This is the proportion of albopictus with and without parasitism

    c = final_population(:,2)./((final_population(:,2) + final_population(:,4)));
    cg = final_populationg(:,2)./((final_populationg(:,2) + final_populationg(:,4)));

    
%%

% This is for Fig.3
% separating by categories
cat_1 = sum(c(:,1) < 0.01);
cat_2 = sum(c(:,1) < 0.20 & c(:,1) > 0.01) ;
cat_3 = sum(c(:,1) < 0.40 & c(:,1) > 0.2);
cat_4 = sum(c(:,1) > 0.40 & c(:,1) < 0.6);
cat_5 = sum(c(:,1) > 0.60 & c(:,1) < 0.80);
cat_6 = sum(c(:,1) > 0.80 & c(:,1) < 0.99);
cat_7 = sum(c(:,1) > 0.99);

sum_cat = (cat_1 + cat_2 + cat_3 + cat_5 + cat_6 + cat_7 + cat_4)/samp;


freq_1 = repelem(1,cat_1);
freq_2 = repelem(2,cat_2);
freq_3 = repelem(3,cat_3);
freq_4 = repelem(4,cat_4);
freq_5 = repelem(5,cat_5);
freq_6 = repelem(6,cat_6);
freq_7 = repelem(7,cat_7);


color_1 = [0.234375,0.1484375,0.65625];
color_2 = [0.2734375,0.25390625,0.8984375];
color_3 = [0.24609375,0.43359375,0.99609375];
color_4 = [0.1015625,0.67578125,0.85546875];
color_5 = [0.55859375,0.7890625,0.30078125];
color_6 = [0.984375,0.78515625,0.19140625];
color_7 = [0.96875,0.9765625,0.07421875];
% Making the histogram
%
figure
subplot(1,2,1)
% Category 1
histogram(freq_1,'FaceColor',color_1)
hold on
% Category 2 
histogram(freq_2,'FaceColor',color_2)
hold on
% Category 3
histogram(freq_3,'FaceColor',color_3)
hold on
% Category 4
histogram(freq_4,'FaceColor',color_4)
hold on
% Category 5
histogram(freq_5,'FaceColor',color_5)
hold on
% Category 6
histogram(freq_6,'FaceColor',color_6)
hold on
% Category 7
histogram(freq_7,'FaceColor',color_7)
% legend('A_a < 0.01', '0.01 < A_a < 0.15','0.15 < A_a < 0.35','0.35 < A_a < 0.65','0.65 < A_a < 0.85', '0.85 < A_a < 0.99', '0.99 < A_a','Location', 'NorthWest')


set(gca, 'XTick',  1:7, 'XTickLabel',  1:7)


%title('Competition Outcome Frequencies')
xlabel('Outcome Category')
ylabel('Number of Occurences')

set(gca, 'FontSize', 20)
axis([0.5,7.5, 0 , 5.8e4])
set(gca,'Position', [ 0.0800    0.1250    0.413    0.8150])
%set(gcf,'Position', [ 252   176   668   522])


cat_1 = sum(cg(:,1) < 0.01);
cat_2 = sum(cg(:,1) < 0.20 & cg(:,1) > 0.01) ;
cat_3 = sum(cg(:,1) < 0.40 & cg(:,1) > 0.2);
cat_4 = sum(cg(:,1) > 0.40 & cg(:,1) < 0.6);
cat_5 = sum(cg(:,1) > 0.60 & cg(:,1) < 0.80);
cat_6 = sum(cg(:,1) > 0.80 & cg(:,1) < 0.99);
cat_7 = sum(cg(:,1) > 0.99);


sum_cat = (cat_1 + cat_2 + cat_3 + cat_5 + cat_6 + cat_7)/samp;


cat_combined = [cat_1, cat_2, cat_3, cat_4, cat_5, cat_6, cat_7];

freq_1 = repelem(1,cat_1);
freq_2 = repelem(2,cat_2);
freq_3 = repelem(3,cat_3);
freq_4 = repelem(4,cat_4);
freq_5 = repelem(5,cat_5);
freq_6 = repelem(6,cat_6);
freq_7 = repelem(7,cat_7);



%figure
subplot(1,2,2)
% Category 1
histogram(freq_1,'FaceColor',color_1)
hold on
% Category 2 
histogram(freq_2,'FaceColor',color_2)
hold on
% Category 3
histogram(freq_3,'FaceColor',color_3)
hold on
% Category 4
histogram(freq_4,'FaceColor',color_4)
hold on
% Category 5
histogram(freq_5,'FaceColor',color_5)
hold on
% Category 6
histogram(freq_6,'FaceColor',color_6)
hold on
% Category 7
histogram(freq_7,'FaceColor',color_7)
 legend('           A_a < 0.01', '0.01 < A_a < 0.20','0.20 < A_a < 0.40','0.40 < A_a < 0.60','0.60 < A_a < 0.80', '0.80 < A_a < 0.99', '0.99 < A_a','Location', 'NorthEast')

set(gca, 'YTick',  2e6, 'YTickLabel',  1)

set(gca, 'XTick',  1:7, 'XTickLabel',  1:7)
%title('Competition Outcome Frequencies')
xlabel('Outcome Category')
%ylabel('Number of Occurences')

set(gca, 'FontSize', 20)

axis([0.5,7.5, 0 ,5.8e4])
%set(gca,'Position', [ 0.0800    0.1200    0.413    0.8150])
set(gca,'Position', [ 0.5503    0.125    0.42    0.8150]) % 0.5703    0.1186    0.3347    0.8064
set(gcf,'Position', [ 196         246        959         452])

% 0.0800    0.1200    0.37    0.8150







%% this is to find the pprc Fig. 4
parb = par200';
parb(:,2) = 1./parb(:,2);
parb(:,6) = 1./parb(:,6);
rho = zeros(12,1);
p = zeros(12,1);
for i = 1:12
  [rho(i),p(i)]=partialcorr(parb(:,i),c,parb(:,[1:(i-1), (i+1):12]),'type','Spearman');
end


%pargb = parg';
rhog = zeros(15,1);
pg = zeros(15,1);
for i = 1:15
  [rhog(i),pg(i)]=partialcorr(parb(:,i),cg,parb(:,[1:(i-1), (i+1):15]),'type','Spearman');
end
% We reorder then so that parameters are in order of parasitism
r1 = rho([ 1:4, 9, 12]);
r2 = rho([ 5:8, 10, 11]);
p1 = p([ 1:4, 9, 12]);
p2 = p([ 5:8, 10, 11]);

s = {'\beta_a', '\delta_a', '\mu_{La}', '\mu_{Aa}', '\beta_t', '\delta_t', '\mu_{Lt}', '\mu_{At}', '\alpha_a', '\alpha_t', '\rho_t', '\rho_a', '\gamma_b', '\gamma_d', '\gamma_m'};

s2 = s([ 5:8, 10, 11]);
s1 = s( [1:4, 9, 12]);

rga = rhog(13:15);
rg1 = rhog([ 1:4, 9, 12]);
rg2 = rhog([ 5:8, 10, 11]);
pga = pg(13:15);
pg1 = pg([ 1:4, 9, 12]);
pg2 = pg([ 5:8, 10, 11]);


sg2 = s([ 5:8, 10, 11]);
sg1 = s( [1:4, 9, 12:15]);
%

[~,b] = sort(rg1);
subplot(2,1,1)
bar(1:6,r1(b), 'FaceColor', [ 0.8500    0.3250    0.0980])
hold on
[~,b2] = sort(rg2);
bar(12:17, r2(b2), 'FaceColor', [ 0    0.4470    0.7410])
axis([0.5, 17.5, -0.8 , 1])
hold on

%set(gca,'Position', [0.100    0.5838    0.3747    0.3412])
%set(gcf,'Position', [ 150   251   710   597])
set(gca, 'XTick',  [1:6, 12:17], 'XTickLabel',[s2(b2), s1(b)] )
set(gca, 'FontSize', 18)
%
p1a = find(p1(b) < 0.00001);
plot(p1a + .1 , r1(b(p1a))+sign(r1(b(p1a))).*.02,'*', 'Color', [ 0    0.4470    0.7410])
plot(p1a - .1 , r1(b(p1a))+sign(r1(b(p1a))).*.02,'*', 'Color', [ 0    0.4470    0.7410])


set(gca, 'FontSize', 18)
p1a = find(p2(b2) < 0.00001);
plot(p1a + .1 +11, r2(b2(p1a))+sign(r2(b2(p1a))).*.02,'*', 'Color', [ 0    0.4470    0.7410])
plot(p1a - .1 + 11, r2(b2(p1a))+sign(r2(b2(p1a))).*.02,'*', 'Color', [ 0    0.4470    0.7410])
%
set(gca,'Position', [0.090    0.5838    0.8350    0.3412])
subplot(2,1,2)

[~,b] = sort(rg1);

bar(1:6, rg1(b), 'FaceColor', [ 0.8500    0.3250    0.0980])

hold on
[~,b2] = sort(rg2);
bar(12:17, rg2(b2), 'FaceColor', [ 0    0.4470    0.7410])
axis([0.5, 17.5, -0.8 , 1])
hold on
bar(8:10, rhog([13, 15, 14]), 'FaceColor', [  0.4940    0.1840    0.5560])
%

%set(gca, 'XTick',  1:9, 'XTickLabel', sg1(b) )
set(gca, 'FontSize', 18)


%p1a = find(pg([13,15,14]) < 0.00001);
plot([8:10] + .1 , rhog([13,15,14])+sign(rhog([13,15,14])).*.02,'*', 'Color', [ 0    0.4470    0.7410])
plot([8:10] - .1, rhog([13,15,14])+sign(rhog([13,15,14])).*.02,'*', 'Color', [ 0    0.4470    0.7410])


%
p1a = find(pg1(b) < 0.00001);
plot(p1a + .1, rg1(b(p1a))+sign(rg1(b(p1a))).*.02,'*', 'Color', [ 0    0.4470    0.7410])
plot(p1a - .1, rg1(b(p1a))+sign(rg1(b(p1a))).*.02,'*', 'Color', [ 0    0.4470    0.7410])
%plot(8, rho(12)+ .01,'*', 'Color', [ 0    0.4470    0.7410])


set(gca, 'XTick', [1:6, 8:10, 12:17], 'XTickLabel', [s2(b2),s([13, 15, 14]), s1(b)] )
set(gca, 'FontSize', 18)

p1a = find(pg2(b2) < 0.00001);
plot(p1a + .1 +11, rg2(b2(p1a))+sign(rg2(b2(p1a))).*.02,'*', 'Color', [ 0    0.4470    0.7410])
plot(p1a - .1+11, rg2(b2(p1a))+sign(rg2(b2(p1a))).*.02,'*', 'Color', [ 0    0.4470    0.7410])
set(gca,'Position', [0.090    0.1100    0.8350    0.3412])
set(gcf,'Position', [ 113   227   839   571])







%% This is for parasitism Fig 5 and 6



gb =1;
gd =1;
gm =1;

tspan = [0, 2000];
IC = [30 0 30 0];


par2 = [32.6, 0.1, 0.067, 0.045, 11, 0.04, 0.011, 0.1, .83,.25,.5, .5,gb, gd, gm, 60];
 %ba, da,mu0, mu2, bt dt mu1 mu3 aa at pt pa gb gd gm K   






final_populationbm = zeros(57,4);
propgbm = zeros(19*3,19);
gd = [ 1 3 5];
tic
for j = 1
gb = 1:0.5:10;
gm = 1:0.5:10;

[gb,gm] = meshgrid(gb,gm);




par2(15) = gd(j);

for i = 1:19
    
     for k = 1:19
         par2(14) = gb(i,k); %use with the mesh
         par2(16) = gm(i,k);
         [t,y] = ode45(@(t,y)competition_ode(t,y,par2),tspan,IC);
         final_populationbm(i + 19*(j-1),:) = y(end,:); 
         propgbm(i+ 19*(j-1),k) = final_populationbm(i+ 19*(j-1),2)/(final_populationbm(i+ 19*(j-1),2) + final_populationbm(i+ 19*(j-1),4));
         %answer should be proportion of one of adults (single values)
        
     end
    
end
end
toc
%
%

IC = [30 0 30 0];
tspan = [0 2000];
%
final_populationdm = zeros(57,4);
propgdm = zeros(19*3,19);
gb = [ 1 3 5];
%
for j = 1:3
gd = 1:0.5:10;
gm = 1:0.5:10;

[gd,gm] = meshgrid(gd,gm);




par2(14) = gb(j);

for i = 1:19
    
     for k = 1:19
         par2(15) = gd(i,k); %use with the mesh
         par2(16) = gm(i,k);
         [t,y] = ode45(@(t,y)competition_ode(t,y,par2),tspan,IC);
         final_populationdm(i + 19*(j-1),:) = y(end,:); 
         propgdm(i+ 19*(j-1),k) = final_populationdm(i+ 19*(j-1),2)/(final_populationdm(i+ 19*(j-1),2) + final_populationdm(i+ 19*(j-1),4));
         %answer should be proportion of one of adults (single values)
        
     end
    
end
end




final_populationdb = zeros(57,4);
propgdb = zeros(19*3,19);
gm = [ 1 3 5];
for j = 1:3
gd = 1:0.5:10;
gb = 1:0.5:10;

[gd,gb] = meshgrid(gd,gb);




par2(16) = gm(j);

for i = 1:19
    
     for k = 1:19
         par2(15) = gd(i,k); %use with the mesh
         par2(14) = gb(i,k);
         [t,y] = ode45(@(t,y)competition_ode(t,y,par2),tspan,IC);
         final_populationdb(i + 19*(j-1),:) = y(end,:); 
         propgdb(i+ 19*(j-1),k) = final_populationdb(i+ 19*(j-1),2)/(final_populationdb(i+ 19*(j-1),2) + final_populationdb(i+ 19*(j-1),4));
         %answer should be proportion of one of adults (single values)
        
     end
    
end
end

%


gb =1;
gd =1;
gm =1;

tspan = [0, 2000];
IC = [30 0 30 0];


par2 = [32.6, 0.1, 0.067, 0.045, 11, 0.04, 0.011, 0.1, .42,.73,.5, .5,gb, gd, gm, 60];







final_populationbm2 = zeros(57,4);
propgbm2 = zeros(19*3,19);
gd = [ 1 3 5];
for j = 1:3
gb = 1:0.5:10;
gm = 1:0.5:10;

[gb,gm] = meshgrid(gb,gm);




par2(15) = gd(j);

for i = 1:19
    
     for k = 1:19
         par2(14) = gb(i,k); %use with the mesh
         par2(16) = gm(i,k);
         [t,y] = ode45(@(t,y)competition_ode(t,y,par2),tspan,IC);
         final_populationbm2(i + 19*(j-1),:) = y(end,:); 
         propgbm2(i+ 19*(j-1),k) = final_populationbm2(i+ 19*(j-1),2)/(final_populationbm2(i+ 19*(j-1),2) + final_populationbm2(i+ 19*(j-1),4));
         %answer should be proportion of one of adults (single values)
        
     end
    
end
end

%
%

final_populationdm2 = zeros(57,4);
propgdm3 = zeros(19*3,19);
gb = [4,5,6];
for j = 1:3
gd = 1:0.5:10;
gm = 1:0.5:10;

[gd,gm] = meshgrid(gd,gm);




par2(13) = gb(j);

for i = 1:19
    
     for k = 1:19
         par2(14) = gd(i,k); %use with the mesh
         par2(15) = gm(i,k);
         [t,y] = ode45(@(t,y)competition_ode(t,y,par2),tspan,IC);
         final_populationdm2(i + 19*(j-1),:) = y(end,:); 
         propgdm3(i+ 19*(j-1),k) = final_populationdm2(i+ 19*(j-1),2)/(final_populationdm2(i+ 19*(j-1),2) + final_populationdm2(i+ 19*(j-1),4));
         %answer should be proportion of one of adults (single values)
        
     end
    
end
end

%
final_populationdm2 = zeros(57,4);
propgdm2 = zeros(19*3,19);
gb = [ 1 3 5];
for j = 1:3
gd = 1:0.5:10;
gm = 1:0.5:10;

[gd,gm] = meshgrid(gd,gm);




par2(13) = gb(j);

for i = 1:19
    
     for k = 1:19
         par2(14) = gd(i,k); %use with the mesh
         par2(15) = gm(i,k);
         [t,y] = ode45(@(t,y)competition_ode(t,y,par2),tspan,IC);
         final_populationdm2(i + 19*(j-1),:) = y(end,:); 
         propgdm2(i+ 19*(j-1),k) = final_populationdm2(i+ 19*(j-1),2)/(final_populationdm2(i+ 19*(j-1),2) + final_populationdm2(i+ 19*(j-1),4));
         %answer should be proportion of one of adults (single values)
        
     end
    
end
end


%


%

final_populationdb2 = zeros(57,4);
propgdb2 = zeros(19*3,19);
gm = [ 1 3 5];
for j = 1:3
gd = 1:0.5:10;
gb = 1:0.5:10;

[gd,gb] = meshgrid(gd,gb);




par2(16) = gm(j);

for i = 1:19
    
     for k = 1:19
         par2(15) = gd(i,k); %use with the mesh
         par2(14) = gb(i,k);
         [t,y] = ode45(@(t,y)competition_ode(t,y,par2),tspan,IC);
         final_populationdb2(i + 19*(j-1),:) = y(end,:); 
         propgdb2(i+ 19*(j-1),k) = final_populationdb2(i+ 19*(j-1),2)/(final_populationdb2(i+ 19*(j-1),2) + final_populationdb2(i+ 19*(j-1),4));
         %answer should be proportion of one of adults (single values)
        
     end
    
end
end



%% This actual plots Fig5

g2v = 1:10;
g1v = 1:10;

h2 = figure;
subplot(1,3,1)
contourf(propgdm(1:19,:),[-0.1 .01 .2 .4 .6 .8 .99],'ShowText','on', 'LineWidth',2)
set(gca,'Position', [0.0800    0.1722    0.25    0.7448]); %    
get(gca, 'Position')
% xticks(2:2:20)
% xticklabels(g1v)
% yticks(2:2:20)
% yticklabels(g2v)
title('a')
ylabel('Effect on Mortality (\gamma_m)')
xlabel('Effect on Development Time (\gamma_d)')


caxis([-0.1 1]) %Keep this as is, it's what makes the seven distinct colors
set(gca, 'XTick',  1:2:20, 'XTickLabel',  1:10)

set(gca, 'YTick',  1:2:20, 'YTickLabel',  1:10)
set(gca, 'FontSize', 18)
axis([1,19,1,19])
%
subplot(1,3,2)
contourf(propgdm(20:38,:),[-0.1 .01 .2 .4 .6 .8 .99],'ShowText','on', 'LineWidth',2)

xlabel('Effect on Development Time (\gamma_d)')

%
title('b')
caxis([-0.1 1]) %Keep this as is, it's what makes the seven distinct colors
set(gca, 'XTick',  1:2:20, 'XTickLabel',  1:10)

%set(gca,'Position', [ 0.39    0.6    0.25    0.35]) % 0.6916   0.6    0.25   0.35
set(gca,'Position', [0.3908    0.1722    0.25    0.7448]); %   0.1722    0.25   0.0800    0.6    0.25    0.35])
%get(gca, 'Position')
set(gca, 'YTick',  25, 'YTickLabel',  1)
set(gca, 'FontSize', 18)
axis([1,19,1,19])
%
subplot(1,3,3)
contourf(propgdm(39:end,:),[-0.1 .01 .2 .4 .6 .8 .99],'ShowText','on', 'LineWidth',2)
set(gca,'Position', [ 0.7    0.1722    0.25    0.7448])%  0.7   0.6    0.25   0.35])
xlabel('Effect on Development Time (\gamma_d)')

title('c')
caxis([-0.1 1]) %Keep this as is, it's what makes the seven distinct colors
set(gca, 'XTick',  1:2:20, 'XTickLabel',  1:10)

set(gca, 'YTick',  25, 'YTickLabel',  1)
set(gca, 'FontSize', 18)
set(h2,'Position',[ 79         451        1201         347])
axis([1,19,1,19])
%

g2v = 1:10;
g1v = 1:10;

h2 = figure;
subplot(1,3,1)
contourf(propgbm(1:19,:),[-0.1 .01 .2 .4 .6 .8 .99],'ShowText','on', 'LineWidth',2)
set(gca,'Position', [0.08    0.1722    0.25    0.7448]); %   0.1722    0.25   0.0800    0.6    0.25    0.35])

xlabel('Effect on Fecundity (\gamma_b)')
ylabel('Effect on Mortality (\gamma_m)')
title('e')
axis([1,19,1,19])

caxis([-0.1 1]) %Keep this as is, it's what makes the seven distinct colors
set(gca, 'XTick',  1:2:20, 'XTickLabel',  1:10)
set(gca, 'YTick',  1:2:20, 'YTickLabel',  1:10)
set(gca, 'FontSize', 18)

subplot(2,3,5)
contourf(propgbm(20:38,:),[-0.1 .01 .2 .4 .6 .8 .99],'ShowText','on', 'LineWidth',2)
title('f')
xlabel('Effect on Fecundity (\gamma_b)')
set(gca,'Position', [0.3908    0.1722    0.25    0.7448]); %   0.1722    0.25   0.0800    0.6    0.25    0.35])

caxis([-0.1 1]) %Keep this as is, it's what makes the seven distinct colors
set(gca, 'XTick',  1:2:20, 'XTickLabel',  1:10)
set(gca, 'YTick',  25, 'YTickLabel',  5)
set(gca, 'FontSize', 18)
axis([1,19,1,19])

subplot(2,3,6)
contourf(propgbm(39:end,:),[-0.1 .01 .2 .4 .6 .8 .99],'ShowText','on', 'LineWidth',2)
set(gca,'Position', [0.7    0.1722    0.25    0.7448]); %   0.1722    0.25   0.0800    0.6    0.25    0.35])
xlabel('Effect on Fecundity (\gamma_b)')
title('g')

caxis([-0.1 1]) %Keep this as is, it's what makes the seven distinct colors
set(gca, 'XTick',  1:2:20, 'XTickLabel',  1:10)
set(gca, 'YTick',  25, 'YTickLabel',  5)
set(gca, 'FontSize', 18)
axis([1,19,1,19])
set(h2,'Position',[79         451        1201         347])
%% This plots Fig6.


g2v = 1:10;
g1v = 1:10;

h2 = figure;
subplot(1,3,1)
contourf(propgdm2(1:19,:),[-0.1 .01 .2 .4 .6 .8 .99],'ShowText','on', 'LineWidth',2)

set(gca,'Position', [0.08    0.1722    0.25    0.7448]);

title('a')
xlabel('Effect on Development Time (\gamma_d)')
ylabel('Effect on Mortality (\gamma_m)')

axis([1,19,1,19])


caxis([-0.1 1]) %Keep this as is, it's what makes the seven distinct colors
set(gca, 'XTick',  1:2:20, 'XTickLabel',  1:10)
set(gca, 'YTick',  1:2:20, 'YTickLabel',  1:10)

set(gca, 'FontSize', 18)

subplot(1,3,2)
contourf(propgdm2(20:38,:),[-0.1 .01 .2 .4 .6 .8 .99],'ShowText','on', 'LineWidth',2)

xlabel('Effect on Development Time (\gamma_d)')

caxis([-0.1 1]) %Keep this as is, it's what makes the seven distinct colors

title('b')
set(gca,'Position', [0.39    0.1722    0.25    0.7448]);
axis([1,19,1,19])
set(gca, 'YTick',  1:2:20, 'YTickLabel',  1:10)
set(gca, 'FontSize', 18)
set(gca, 'XTick',  1:2:20, 'XTickLabel',  1:10)
set(gca, 'YTick',  25, 'YTickLabel',  1)

subplot(1,3,3)
contourf(propgdm2(39:end,:),[-0.1 .01 .2, .4 .6 .8 .99],'ShowText','on', 'LineWidth',2)
%set(gca,'Position', [ .7, 0.15, 0.27, 0.8])
%xlabel('Parasite Effect on Fecundity (\gamma_b)')
xlabel('Effect on Development Time (\gamma_d)')
title('c')
axis([1,19,1,19])

caxis([-0.1 1]) %Keep this as is, it's what makes the seven distinct colors
set(gca, 'XTick',  1:2:20, 'XTickLabel',  1:10)
set(gca, 'YTick',  25, 'YTickLabel',  1)
set(gca, 'FontSize', 18)
set(gca,'Position', [0.7    0.1722    0.25    0.7448]);
set(h2,'Position',[79         451        1201         347])
%



g2v = 1:10;
g1v = 1:10;

h2 = figure;
subplot(1,3,1)
contourf(propgbm2(1:19,:),[-0.1 .01 .2 .4 .6 .8 .99],'ShowText','on', 'LineWidth',2)
set(gca,'Position', [0.08    0.1722    0.25    0.7448]);
% xticks(2:2:20)
% xticklabels(g1v)
% yticks(2:2:20)
% yticklabels(g2v)
xlabel('Effect on Fecundity (\gamma_b)')
ylabel('Effect on Mortality (\gamma_m)')
title('e')
axis([1,19,1,19])

caxis([-0.1 1]) %Keep this as is, it's what makes the seven distinct colors
set(gca, 'XTick',  1:2:20, 'XTickLabel',  1:10)
set(gca, 'YTick',  1:2:20, 'YTickLabel',  1:10)
set(gca, 'FontSize', 18)

subplot(1,3,2)
contourf(propgbm2(20:38,:),[-0.1 .01 .2 .4 .6 .8 .99],'ShowText','on', 'LineWidth',2)
title('f')
xlabel('Effect on Fecundity (\gamma_b)')
set(gca,'Position', [0.39    0.1722    0.25    0.7448]);

caxis([-0.1 1]) %Keep this as is, it's what makes the seven distinct colors
set(gca, 'XTick',  1:2:20, 'XTickLabel',  1:10)
set(gca, 'YTick',  25, 'YTickLabel',  5)
set(gca, 'FontSize', 18)
axis([1,19,1,19])

subplot(1,3,3)
contourf(propgbm2(39:end,:),[-0.1 .01 .2 .4 .6 .8 .99],'ShowText','on', 'LineWidth',2)
set(gca,'Position', [0.7    0.1722    0.25    0.7448]);
%xlabel('Parasite Effect on Fecundity (\gamma_b)')
xlabel('Effect on Fecundity (\gamma_b)')
title('g')

caxis([-0.1 1]) %Keep this as is, it's what makes the seven distinct colors
set(gca, 'XTick',  1:2:20, 'XTickLabel',  1:10)
set(gca, 'YTick',  25, 'YTickLabel',  5)
set(gca, 'FontSize', 18)
axis([1,19,1,19])

set(h2,'Position',[79         451        1201         347])




%% This is the other one left out


g2v = 1:10;
g1v = 1:10;

h2 = figure;
subplot(2,3,1)
contourf(propgdb(1:19,:),[-0.1 .01 .2 .4 .6 .8 .99],'ShowText','on', 'LineWidth',2)
set(gca,'Position', [ 0.0800    0.6    0.25    0.35])
% xticks(2:2:20)
% xticklabels(g1v)
% yticks(2:2:20)
% yticklabels(g2v)
title('a')
ylabel('Effect on Fecundity (\gamma_b)')
%xlabel('Effect on Development Time (\gamma_d)')


caxis([-0.1 1]) %Keep this as is, it's what makes the seven distinct colors
set(gca, 'XTick',  1:2:20, 'XTickLabel',  1:10)

set(gca, 'YTick',  1:2:20, 'YTickLabel',  1:10)
set(gca, 'FontSize', 18)
axis([1,19,1,19])
%
subplot(3,3,2)
contourf(propgdb(20:38,:),[-0.1 .01 .2 .4 .6 .8 .99],'ShowText','on', 'LineWidth',2)

%xlabel('Effect on Development Time (\gamma_d)')


title('b')
caxis([-0.1 1]) %Keep this as is, it's what makes the seven distinct colors
set(gca, 'XTick',  1:2:20, 'XTickLabel',  1:10)

set(gca,'Position', [ 0.39    0.6    0.25    0.35]) % 0.6916   0.6    0.25   0.35
set(gca, 'YTick',  25, 'YTickLabel',  1)
set(gca, 'FontSize', 18)
axis([1,19,1,19])
%
subplot(2,3,3)
contourf(propgdb(39:end,:),[-0.1 .01 .2 .4 .6 .8 .99],'ShowText','on', 'LineWidth',2)
set(gca,'Position', [0.7   0.6    0.25   0.35])
%xlabel('Effect on Development Time (\gamma_d)')

title('c')
caxis([-0.1 1]) %Keep this as is, it's what makes the seven distinct colors
set(gca, 'XTick',  1:2:20, 'XTickLabel',  1:10)

set(gca, 'YTick',  25, 'YTickLabel',  1)
set(gca, 'FontSize', 18)
%set(h2,'Position',[49      197        1095         588])
axis([1,19,1,19])
%
% %%


g2v = 1:10;
g1v = 1:10;

%h2 = figure;
subplot(2,3,4)
contourf(propgdb2(1:19,:),[-0.1 .01 .2 .4 .6 .8 .99],'ShowText','on', 'LineWidth',2)
set(gca,'Position', [ 0.0800    0.1    0.25    0.35])
% xticks(2:2:20)
% xticklabels(g1v)
% yticks(2:2:20)
% yticklabels(g2v)
title('d')
ylabel('Effect on Fecundity (\gamma_b)')
xlabel('Effect on Development Time (\gamma_d)')


caxis([-0.1 1]) %Keep this as is, it's what makes the seven distinct colors
set(gca, 'XTick',  1:2:20, 'XTickLabel',  1:10)

set(gca, 'YTick',  1:2:20, 'YTickLabel',  1:10)
set(gca, 'FontSize', 18)
axis([1,19,1,19])
%
subplot(2,3,5)
contourf(propgdb2(20:38,:),[-0.1 .01 .2 .4 .6 .8 .99],'ShowText','on', 'LineWidth',2)

xlabel('Effect on Development Time (\gamma_d)')


title('e')
caxis([-0.1 1]) %Keep this as is, it's what makes the seven distinct colors
set(gca, 'XTick',  1:2:20, 'XTickLabel',  1:10)

set(gca,'Position', [ 0.39    0.1    0.25    0.35]) % 0.6916   0.6    0.25   0.35
set(gca, 'YTick',  25, 'YTickLabel',  1)
set(gca, 'FontSize', 18)
axis([1,19,1,19])
%
subplot(2,3,6)
contourf(propgdb2(39:end,:),[-0.1 .01 .2 .4 .6 .8 .99],'ShowText','on', 'LineWidth',2)
set(gca,'Position', [0.7   0.1    0.25   0.35])
xlabel('Effect on Development Time (\gamma_d)')

title('f')
caxis([-0.1 1]) %Keep this as is, it's what makes the seven distinct colors
set(gca, 'XTick',  1:2:20, 'XTickLabel',  1:10)

set(gca, 'YTick',  25, 'YTickLabel',  1)
set(gca, 'FontSize', 18)
%set(h2,'Position',[49      197        1095         588])
axis([1,19,1,19])
%
% %%

set(h2,'Position',[49      117        1095         588])



%% Using the proportion value found in the text, we solved numerically for gm using mathematica
% and the values in table 1. If you compare these to the contour plots they
% agree.




% 0.8
gm = -((143.491 + 2.09655.*gb - 20.0694.*gd - 1.49254.*gb.*gd)./(1.40469.*gb.*gd - gb.*gd.^2));

figure
c = zeros(41,46,3);
c(:,:,1) =  0.9871; %0.4660  ;
     c(:,:,2) =  .7348; % 0.6740 ;
     c(:,:,3) =   .2438; %  0.1880;
surf(gb,gd,gm, c)

hold on
%axis([1,10,1, 10, 1, 10])

%
%0.6
gm = -((302.873 + 1.65948.*gb - 15.8855.*gd - 1.49254.*gb.*gd)./(1.11185.*gb.*gd -  gb.*gd.^2));


c = zeros(41,46,3);
c(:,:,1) =  0.5044; %0.4660  ;
     c(:,:,2) =  .7993; % 0.6740 ;
     c(:,:,3) =   .348; %  0.1880;

surf(gb,gd,gm, c)
hold on

%
% 0.4
gm = -((480.939 + 1.17117.*gb - 11.2111.*gd - 1.49254.*gb.*gd)./(0.784683.*gb.*gd -  gb.*gd.^2));


c = zeros(41,46,3);
c(:,:,1) = 0.0704 ;
     c(:,:,2) =  0.7457 ;
     c(:,:,3) =  0.7258;
surf(gb,gd,gm, c)
hold on
%
% 0.2
gm = -((681.18 + 0.622046.*gb - 5.95456.*gd - 1.49254.*gb.*gd)./(0.416771.*gb.*gd - gb.*gd.^2));

c = zeros(41,46,3);
c(:,:,1) = 0.154 ;
     c(:,:,2) =  0.5902;
     c(:,:,3) =  0.9218;
surf(gb,gd,gm, c)
hold on
%

% 0.01
gm = -((895.958 + 0.0330577.*gb - 0.316447.*gd - 1.49254.*gb.*gd)./(0.0221487.*gb.*gd -  gb.*gd.^2));

c = zeros(41,46,3);
c(:,:,1) = .278 ;
     c(:,:,2) =   0.3556;
     c(:,:,3) =   0.9777;
surf(gb,gd,gm, c)

hold on
axis([1,10,1, 10, 1, 10])
ylabel('\gamma_d')
xlabel('\gamma_b')
zlabel('\gamma_m')
set(gca, 'FontSize', 20)


legend('0.8', '0.6','0.4','0.2','0.01')


%%
figure
gb = 1:.2:10;
gd = 1.4:.2:10;
[gb, gd] =meshgrid(gb,gd);


%0.6
gm = -((125.258 + 1.80782.*gb - 65.6321.*gd - 1.49254.*gb.*gd)./(1.21124.*gb.*gd - 1.*gb.*gd.^2));


c = zeros(44,46,3);
c(:,:,1) =  0.5044; %0.4660  ;
     c(:,:,2) =  .7993; % 0.6740 ;
     c(:,:,3) =   .348; %  0.1880;

surf(gb,gd,gm, c)
hold on

%0.4
gm = -((193.802 + 1.24316.*gb - 45.1323.*gd - 1.49254.*gb.*gd)./(0.832917.*gb.*gd - gb.*gd.^2));
c = zeros(44,46,3);
c(:,:,1) = 0.0704 ;
     c(:,:,2) =  0.7457 ;
     c(:,:,3) =  0.7258;
surf(gb,gd,gm, c)
%0.2
gm = -((266.803 + 0.641786.*gb - 23.2997.*gd - 1.49254.*gb.*gd)./( 0.429997.*gb.*gd - gb.*gd.^2));

%c = zeros(41,46,3);
c(:,:,1) = 0.154 ;
     c(:,:,2) =  0.5902;
     c(:,:,3) =  0.9218;
surf(gb,gd,gm, c)

%0.01
gm = -((340.69 + 0.0331119.*gb - 1.20211.*gd - 1.49254.*gb.*gd)./( 0.022185.*gb.*gd -gb.*gd.^2));

%c = zeros(46,46,3);
c(:,:,1) = .278 ;
     c(:,:,2) =   0.3556;
     c(:,:,3) =   0.9777;
surf(gb,gd,gm, c)

hold on
axis([1,10,1, 10, 1, 10])
ylabel('\gamma_d')
xlabel('\gamma_b')
zlabel('\gamma_m')
set(gca, 'FontSize', 20)

set(gcf, 'Position', [265   258   735   540])
set(gca, 'CameraPosition', [-39.3928 47.2432 53.6367])
legend( '0.6','0.4','0.2','0.01')