%ATP utilization is insensitive to energy state in this simple 4 tank simulator: Fuel, Mito (includes Redox+DY), ATP, PCr 
clear; close all;
global iFuel iMito iATP iPCR

Indices();

% Load Initial Conditions
pathname = 'C:\Users\Ben\Desktop\ANT Model Files\Summer 2017 ANT Model Files\Hydraulic Model\';
filename = 'HydraulicLoadingSimple.xlsx';
[d] = xlsread([filename],'C3:C6');

xo = d;

% Load Parameters
[p] = xlsread([filename],'H3:H11');
param = p;

indexC = [1 3:6];
options = odeset('MaxStep',5e-2);

% Below are the intial water levels
xo(iFuel) = 8;
xo(iMito) = 8;
xo(iATP) = 8;
xo(iPCR) = 8;


% Below are the conductances of each pipe
CDH = param(1);  % this is the conductance of the fuel oxidation pathway
param(1) = 1*CDH;

COxPhos = param(2);  % this is the conductance of the ox phos pathway (includes ETC + CV & ANT)
param(2) = 1*COxPhos;

param(3)= 0; %the energy sensitive ATP breakdown rxn is zero in this script

vCONST = param(9); %this is the velocity of the energy insensitive ATP breakdown reaction

CCK = param(4);  % this is the conductance of the CKase pipe
param(4) = 1*CCK;

Xpcr = param(8); %cross-sectional area of PCr tank. Each x-fold adjusts the PCr tank by that factor
param(8) = 1*Xpcr;
Xatp = param(9); %cross-sectional area of ATP tank. Each x-fold adjusts the ATP tank by that factor
param(9) = 1*Xatp;

param(9) = 0.01*vCONST; %resting metabolic rate

tic
t_ib1 = 300; %rest
[t,x] = ode15s(@Metabolites,[0 t_ib1],xo,options,param);
toc

J1 = Fluxes(x(end,:),param);
xo = x(end,:);
x1 = x;
t1 = t;


param(9) = 0.33*vCONST;

% param(1) = 1*CDH;
% param(2) = 1*COxPhos;
% param(4) = 1*CCK;

tic
t_ib2 = 300; % exercise time
[t,x] = ode15s(@Metabolites,[0 t_ib2],xo,options,param);
toc

J2 = Fluxes(x(end,:),param);
xo = x(end,:);
x2 = x;
t2 = t;

param(9) = 0.01*vCONST;

% param(1) = 1*CDH;
% param(2) = 1*COxPhos;
% param(4) = 1*CCK;

tic
t_ib3 = 300;
[t,x] = ode15s(@Metabolites,[0 t_ib3],xo,options,param);
toc

J3 = Fluxes(x(end,:),param);
xo = x(end,:);
x3 = x;
t3 = t;

% param(1) = 1*CDH;
% param(2) = 1*COxPhos;
% param(4) = 1*CCK;

param(8) = 1*Xpcr;

CCK = param(4);  % this is the conductance of CK equilibration of ATP with PCr
param(4) = 0*CCK;

param(9) = 0.33*vCONST;

tic
t_ib4 = 300;
[t,x] = ode15s(@Metabolites,[0 t_ib4],xo,options,param);
toc

J4 = Fluxes(x(end,:),param);
xo = x(end,:);
x4 = x;
t4 = t;

param(9) = .01*vCONST;

tic
t_ib5 = 300;
[t,x] = ode15s(@Metabolites,[0 t_ib5],xo,options,param);
toc

J5 = Fluxes(x(end,:),param);
xo = x(end,:);
x5 = x;
t5 = t;

x = [x1' x2' x3' x4' x5']';
t = [t1' (t2+t_ib1)' (t3+t_ib2+t_ib1)' (t4+t_ib3+t_ib2+t_ib1)' (t5+t_ib4+t_ib3+t_ib2+t_ib1)'];

%Defining array for download into Excel

tt1=t';
y=tt1(1:30:end,1);
xlswrite('kindataNEW.xlsx',y,1,'A3:A1004');

variables=x(1:30:end,[4 3]);  %variables=x(1:200:end,[25 35 36 62 92]); %this is the full list of variables if wanted and expand Column letters accordingly e.g., 'B3:C410' for one additional variable
xlswrite('kindataNEW.xlsx',variables,1,'B3:C1004');

%Plotting
Xpcr = param(8);
Xatp = param(7);

figure('NumberTitle','off','Name','CKase & OxPhos flux')
xx=cat(1,0,-diff(x(:,4)));
tt=cat(1,0,diff(t'));
subplot(2,2,1),plot(t,x(:,4))
ylabel('[PCr]')
set(gca,'xlim',[0 (t_ib5+t_ib4+t_ib3+t_ib2+t_ib1)],'ylim',[5  8]);
subplot(2,2,2),plot(t,xx./tt*1e6*Xpcr)
ylabel('d[PCr]/dt (ml H_2O/sec)')
set(gca,'xlim',[0 (t_ib5+t_ib4+t_ib3+t_ib2+t_ib1)],'ylim',[-800  800]);
xx=cat(1,0,-diff(x(:,3)));
tt=cat(1,0,diff(t'));
subplot(2,2,3),plot(t,x(:,3))
ylabel('[ATP]')
set(gca,'xlim',[0 (t_ib5+t_ib4+t_ib3+t_ib2+t_ib1)],'ylim',[5  8]);
subplot(2,2,4),plot(t,xx./tt*1e6*Xatp)
ylabel('d[ATP]/dt (ml H_2O/sec)')
set(gca,'xlim',[0 (t_ib5+t_ib4+t_ib3+t_ib2+t_ib1)],'ylim',[-500  500]);

figure %Water levels in each tank
subplot(1,4,1),plot(t,x(:,1));
ylabel('Fuel')
set(gca,'xlim',[0 (t_ib5+t_ib4+t_ib3+t_ib2+t_ib1)],'ylim',[5 8.5]);
subplot(1,4,2),plot(t,x(:,2));
ylabel('Mito')
set(gca,'xlim',[0 (t_ib5+t_ib4+t_ib3+t_ib2+t_ib1)],'ylim',[5 8.5]);
subplot(1,4,3),plot(t,x(:,3));
ylabel('ATP')
set(gca,'xlim',[0 (t_ib5+t_ib4+t_ib3+t_ib2+t_ib1)],'ylim',[5 8.5]);
subplot(1,4,4),plot(t,x(:,4));
ylabel('PCr')
set(gca,'xlim',[0 (t_ib5+t_ib4+t_ib3+t_ib2+t_ib1)],'ylim',[5 8.5]);


figure %Fluxes
subplot(1,4,1),plot([t_ib1 (t_ib2+t_ib1) (t_ib3+t_ib2+t_ib1) (t_ib4+t_ib3+t_ib2+t_ib1) (t_ib5+t_ib4+t_ib3+t_ib2+t_ib1)],[J1(1) J2(1) J3(1) J4(1) J5(1)],'ko');
ylabel('DH (mL H2O/sec)'),xlabel('Time (Seconds)')
subplot(1,4,2),plot([t_ib1 (t_ib2+t_ib1) (t_ib3+t_ib2+t_ib1) (t_ib4+t_ib3+t_ib2+t_ib1) (t_ib5+t_ib4+t_ib3+t_ib2+t_ib1)],[J1(2) J2(2) J3(2) J4(2) J5(2)],'ko');
ylabel('OxPhos (mL H2O/sec)'),xlabel('Time (Seconds)')

subplot(1,4,3),plot([t_ib1 (t_ib2+t_ib1) (t_ib3+t_ib2+t_ib1) (t_ib4+t_ib3+t_ib2+t_ib1) (t_ib5+t_ib4+t_ib3+t_ib2+t_ib1)],[J1(5) J2(5) J3(5) J4(5) J5(5)],'ko');
ylabel('ATP Util (mL H2O/sec)'),xlabel('Time (Seconds)')
subplot(1,4,4),plot([t_ib1 (t_ib2+t_ib1) (t_ib3+t_ib2+t_ib1) (t_ib4+t_ib3+t_ib2+t_ib1) (t_ib5+t_ib4+t_ib3+t_ib2+t_ib1)],[J1(4) J2(4) J3(4) J4(4) J5(4)],'ko');
ylabel('CKase(mL H2O/sec)'),xlabel('Time (Seconds)')

% subplot(2,4,8),plot([t_ib1 (t_ib2+t_ib1) (t_ib3+t_ib2+t_ib1) (t_ib4+t_ib3+t_ib2+t_ib1)],[J1(4)/1000*x(t_ib1, iGh-iGA)*9806 J2(4)/1000*x((t_ib2+t_ib1), iGh-iGA)*9806 J3(4)/1000*x((t_ib3 + t_ib2+t_ib1), iGh-iGA)*9806 J4(4)/1000*x((t_ib4 + t_ib3 + t_ib2+t_ib1), iGh-iGA)*9806],'ko')
% ylabel('Po (watts)'),xlabel('Time (Seconds)');
% subplot(2,4,8),plot([t_ib1 (t_ib2+t_ib1) (t_ib3+t_ib2+t_ib1) (t_ib4+t_ib3+t_ib2+t_ib1)],[J1(Po) J2(Po) J3(Po) J4(Po)],'ko')
% ylabel('Po (watts)'),xlabel('Time (Seconds)');
