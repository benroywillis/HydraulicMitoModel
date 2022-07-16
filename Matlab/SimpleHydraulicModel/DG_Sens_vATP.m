%The vATP rxn is Delta G sensitive in th simple 4 tank hydraulic model simulator: Fuel, Mito (includes Redox+DY), ATP, PCr 
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
vATP = param(3); %this is the energy sensitive ATPase, which is active in this script
param(3) = 0.01*vATP; %v = 0.01 would simulate muscle resting metabolic rate

vCONST = param(9); %this energy insensitive ATP breakdown reaction is always zero in this script
param(9) = 0*vCONST; %

CDH = param(1);  % this is the conductance of the substrate (generic) ox step
param(1) = 1*CDH;

COxPhos = param(2);  % this is the conductance of the ox phos pathway (includes both ETC + CV-ANT)
param(2) = 1*COxPhos;

CCK = param(4);  % conductance of CKase
param(4) = 1*CCK;

Xpcr = param(8); %cross-sectional area of PCr tank. Each x-fold adjusts [TCr] by that factor
param(8) = 1*Xpcr;

param(3) = 0.01*vATP; %Resting metabolic rate

tic
t_ib1 = 300; % incubation period of 100 seconds
[t,x] = ode15s(@Metabolites,[0 t_ib1],xo,options,param);
toc

J1 = Fluxes(x(end,:),param);
xo = x(end,:);
x1 = x;
t1 = t;

param(3) = 0.25*vATP; % onset of first level of ATP breakdown

tic
t_ib2 = 300; % exercise time
[t,x] = ode15s(@Metabolites,[0 t_ib2],xo,options,param);
toc

J2 = Fluxes(x(end,:),param);
xo = x(end,:);
x2 = x;
t2 = t;

param(3) = 0.5*vATP; % onset of 2nd level of ATP breakdown

% xo(iPCR) = 8;
% param(1) = 4.17*Cdh;
% param(2) = 1*Cetc;
% param(3) = 2.01*Cva;


tic
t_ib3 = 300; % exercise time
[t,x] = ode15s(@Metabolites,[0 t_ib3],xo,options,param);
toc

J3 = Fluxes(x(end,:),param);
xo = x(end,:);
x3 = x;
t3 = t;

COxPhos = param(2);  % this is the conductance of the ox phos pathway
param(2) = 1*COxPhos;

% XCap = param(8); %cross-sectional area of PCr tank. Each x-fold adjusts [TCr] by that factor
% param(8) = 2*XCap;

% CCK = param(4);  % this is the conductance of CK equilibration of ATP with PCr
% param(4) = 0*CCK;

param(3) = 0.75*vATP; % onset of third level of ATP breakdown

tic
t_ib4 = 300; % exercise time
[t,x] = ode15s(@Metabolites,[0 t_ib4],xo,options,param);
toc

J4 = Fluxes(x(end,:),param);
xo = x(end,:);
x4 = x;
t4 = t;

param(3) = 1*vATP; %

tic
t_ib5 = 300; % exercise time
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

variables=x(1:30:end,4);  %variables=x(1:200:end,[25 35 36 62 92]); %this is the full list of variables if wanted and expand Column letters accordingly e.g., 'B3:C410' for one additional variable
xlswrite('kindataNEW.xlsx',variables,1,'B3:B1004');

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
ylabel('J_C_K (ml H_2O/sec)')
set(gca,'xlim',[0 (t_ib5+t_ib4+t_ib3+t_ib2+t_ib1)],'ylim',[-1000  1000]);
xx=cat(1,0,-diff(x(:,3)));
tt=cat(1,0,diff(t'));
subplot(2,2,3),plot(t,x(:,3))
ylabel('[ATP]')
set(gca,'xlim',[0 (t_ib5+t_ib4+t_ib3+t_ib2+t_ib1)],'ylim',[5  8]);
subplot(2,2,4),plot(t,xx./tt*1e6*Xatp)
ylabel('Net ATP change (ml H_2O/sec)')
set(gca,'xlim',[0 (t_ib5+t_ib4+t_ib3+t_ib2+t_ib1)],'ylim',[-1000  1000]);

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

subplot(1,4,3),plot([t_ib1 (t_ib2+t_ib1) (t_ib3+t_ib2+t_ib1) (t_ib4+t_ib3+t_ib2+t_ib1) (t_ib5+t_ib4+t_ib3+t_ib2+t_ib1)],[J1(3) J2(3) J3(3) J4(3) J5(3)],'ko');
ylabel('ATPase (mL H2O/sec)'),xlabel('Time (Seconds)')
subplot(1,4,4),plot([t_ib1 (t_ib2+t_ib1) (t_ib3+t_ib2+t_ib1) (t_ib4+t_ib3+t_ib2+t_ib1) (t_ib5+t_ib4+t_ib3+t_ib2+t_ib1)],[J1(4) J2(4) J3(4) J4(4) J5(4)],'ko');
ylabel('CKase(mL H2O/sec)'),xlabel('Time (Seconds)')

% subplot(2,4,8),plot([t_ib1 (t_ib2+t_ib1) (t_ib3+t_ib2+t_ib1) (t_ib4+t_ib3+t_ib2+t_ib1)],[J1(4)/1000*x(t_ib1, iGh-iGA)*9806 J2(4)/1000*x((t_ib2+t_ib1), iGh-iGA)*9806 J3(4)/1000*x((t_ib3 + t_ib2+t_ib1), iGh-iGA)*9806 J4(4)/1000*x((t_ib4 + t_ib3 + t_ib2+t_ib1), iGh-iGA)*9806],'ko')
% ylabel('Po (watts)'),xlabel('Time (Seconds)');
% subplot(2,4,8),plot([t_ib1 (t_ib2+t_ib1) (t_ib3+t_ib2+t_ib1) (t_ib4+t_ib3+t_ib2+t_ib1)],[J1(Po) J2(Po) J3(Po) J4(Po)],'ko')
% ylabel('Po (watts)'),xlabel('Time (Seconds)');
