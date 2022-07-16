%Hydraulic model simulator includes 2 fuel sources, 6 tanks (fat, cho, Gr, Gh, ATP, PCR), and vATP that is energy sensitive 
clear; close all;
global iFAT iCHO iGr iGh iGA iPCR

Indices();

% Load Initial Conditions
pathname = 'C:\Users\Ben\Desktop\ANT Model Files\Summer 2017 ANT Model Files\Hydraulic Model\';
filename = 'HydraulicLoading.xlsx';
[d] = xlsread([filename],'C3:C8');

xo = d;

% Load Parameters
[p] = xlsread([filename],'H3:H17');
param = p;

indexC = [1 3:6];
options = odeset('MaxStep',5e-2);

% Below are the intial water levels
xo(iFAT) = 8;
xo(iCHO) = 11.5;
xo(iGr) = 8;
xo(iGh) = 8;
xo(iGA) = 8;
xo(iPCR) = 8;

vCONST = param(15); %Constant ATP breakdown reaction is zero in this script
param(15) = 0*vCONST;

CHOSat = param(14); %CHO status (values 0 - 1.0) gives glucose & glycogen availability
param(14) = 1*CHOSat;

Vpdh = param(13); % conductance of the glycolytic pathway
param(13) = 1*Vpdh;
% Below are the conductances of each pipe
Cbox = param(1);  % this is the conductance of the beta oxidation pathway
param(1) = 1*Cbox;

Cetc = param(2);  % this is the conductance of ETC
param(2) = 1*Cetc;

Cva = param(3);  % this is the conductance of Complex V + ANT
param(3) = 1*Cva;

CCK = param(7);  % this is the CKase pipe
param(7) = 1*CCK;

CHL = param(5);  % this is conductance factor of proton leak
param(5) = 1*CHL;

CSO = param(6);  % this is conductance factor of superoxide electron leak
param(6) = 1*CSO;

Xatp = param(11); %cross-sectional area of PCr tank. Each x-fold adjusts [TCr] by that factor
param(11) = 1*Xatp;

Xpcr = param(12); %cross-sectional area of PCr tank. Each x-fold adjusts [TCr] by that factor
param(12) = 1*Xpcr;


vatp = param(4); %Conductance of the energy sensitive ATP breakdown reaction
param(4) = 0.05*vatp;

tic
t_ib1 = 600; % resting metabolism
[t,x] = ode15s(@Metabolites,[0 t_ib1],xo,options,param);
toc

J1 = Fluxes(x(end,:),param);
xo = x(end,:);
x1 = x;
t1 = t;

param(4) = 0.5*vatp; % onset of first level of ATP breakdown
% xo(iPCR) = 8;
% param(1) = 4.64*Cdh;
% param(2) = 1*Cetc;
% param(3) = 1.8*Cva;
tic
t_ib2 = 300; % exercise time
[t,x] = ode15s(@Metabolites,[0 t_ib2],xo,options,param);
toc

J2 = Fluxes(x(end,:),param);
xo = x(end,:);
x2 = x;
t2 = t;

param(4) = .75*vatp; % onset of 2nd level of ATP breakdown
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

% param(1) = 4.64*Cdh;
% param(2) = .9*Cetc;
% param(3) = 2.01*Cva;
param(4) = 1*vatp; % onset of third level of ATP breakdown
% xo(iPCR) = 8;


tic
t_ib4 = 300; % exercise time
[t,x] = ode15s(@Metabolites,[0 t_ib3],xo,options,param);
toc

J4 = Fluxes(x(end,:),param);
xo = x(end,:);
x4 = x;
t4 = t;

param(4) = 1.5*vatp; %

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


%Po = zeros(length(x));
%P%a = 9806;
%f%or i = 1:length(x)-1
%    Po(i) = Pa*param(4)*x(iGA,i+1)-Pa*param(4)*x(iGA,i);
%end

%Defining array for download into Excel

tt1=t';
y=tt1(1:30:end,1);
xlswrite('kindataNEW.xlsx',y,1,'A3:A811');

variables=x(1:30:end,4);  %variables=x(1:200:end,[25 35 36 62 92]); %this is the full list of variables if wanted and expand Column letters accordingly e.g., 'B3:C410' for one additional variable
xlswrite('kindataNEW.xlsx',variables,1,'B3:B1600');

%Plotting
Xatp = param(11);
Xpcr = param(12);

figure('NumberTitle','off','Name','CKase & OxPhos flux')
xx=cat(1,0,-diff(x(:,6)));
tt=cat(1,0,diff(t'));
subplot(2,2,1),plot(t,x(:,6))
ylabel('[PCr]')
set(gca,'xlim',[0 (t_ib5+t_ib4+t_ib3+t_ib2+t_ib1)],'ylim',[5  8]);
subplot(2,2,2),plot(t,xx./tt*1e6*Xpcr)
ylabel('J_C_K (ml H_2O/sec)')
set(gca,'xlim',[0 (t_ib5+t_ib4+t_ib3+t_ib2+t_ib1)],'ylim',[-1000  1000]);
xx=cat(1,0,-diff(x(:,5)));
tt=cat(1,0,diff(t'));
subplot(2,2,3),plot(t,x(:,5))
ylabel('[ATP]')
set(gca,'xlim',[0 (t_ib5+t_ib4+t_ib3+t_ib2+t_ib1)],'ylim',[5  8]);
subplot(2,2,4),plot(t,xx./tt*1e6*Xatp)
ylabel('Net ATP change (ml H_2O/sec)')
set(gca,'xlim',[0 (t_ib5+t_ib4+t_ib3+t_ib2+t_ib1)],'ylim',[-1000  1000]);


figure %Water levels in each tank
subplot(1,5,1),plot(t,x(:,1));
ylabel('Lipid')
set(gca,'xlim',[0 (t_ib5+t_ib4+t_ib3+t_ib2+t_ib1)],'ylim',[5 8.5]);
subplot(1,5,2),plot(t,x(:,3));
ylabel('Redox')
set(gca,'xlim',[0 (t_ib5+t_ib4+t_ib3+t_ib2+t_ib1)],'ylim',[5 8.5]);
subplot(1,5,3),plot(t,x(:,4));
ylabel('Delta p')
set(gca,'xlim',[0 (t_ib5+t_ib4+t_ib3+t_ib2+t_ib1)],'ylim',[5 8.5]);
subplot(1,5,4),plot(t,x(:,5));
ylabel('ATP')
set(gca,'xlim',[0 (t_ib5+t_ib4+t_ib3+t_ib2+t_ib1)],'ylim',[5 8.5]);
subplot(1,5,5),plot(t,x(:,6));
ylabel('PCr')
set(gca,'xlim',[0 (t_ib5+t_ib4+t_ib3+t_ib2+t_ib1)],'ylim',[5 8.5]);


figure %Fluxes
subplot(2,4,1),plot([t_ib1 (t_ib2+t_ib1) (t_ib3+t_ib2+t_ib1) (t_ib4+t_ib3+t_ib2+t_ib1) (t_ib5+t_ib4+t_ib3+t_ib2+t_ib1)],[J1(2) J2(2) J3(2) J4(2) J5(2)],'ko');
ylabel('BetaOx (mL H2O/sec)'),xlabel('Time (Seconds)')
subplot(2,4,2),plot([t_ib1 (t_ib2+t_ib1) (t_ib3+t_ib2+t_ib1) (t_ib4+t_ib3+t_ib2+t_ib1) (t_ib5+t_ib4+t_ib3+t_ib2+t_ib1)],[J1(3) J2(3) J3(3) J4(3) J5(3)],'ko');
ylabel('PDH (mL H2O/sec)'),xlabel('Time (Seconds)')

subplot(2,4,3),plot([t_ib1 (t_ib2+t_ib1) (t_ib3+t_ib2+t_ib1) (t_ib4+t_ib3+t_ib2+t_ib1) (t_ib5+t_ib4+t_ib3+t_ib2+t_ib1)],[J1(4) J2(4) J3(4) J4(4) J5(4)],'ko');
ylabel('ETC (mL H2O/sec)'),xlabel('Time (Seconds)')
subplot(2,4,4),plot([t_ib1 (t_ib2+t_ib1) (t_ib3+t_ib2+t_ib1) (t_ib4+t_ib3+t_ib2+t_ib1) (t_ib5+t_ib4+t_ib3+t_ib2+t_ib1)],[J1(5) J2(5) J3(5) J4(5) J5(5)],'ko');
ylabel('Comp V + ANT (mL H2O/sec)'),xlabel('Time (Seconds)')
subplot(2,4,5),plot([t_ib1 (t_ib2+t_ib1) (t_ib3+t_ib2+t_ib1) (t_ib4+t_ib3+t_ib2+t_ib1) (t_ib5+t_ib4+t_ib3+t_ib2+t_ib1)],[J1(9) J2(9) J3(9) J4(9) J5(9)],'ko');
ylabel('CKase(mL H2O/sec)'),xlabel('Time (Seconds)')


subplot(2,4,6),plot([t_ib1 (t_ib2+t_ib1) (t_ib3+t_ib2+t_ib1) (t_ib4+t_ib3+t_ib2+t_ib1) (t_ib5+t_ib4+t_ib3+t_ib2+t_ib1)],[J1(6) J2(6) J3(6) J4(6) J5(6)],'ko');
ylabel('H Leak (mL H2O/sec)'),xlabel('Time (Seconds)')
subplot(2,4,7),plot([t_ib1 (t_ib2+t_ib1) (t_ib3+t_ib2+t_ib1) (t_ib4+t_ib3+t_ib2+t_ib1) (t_ib5+t_ib4+t_ib3+t_ib2+t_ib1)],[J1(7) J2(7) J3(7) J4(7) J5(7)],'ko');
ylabel('Superoxide (mL H2O/sec)'),xlabel('Time (Seconds)')
subplot(2,4,8),plot([t_ib1 (t_ib2+t_ib1) (t_ib3+t_ib2+t_ib1) (t_ib4+t_ib3+t_ib2+t_ib1) (t_ib5+t_ib4+t_ib3+t_ib2+t_ib1)],[J1(8) J2(8) J3(8) J4(8) J5(8)],'ko');
ylabel('ATP(mL H2O/sec)'),xlabel('Time (Seconds)');
% subplot(2,4,8),plot([t_ib1 (t_ib2+t_ib1) (t_ib3+t_ib2+t_ib1) (t_ib4+t_ib3+t_ib2+t_ib1)],[J1(4)/1000*x(t_ib1, iGh-iGA)*9806 J2(4)/1000*x((t_ib2+t_ib1), iGh-iGA)*9806 J3(4)/1000*x((t_ib3 + t_ib2+t_ib1), iGh-iGA)*9806 J4(4)/1000*x((t_ib4 + t_ib3 + t_ib2+t_ib1), iGh-iGA)*9806],'ko')
% ylabel('Po (watts)'),xlabel('Time (Seconds)');
% subplot(2,4,8),plot([t_ib1 (t_ib2+t_ib1) (t_ib3+t_ib2+t_ib1) (t_ib4+t_ib3+t_ib2+t_ib1)],[J1(Po) J2(Po) J3(Po) J4(Po)],'ko')
% ylabel('Po (watts)'),xlabel('Time (Seconds)');
