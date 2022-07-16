function [J] = Fluxes(x,param)

% This function outputs the discrete time flux values for the simple 4 Tank Hydraulic Model: Fuel - Mito (Redox + DY) - ATP - PCr

% Defining Global Variables
global iFuel iMito iATP iPCR

% Defining local variables
Fuel     = x(iFuel);
Mito     = x(iMito);
ATP      = x(iATP);
PCR      = x(iPCR);


%Defining parameters (Vmax of DH & proton leak; conductances of DH, ETC, VA; velocity of ATPase)
DH  = param(1); %conductance of the fuel oxidation pathway into mito tank
OxPhos = param(2); %conductance of the ETC (= J/(Gr-Gh))
ATPase  = param(3); %conductance of ATP breakdown step)
CK      = param(4); %high conductance pipe (CKase) connecting PCr and ATP tanks
ATPconst = param(9); %ATP breakdown at constant rate insensitive to energy state

Pa = 9806; %converts tank height to Pa of driving pressure
% Computing flux into fuel source

% Computing fuel oxidation flux into mito tank
J_DH = Pa*(Fuel-Mito)*DH;

% Computing flux from mito to cyto ATP via OxPhos path
J_OxPhos = Pa*(Mito-ATP)*OxPhos;

% Computing CKase flux of PCr into ATP;
J_CK = Pa*(PCR-ATP)*CK;

% Computing ATP utilization rate by energetically sensitive ATPase
% nH = 0.2;
% ki = ((GA-5)/(8-5))^nH;

J_ATPase = Pa*ATP*ATPase; %insert ki factor in rate equation if needed

% Computing ATP utilization rate by insensitive (constant flux) ATP
% breakdown reaction

J_ATPconst = ATPconst; %This is the Vmax of the healthy 1 kg muscle (1220 ml H2O/sec)


Af = param(5); %area of the fat fuel tank in m^2
Am = param(6); %area of the redox tank in m^2
Aatp = param(7); %area of Dp tank in m^2
Apcr = param(8); %area of PCr tank in m^2

% Computing Time Derivatives of Global Variables
f(iFuel)    = -J_DH/Af;
f(iMito)    = (J_DH - J_OxPhos)/Am; %
f(iATP)     = (J_OxPhos - J_ATPase - J_ATPconst + J_CK)/Aatp; %
f(iPCR)     = (-J_CK)/Apcr;

f = f';

% Outputting Fluxes
J(1) = J_DH*1000000; %*1000000 converts m^3 H2O to mL H2O
J(2) = J_OxPhos*1000000;
J(3) = J_ATPase*1000000;
J(4) = J_CK*1000000;
J(5) = J_ATPconst*1000000;


