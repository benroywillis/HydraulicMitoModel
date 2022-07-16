function [J] = Fluxes(x,param)

% This function outputs the discrete time flux values for the hydraulic model.

% Defining Global Variables
global iFAT iCHO iGr iGh iGA iPCR

% Defining local variables
FAT = x(iFAT);
CHO = x(iCHO);
Gr = x(iGr);
Gh = x(iGh);
GA = x(iGA);
PCR = x(iPCR);


%Defining conductances or Vmax (Vmax for PDH only, see below)
BOX = param(1); %Beta oxidation pathway conductance from FAT to matrix redox
Vpdh = param(13); %PDH Vmax expressed in m^3 H2O per sec
ETC     = param(2);
VA      = param(3);
ATPase  = param(4); %this ATP breakdown reaction is inhibited by falling DGatp
pleak   = param(5);
eleak   = param(6);
CK      = param(7);
vCONST  = param(15); %this ATP breakdown reaction has one constant rate of ATP utilization (1220 mlH2O/sec)

Pa = 9806; %converts tank height to Pa of driving pressure
% Computing flux into fuel source
J_fuel = 0;

% Computing beta oxidation flux of lipid into matrix redox potential
J_BOX = Pa*(FAT-Gr)*BOX;

% Computing PDH flux of CHO into matrix redox potential
Cpdh = param(13); %Conductance of glycolytic carbon entry into matrix redox
CHOSat = param(14); %Gives status of CHO availability (0 - 1.0)
nH = 3;
Act = ((8-GA)/(8-5))^nH; %Delta Gatp related activation term of glycolysis

% Inhib = (8-Gr)^0.25/(8-5)^0.25; %Delta G redox related feedback inhibition of PDH

J_pdh = Pa*(CHO-Gr)* Cpdh * Act * CHOSat;

% Computing movement of electrons into proton motive force
J_ETC = Pa*(Gr-Gh)*ETC;

% Computing flux of protons into ATP Synthesis & Export
KmADP = 1;
nH = 2;
VAkin = (8-GA)^nH/((8-GA)^nH+KmADP^nH); %M-M relation with Hill nH and "KmADP" at 6.5 (~12 uM ADP)

J_VA = Pa*(Gh-GA)*VA*VAkin;

% Computing flux of PCr into ATP
J_CK = Pa*(PCR-GA)*CK;

% Computing ATP breakdown rate by energy sensitive ATPase
nH = 0.2;
ki = ((GA-5)/(8-5))^nH;

J_ATPase = Pa*GA*ATPase*ki;

% Computing ATP breakdown by energy insensitive ATPase
J_ATPconst = vCONST; %This is the Vmax of the healthy 1 kg muscle (1220 ml H2O/sec)


% Computing proton leak flux

J_HL = pleak*exp(Gh);
    
% Computing superoxide leak flux
J_SO = Pa*eleak*Gr;

% Outputting Fluxes
J(1)  = J_fuel*1000000; %*1000000 converts m^3 H2O to mL H2O
J(2)  = J_BOX*1000000; %*1000000 converts m^3 H2O to mL H2O
J(3)  = J_pdh*1000000;
J(4)  = J_ETC*1000000;
J(5)  = J_VA*1000000;
J(6)  = J_HL*1000000;
J(7)  = J_SO*1000000;
J(8)  = J_ATPase*1000000;
J(9)  = J_CK*1000000;
J(10) = J_ATPconst*1000000;

