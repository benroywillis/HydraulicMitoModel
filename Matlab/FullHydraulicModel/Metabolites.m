function [f] = Metabolites(t,x,param)

% Rate Equations and Time Derivatives for the "simple hydraulic" model.

% Defining Global Variables
global iFAT iCHO iGr iGh iGA iPCR

% Defining local variables
FAT     = x(iFAT);
CHO     = x(iCHO);
Gr      = x(iGr);
Gh      = x(iGh);
GA      = x(iGA);
PCR     = x(iPCR);


%Defining parameters (Vmax of DH & proton leak; conductances of DH, ETC, VA; velocity of ATPase)
BOX     = param(1); %conductance of the Beta oxidation pathway into matrix redox tank
ETC     = param(2); %conductance of the ETC (= J/(Gr-Gh))
VA      = param(3); %conductance of Comp V + ANT (= J/(Gh-GA))
ATPase  = param(4); %Rate of ATP breakdown (in amps and 1.0 amp ~ 10.4 umol ATP/sec)
pleak   = param(5); %conductance of proton leak
eleak   = param(6); %conductance of electron leak (= J/Gr)
CK      = param(7); %high conductance pipe (CKase) connecting PCr and ATP tanks
Cpdh    = param(13); %Vmax in m^3 H2O per sec
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

% Computing ATP breakdown by energy sensitive ATPase
nH = 0.2;
ki = ((GA-5)/(8-5))^nH;

J_ATPase = Pa*GA*ATPase*ki;

% Computing ATP breakdown by energy insensitive ATPase
J_ATPconst = vCONST; %This is the Vmax of the healthy 1 kg muscle (1220 ml H2O/sec)


% Computing proton leak flux

J_HL = pleak*exp(Gh);

% Computing superoxide leak flux

J_SO = Pa*eleak*Gr;

Af = param(8); %area of the fat fuel tank in m^2
Ar = param(9); %area of the redox tank in m^2
Ap = param(10); %area of Dp tank in m^2
Aa = param(11); %area of ATP tank in m^2
Apcr = param(12); %area of PCr tank in m^2

% Computing Time Derivatives of Global Variables
f(iFAT)     = -J_BOX/Af;
f(iGr)      = (J_BOX + J_pdh - J_ETC - J_SO)/Ar; %(J_DH - J_ETC - J_SO)/param(11);
f(iGh)      = (J_ETC - J_VA - J_HL)/Ap; %(J_ETC - J_VA - J_HL)/param(12);
f(iGA)      = (J_VA + J_CK - J_ATPase - J_ATPconst)/Aa;
f(iPCR)     = (-J_CK)/Apcr;

f = f';


