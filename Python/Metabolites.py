import Global as g
import math

def metabolites(t, x, param):
	FAT = x[g.iFAT]
	CHO = x[g.iCHO]
	Gr  = x[g.iGr]
	Gh  = x[g.iGh]
	GA  = x[g.iGA]
	PCR = x[g.iPCR]

	# conductance of Beta oxidation pathway into matrix redox tank
	BOX    = param[0]
	# conductance of ETC (J/(Gr-Gh))
	ETC    = param[1]
	# conductance of Comp V + ANT (J/(Gh-GA))
	VA     = param[2]
	# Rate of ATP Breakdown (in amps, 1A ~ 10.4 umol ATP/sec)
	ATPase = param[3]
	# conductance of proton leak
	pleak  = param[4]
	# conductance of electron leak (J/Gr)
	eleak  = param[5]
	# high conductance pipe (CKase) connecting PCr and ATP tanks
	CK     = param[6]
	# VMax in m^3 H2O per sec
	Cpdh   = param[12]
	# ATP breakdown reaction of fixed ATP utilization (1220 ml H2O per sec)
	vCONST = param[14]

	# converts tank height to Pascals (N / m^3 water on earth)
	Pa = 9806
	# computes flux into fuel source
	J_fuel = 0
	# computes beta oxidation flux of lipid into matrix redox potential
	J_BOX = Pa*(FAT-Gr)*BOX
	# computes PDH flux of SHO into matrix redox potential
	Cpdh = param[12] #Conductance of glycolytic carbon entry into matrix redox
	CHOSat = param[13] #Gives status of CHO availability (0 - 1.0)
	nH = 3

	Act = ((8-GA)/(8-5))**nH #Delta Gatp related activation term of glycolysis

	# Inhib = (8-Gr)^0.25/(8-5)^0.25 #Delta G redox related feedback inhibition of PDH
	J_pdh = Pa*(CHO-Gr)* Cpdh * Act * CHOSat

	# Computing movement of electrons into proton motive force
	J_ETC = Pa*(Gr-Gh)*ETC

	# Computing flux of protons into ATP Synthesis & Export
	KmADP = 1
	nH = 2
	VAkin = (8-GA)**nH/((8-GA)**nH+KmADP**nH) #M-M relation with Hill nH and "KmADP" at 6.5 (~12 uM ADP)

	J_VA = Pa*(Gh-GA)*VA*VAkin

	# Computing flux of PCr into ATP
	J_CK = Pa*(PCR-GA)*CK

	# Computing ATP breakdown by energy sensitive ATPase
	nH = 0.2
	ki = ((GA-5)/(8-5))**nH

	J_ATPase = Pa*GA*ATPase*ki

	# Computing ATP breakdown by energy insensitive ATPase
	J_ATPconst = vCONST #This is the Vmax of the healthy 1 kg muscle (1220 ml H2O/sec)


	# Computing proton leak flux
	J_HL = pleak*math.exp(Gh)

	# Computing superoxide leak flux
	J_SO = Pa*eleak*Gr

	Af = param[7] #area of the fat fuel tank in m^2
	Ar = param[8] #area of the redox tank in m^2
	Ap = param[9] #area of Dp tank in m^2
	Aa = param[10] #area of ATP tank in m^2
	Apcr = param[11] #area of PCr tank in m^2

	# Computing Time Derivatives of Global Variables
	vector = []
	# fat
	vector.append( -J_BOX/Af )
	# CHO
	vector.append( 0 )
	# Gr
	vector.append( (J_BOX + J_pdh - J_ETC - J_SO)/Ar ) #(J_DH - J_ETC - J_SO)/param(11)
	# Gh
	vector.append( (J_ETC - J_VA - J_HL)/Ap ) #(J_ETC - J_VA - J_HL)/param(12)
	# GA
	vector.append( (J_VA + J_CK - J_ATPase - J_ATPconst)/Aa )
	# PCR
	vector.append( (-J_CK)/Apcr )

	return vector
