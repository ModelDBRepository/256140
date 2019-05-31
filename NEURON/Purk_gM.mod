TITLE Purkinje simplified Cell model

COMMENT
        kM channel
        
        Authors: Niceto R. Luque & Richard R. Carrillo
		Last revised: 13.9.2018
		Paper related:Burst-pause Purkinje dynamics regulate motor adaptation
   
ENDCOMMENT
 
NEURON { 
	SUFFIX Purk_gKM
	USEION k READ ek WRITE ik
	RANGE gkmbar, ik, g


	RANGE M_inf, tau_M, M
} 
 
UNITS { 
	(mA) = (milliamp) 
	(mV) = (millivolt) 
	(molar) = (1/liter)
	(mM) = (millimolar)
} 
 
PARAMETER { 

	v (mV) 
	gkmbar= 0.001 (mho/cm2)
	ek = 125 (mV)
} 

STATE { 
	M
} 

ASSIGNED { 
	ik (mA/cm2) 

	M_inf 
	tau_M (ms) 
	g (mho/cm2)
	alpha_M (/ms) 
	beta_M (/ms)
} 
 
INITIAL { 
	rate(v) 
	M = M_inf 
} 
 
BREAKPOINT { 
	SOLVE states METHOD derivimplicit 
	g = gkmbar*M
	ik = g*(v - ek)
	alpha_M = alp_M(v) 
	beta_M = bet_M(v)
} 
 
DERIVATIVE states { 
	rate(v) 
	M' =(M_inf - M)/tau_M
} 
 
FUNCTION alp_M(v(mV))(/ms) {
	alp_M = 0.02/(1 + exp((-v - 20)/5))
} 
 
FUNCTION bet_M(v(mV))(/ms) {
	bet_M = 0.01*exp((-v - 43)/18)
} 

PROCEDURE rate(v (mV)) { LOCAL a_M, b_M
	a_M = alp_M(v)  
	b_M = bet_M(v) 
	tau_M = 1/(a_M + b_M)
	M_inf = a_M/(a_M + b_M) 
}
