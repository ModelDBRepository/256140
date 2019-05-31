TITLE Purkinje simplified Cell model

COMMENT
        Ca channel
        
        Authors: Niceto R. Luque & Richard R. Carrillo
		Last revised: 13.9.2018
		Paper related:Burst-pause Purkinje dynamics regulate motor adaptation
   
ENDCOMMENT
 
NEURON { 
	SUFFIX Purk_gCa
	USEION ca READ eca WRITE ica
	RANGE gcabar, ica, g


	RANGE c_inf, tau_c, c
} 
 
UNITS { 
	(mA) = (milliamp) 
	(mV) = (millivolt) 
	(molar) = (1/liter)
	(mM) = (millimolar)
} 
 
PARAMETER { 

	v (mV) 
	gcabar= 0.001 (mho/cm2)
	eca = 125 (mV)
} 

STATE { 
	c
} 

ASSIGNED { 
	ica (mA/cm2) 

	c_inf 
	tau_c (ms) 
	g (mho/cm2)
	alpha_c (/ms) 
	beta_c (/ms)
} 
 
INITIAL { 
	rate(v) 
	c = c_inf 
} 
 
BREAKPOINT { 
	SOLVE states METHOD derivimplicit 
	g = gcabar*c*c
	ica = g*(v - eca)
	alpha_c = alp_c(v) 
	beta_c = bet_c(v)
} 
 
DERIVATIVE states { 
	rate(v) 
	c' =(c_inf - c)/tau_c
} 
 
FUNCTION alp_c(v(mV))(/ms) {
	alp_c = 1.6/(1 + exp(-0.072(v - 5)))
} 
 
FUNCTION bet_c(v(mV))(/ms) {
	bet_c = 0.02 * (v + 8.9)/(exp((v + 8.9)/5) - 1)
} 

PROCEDURE rate(v (mV)) { LOCAL a_c, b_c
	a_c = alp_c(v)  
	b_c = bet_c(v) 
	tau_c = 1/(a_c + b_c) 
	c_inf = a_c/(a_c + b_c) 
}
