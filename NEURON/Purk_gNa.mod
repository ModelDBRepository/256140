TITLE Purkinje simplified Cell model

COMMENT
        Na channel
        
        Authors: Niceto R. Luque & Richard R. Carrillo
		Last revised: 13.9.2018
		Paper related:Burst-pause Purkinje dynamics regulate motor adaptation
   
ENDCOMMENT
 
NEURON { 
	SUFFIX Purk_gNa
	USEION na READ ena WRITE ina
	RANGE gnabar, ina, g


	RANGE h_inf, m_inf, tau_h, h
} 
 
UNITS { 
	(mA) = (milliamp) 
	(mV) = (millivolt) 
	(molar) = (1/liter)
	(mM) = (millimolar)
} 
 
PARAMETER { 

	v (mV) 
	gnabar= 0.125 (mho/cm2)
	ena = 50 (mV)
} 

STATE { 
	h
} 

ASSIGNED { 
	ina (mA/cm2) 

	h_inf
	m_inf
	tau_h (ms)
	g (mho/cm2)
} 
 
INITIAL { 
	rate(v) 
	h = h_inf 
} 
 
BREAKPOINT { 
	SOLVE states METHOD derivimplicit 
	m_inf = 1/(1 + exp((-v - 48)/10))
	g = gnabar*m_inf*m_inf*m_inf*h 
	ina = g*(v - ena)	
} 
 
DERIVATIVE states { 
	rate(v) 
	h' =(h_inf - h)/tau_h
} 
 
PROCEDURE rate(v (mV)) {
   tau_h = 0.15 + 1.15/(1 + exp((v + 33.5)/15))
   h_inf = 1/(1 + exp((v + 59.4)/10.7))
}
