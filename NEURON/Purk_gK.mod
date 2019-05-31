TITLE Purkinje simplified Cell model

COMMENT
        K channel
        
        Authors: Niceto R. Luque & Richard R. Carrillo
		Last revised: 13.9.2018
		Paper related:Burst-pause Purkinje dynamics regulate motor adaptation
   
ENDCOMMENT
 
NEURON { 
	SUFFIX Purk_gK
	USEION k READ ek WRITE ik 
	RANGE gkbar, ik, g


	RANGE n_inf, tau_n, n
} 
 
UNITS { 
	(mA) = (milliamp) 
	(mV) = (millivolt) 
	(molar) = (1/liter)
	(mM) = (millimolar)
} 
 
PARAMETER { 

	v (mV) 
	gkbar= 0.01 (mho/cm2)
	ek = -95 (mV)
} 

STATE { 
	n 
} 

ASSIGNED { 
	ik (mA/cm2) 

	n_inf 
	tau_n (ms) 
	g (mho/cm2)
} 
 
INITIAL { 
	rate(v) 
	n = n_inf 
} 
 
BREAKPOINT { 
	SOLVE states METHOD derivimplicit 
	g = gkbar*n*n*n*n
	ik = g*(v - ek)
} 
 
DERIVATIVE states { 
	rate(v) 
	n' =(n_inf - n)/tau_n
} 
 
 
PROCEDURE rate(v (mV)) {
   if (v<=-10) {
      tau_n = 0.25 + 4.35 * exp((v + 10)/10)
   } else {
      tau_n = 0.25 + 4.35 * exp((-v - 10)/10)
   }

   n_inf = 1/(1 + exp((-v - 29.5)/10))
}
