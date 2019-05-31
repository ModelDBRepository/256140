TITLE Purkinje simplified Cell model

COMMENT
        Leaky current
        
        Authors: Niceto R. Luque & Richard R. Carrillo
		Last revised: 13.9.2018
		Paper related:Burst-pause Purkinje dynamics regulate motor adaptation
ENDCOMMENT

NEURON {
	SUFFIX Purk_lkg
	NONSPECIFIC_CURRENT i
	RANGE el, glbar, i
}

UNITS {
	(mA) = (milliamp)
	(mV) = (millivolt)
}

PARAMETER {
	v (mV)
	glbar = 0.002 (mho/cm2)
	el = -70 (mV)
}

ASSIGNED { i (mA/cm2) }

BREAKPOINT { i = glbar * (v - el ) }
