TITLE High threshold calcium current
:
:   Ca++ current, presumably dominated by L type channels
:
:   Adapted from Model of Huguenard & McCormick, J Neurophysiol, 1992
:   Formalism of Goldman-Hodgkin-Katz
:
:   Fitted to data from gonadotrope cells in Medaka
:   Written by Geir Halnes, Norwegian University of Life Sciences, June 2018


INDEPENDENT {t FROM 0 TO 1 WITH 1 (ms)}

NEURON {
	SUFFIX ihva
	USEION Ca READ Cai, Cao WRITE iCa VALENCE 2
      RANGE pcabar, g
	GLOBAL 	m_inf, taum
}


UNITS {
	(mA) = (milliamp)
	(mV) = (millivolt)
	(molar) = (1/liter)
	(mM) = (millimolar)
	FARADAY = (faraday) (coulomb)
	R = (k-mole) (joule/degC)
}


PARAMETER {
	v		(mV)
	celsius	= 36	(degC)
	eCa     = 60		(mV)
	Cai 	= .00005	(mM)	: initial [Ca]i = 50 nM
	Cao 	= 2		(mM)	: [Ca]o = 2 mM
	pcabar	= 2e-4	(mho/cm2)
    mmin=0.02
    hmin=0.1

    :  Cell B fitting traub
    pm1 = -0.1278
    pm2 = -46.7 (mV)
    pm3 = 19.0 (mV)
    pm4 = -101.54
    pm5 = 535.1 (mV)
    pm6 = -60.0 (mV)
    pm7 = -6.79 (mV)
    pm8 = 6.57 (mV)
    vlj = 15 (mV)
}

STATE {
	m
}

INITIAL {
:	tadj = 3 ^ ((celsius-21.0)/10)
	trates(v+vlj)
	m = m_inf
}


ASSIGNED {
	iCa	(mA/cm2)
	g       (mho/cm2)
	m_inf
	taum	(ms)
}

BREAKPOINT {
	SOLVE states METHOD cnexp
	g = pcabar * m * m
	iCa = g * ghk(v, Cai, Cao)
}

DERIVATIVE states {
	trates(v+vlj)
	m' = (m_inf - m) / taum
}


UNITSOFF

PROCEDURE trates(vm) {
        LOCAL  a, b, qt
:        qt=q10^((celsius-24)/10)

        a = trap1(vm, pm1, pm2, pm3)
        b = trap1(-vm, pm4, -pm5, pm6)
        taum = 1/(a+b)/100000
        if (taum<mmin) {taum=mmin}
        m_inf = 1/(1+exp((pm7-vm)/pm8))
}


FUNCTION ghk(v(mV), ci(mM), co(mM)) (.001 coul/cm3) {
	LOCAL z, eci, eco
	z = (1e-3)*2*FARADAY*v/(R*(celsius+273.15))
	eco = co*efun(z)
	eci = ci*efun(-z)
	:high co charge moves inward
	:negative potential charge moves inward
	ghk = (.001)*2*FARADAY*(eci - eco)
}

FUNCTION efun(z) {
	if (fabs(z) < 1e-4) {
		efun = 1 - z/2
	}else{
		efun = z/(exp(z) - 1)
	}
}

FUNCTION trap1(v,pa,pb,pc) {
	if (fabs(v-pb) > 1e-6) {
	        trap1 = pa * (pb-v) / (exp((pb - v)/pc)-1)
	} else {
	        trap1 = pa*pc
 	}
}
UNITSON


UNITSON
