TITLE BK-near channel
: Model from Fletcher et al. 2017 J Neurophysiol 117 (2298-2311)
: Originally to model pituitary corticotrophins in mice
: Used as BK channel for gonandotrope cell in pituitary gland of medaka.
: Implemented by G.Halnes March 2019

NEURON {
    SUFFIX bk
    USEION k READ ek WRITE ik VALENCE 1
    USEION Ca READ iCa VALENCE 2
    GLOBAL kninf
    RANGE gbk, ntau, AA
}

PARAMETER {
        gbk = 0.00 (mho/cm2)
        ek (mV)		: must be explicitely def. in hoc
        v (mV)
	iCa (mA)
        ntau = 20 (ms)
	AA = 1.21
}


UNITS {
        (mA) = (milliamp)
        (mV) = (millivolt)
}


ASSIGNED {
        ik (mA/cm2)
        kninf
	cdom
	vbk
}

STATE {n}

BREAKPOINT {
        SOLVE states METHOD cnexp
        ik = gbk*n*(v-ek)
}

INITIAL {
	trates(v, iCa)
}

DERIVATIVE states {
        trates(v,iCa)
        n' = (kninf-n)/ntau
}


PROCEDURE trates(v (mV), iCa (mM)) { :callable from hoc
	cdom = -AA*iCa
	vbk = 0.1-18*log(cdom/0.002)
        kninf = 1/(1+exp((vbk -v)/3) )
}