TITLE naxm
: Na current for medaka.
: G.Halnes Jun. 2017

NEURON {
	SUFFIX naxm
	USEION na READ ena WRITE ina
	RANGE  gbar
	GLOBAL minf, hinf, mtau, htau
}

PARAMETER {
	gbar = 0.07   	(mho/cm2)
	mmin=0.02
	hmin=0.1
	q10=2
	ena		(mV)            : must be explicitly def. in hoc
	celsius
	v 		(mV)
    vlj = 9 (mV)

:  Cell B fitting traub
    pm1 = 0.0382
    pm2 = -60.32 (mV)
    pm3 = 5.77 (mV)
    pm4 = 0.135
    pm5 = -26.17 (mV)
    pm6 = 2.9624e-5 (mV)
    pm7 = 28.84 (mV)
    pm8 = 4.55 (mV)

: Cell b comes from Matlab-file (slow recovery)
    ph1 = 0.0401
    ph2 = -32.4 (mV)
    ph3 = 3.29 (mV)
    ph4 = 2.65
    ph5 = -2145 (mV)
    ph6 = 139.3 (mV)
    ph7 = 55.0 (mV)
    ph8 = 5.07 (mV)
}

UNITS {
	(mA) = (milliamp)
	(mV) = (millivolt)
	(pS) = (picosiemens)
	(um) = (micron)
}

ASSIGNED {
	ina 		(mA/cm2)
	thegna		(mho/cm2)
	minf 		hinf
	mtau (ms)	htau (ms)
}


STATE { m h}

BREAKPOINT {
        SOLVE states METHOD cnexp
        thegna = gbar*m*m*m*h
	ina = thegna * (v - ena)
}

INITIAL {
	trates(v+vlj)
	m=minf
	h=hinf
}

DERIVATIVE states {
        trates(v+vlj)
        m' = (minf-m)/mtau
        h' = (hinf-h)/htau
}

PROCEDURE trates(vm) {
        LOCAL  a, b, qt
:        qt=q10^((celsius-24)/10)

        a = trap1(vm, pm1, pm2, pm3)
        b = trap1(-vm, pm4, -pm5, pm6)
        mtau = 1/(a+b)
        if (mtau<mmin) {mtau=mmin}

        minf = 1/(1+exp((-pm7-vm)/pm8))

        a = trap1(vm, ph1, ph2, ph3)
        b = trap1(-vm, ph4, -ph5, ph6)
        htau = 1/(a+b)
        if (htau<hmin) {htau=hmin}

        hinf = 1/(1+exp((vm+ph7)/ph8))

}


FUNCTION trap1(v,pa,pb,pc) {
	if (fabs(v-pb) > 1e-6) {
	        trap1 = pa * (pb-v) / (exp((pb - v)/pc)-1)
	} else {
	        trap1 = pa*pc
 	}
}







