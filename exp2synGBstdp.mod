COMMENT
NEURON module of Graupner&Brunel_PNAS_2012 STDP learning rule, implemented
for double exponential synapse.

Midnight project of
Ruben A. Tikidji-Hamburyan <rth@nisms.krinc.ru>
ENDCOMMENT

NEURON {
	POINT_PROCESS exp2synGBstdp
	: Synapse
	RANGE tau1, tau2, e, g, i 
	: STDP
	RANGE gammad, gammap, thetap, thetad, taurho, rho12, rhoinit
	: Ca
	RANGE tauca, cpre, cpost, cdelay
	: STDP FLAGs
	RANGE bistable, plastic
	NONSPECIFIC_CURRENT i
}

UNITS {
	(nA) = (nanoamp)
	(mV) = (millivolt)
	(uS) = (microsiemens)
}

PARAMETER {
	:Synaptic Parameters
	tau1 = 0.8 (ms) <1e-9,1e9>
	tau2 = 2 (ms) <1e-9,1e9>
	e = 0 (mV)
	: STDP parameters (see Table S1 in supplementary Tables G&B 2012) 
	gammad = 200 <5,50000>				: gain of depression
	gammap = 321.808 <5,25000>			: gain of potentiation
	thetad = 1 <1e-9,1e9>				: threshold of depression
	thetap = 1.3 <1e-9,1e9>				: threshold of potentiation
	taurho = 150e3 (ms) <25e3,2500e3>	: time constant of STDP 
	rho12    = 0.5						: middle point of cubic polynomial term
	rhoinit  = 1						: initial condition for STDP
	: Ca dynamic parameters
	tauca  = 20 (ms) <1,100>			: calcium time constant 
	cpre   = 1 <0.1,20>					: calcium influx evoked by presynaptic spike
	cpost  = 2 <0.1,50>					: calcium influx evoked by postsynaptic spike
	cdelay = 13.7 (ms) <0,50>			: delay of presynapticall evoked calcium influx 
	: STDP FLAGs
	bistable = 1 <0,1>					: FLAG, if 0 - STDP is not bistable, if 1 - STDP has bistability.
	plastic = 1 <0,1>					: FLAG if 1 - synapse is plastic, i.e. STDP is active.
}

ASSIGNED {
	v (mV)
	i (nA)
	g (nS)
	factor
	lud (ms) : last update variable keep track of time intervals between
			 : STDP updates.
}

STATE {
	A (uS)
	B (uS)
}

INITIAL {
	LOCAL tp
	if (tau1/tau2 > .9999) {
		tau1 = .9999*tau2
	}
	A = 0
	B = 0
	tp = (tau1*tau2)/(tau2 - tau1) * log(tau2/tau1)
	factor = -exp(-tp/tau1) + exp(-tp/tau2)
	factor = 1/factor
	net_send(0, 1)
}

BREAKPOINT {
	SOLVE state METHOD cnexp
	g = B -A
	i = g*(v - e)
	STDPupdate( dt )
}

DERIVATIVE state {
  A' = -A/tau1
  B' = -B/tau2 
}



: to catch the list of NetCon(s) to update variables dynamically
VERBATIM
	int nc_cnt, nc_number = 0;
	double* nc_vector, **nc_vectorhandler = NULL;
	double t_lock=-1e9;
ENDVERBATIM


PROCEDURE STDPupdate( time_step_dt (ms) ){ : calculated dynamical variables in netcon(s)
	VERBATIM
	if ( nc_vectorhandler != NULL && nc_number > 0 && fabs(t_lock-t)>1e-6) {
		//DB>>
		fprintf(stderr,"nc_vectorhandler is not zero. Number of netcon(s) is %d,time step is %g, time %g\n",nc_number,_ltime_step_dt,t);
		//<<DB
		for(nc_cnt=0, nc_vector = nc_vectorhandler[0]; nc_cnt < nc_number; nc_vector = nc_vectorhandler[++nc_cnt]){
			double rho = nc_vector[1], conc = nc_vector[2];
			nc_vector[1] += (
				-bistable*rho*(1.-rho)*(rho12-rho)			//STDP bistability
				+gammap*(1.-rho)*(float)((conc-thetap)>=0)	//Potentiation
				-gammad*rho*(float)((conc-thetad)>=0)		//Depression
			)*_ltime_step_dt/taurho;
			nc_vector[2]  = conc * exp(-_ltime_step_dt/tauca);
		}
		t_lock = t;
	}
	ENDVERBATIM

}

: w    intrinsic synaptic weight
: rho  plasticity factor (see original paper)
: c    calcium concentration
NET_RECEIVE(w (uS), rho, c) {
	INITIAL { rho = rhoinit  c = 0 }
VERBATIM
		if ( nc_vectorhandler == NULL || nc_number <= 0 ){
			nc_number = _nrn_netcon_args(_ppvar[_fnc_index]._pvoid, &nc_vectorhandler);
			//DB>>
			fprintf(stderr,"Set up nc_number in %d\n",nc_number);
			//<<DB
		}
ENDVERBATIM
	if (flag == 0) { : presynaptic spike (after last post so depress)
		net_send(cdelay, 3) : wait delay period
		if (plastic) {
			A = A + factor*w*rho
			B = B + factor*w*rho
		} else { : plasticity is turned off
			A = A + factor*w
			B = B + factor*w
		} 
	}else if (flag == 2) { : postsynaptic spike (after last pre so potentiate)
		FOR_NETCONS(w1, rho1, c1) {
			c1 = c1 + cpost
		}
	}else if (flag == 3) { : delay after presynaptic spike
		c = c + cpre
	} else { : flag == 1 from INITIAL block
		WATCH (v > -20) 2
	}
}

