COMMENT
NEURON module of Calcium-based plasticity model.
Graupner and Brunel PNAS 2012

Uses double exponential synapse model as a basis.

WARNING! It is a pre-alpha version. The module  does not produce curves 
         from the paper. The shape of curves are similar, but amplitudes
         are smaller than in the paper.
         The module does not have noise component. Noise will be added 
         later, when deterministic part will be debugged.
         The module will be revise after author's response.

RTH Midnight Project by
Ruben A. Tikidji-Hamburyan <rth@nisms.krinc.ru>
May 2015
ENDCOMMENT

NEURON {
    POINT_PROCESS exp2synGBstdp
    NONSPECIFIC_CURRENT i
    : Synapse
    RANGE tau1, tau2, e, g, i 
    : STDP
    RANGE gammad, gammap, thetap, thetad, taurho, rho12, rhoinit
    : Ca
    RANGE tauca, cpre, cpost, cdelay
    : STDP FLAGs
    RANGE bistable, plastic, learningdependence
    : Time step for STDP equations
    RANGE STDPstep
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
    gammad = 200 <5,50000>              : gain of depression
    gammap = 321.808 <5,25000>          : gain of potentiation
    thetad = 1 <1e-9,1e9>               : threshold of depression
    thetap = 1.3 <1e-9,1e9>             : threshold of potentiation
    taurho = 150e3 (ms) <25e3,2500e3>   : time constant of STDP 
    rho12    = 0.5                      : middle point of cubic polynomial term
    : Ca dynamic parameters
    tauca  = 20 (ms) <1,100>            : calcium time constant 
    cpre   = 1 <0.1,20>                 : calcium influx evoked by presynaptic spike
    cpost  = 2 <0.1,50>                 : calcium influx evoked by postsynaptic spike
    cdelay = 13.7 (ms) <0,50>           : delay of presynapticall evoked calcium influx 
    : STDP FLAGs
    bistable = 1 <0,1>                  : FLAG if 0 - STDP is not bistable, if 1 - STDP has bistability.
    plastic = 1 <0,1>                   : FLAG if 1 - synapse is plastic, i.e. STDP is active.
    learningdependence = 0 <0,1>        : FLAG if 1 - presynaptic calcium influx will depends on rho.
                                        : The reasons is simple: learning modulates current through NMDA and AMPA
                                        : and change calcium influx proportionally. In this case cpre is maximal 
                                        : possible calcium influx for this synapse. If 0 works as in the paper
    : Time step for STDP equations
    STDPstep = 0.25 (ms) <1e-9,1e9>     : time step for solution of STDP equations
}

ASSIGNED {
    v (mV)
    i (nA)
    g (nS)
    factor (1)
    lndt (ms)
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
    net_send(0, 1)          :init now
    net_send(STDPstep, 4)   :solve equation in fixed time step
    lndt = t
}

BREAKPOINT {
    SOLVE state METHOD cnexp
    g = B-A
    i = g*(v - e)
}


DERIVATIVE state {
  A' = -A/tau1
  B' = -B/tau2 
}


FUNCTION heav( x (1) ) {
    if ( x >= 0) {
        heav = 1
    } else {
        heav = 0
    }
}


:====== Netcon vector consists 3 variables ======: 
: w    intrinsic synaptic weight
: rho  plasticity factor (see original paper)
: c    calcium concentration
:====== ================================== ======: 
NET_RECEIVE(w (uS), rho, c ) {
    LOCAL timestep
    : YOU SHOULD SET UP ALL VECTOR VARIABLE IN HOC FILE!!!!!!
    INITIAL {  }
    if (flag == 0) {          : presynaptic spike 
        net_send(cdelay, 3)   : wait delay after presynapric event
        if (plastic) {
            A = A + factor*w*rho
            B = B + factor*w*rho
        } else {              : plasticity is turned off
            A = A + factor*w
            B = B + factor*w
        } 
    } else if (flag == 2) {   : postsynaptic spike
        net_send(0, 4) : update first
        FOR_NETCONS(wX, rhoX, cX) {
            cX = cX + cpost   : update all calcium concentration
        }
    } else if (flag == 3) {   : delay after presynaptic spike
        net_send(0, 4) : update first
        if ( learningdependence) {
            c = c + cpre*rho: update only local calcium concentration
        } else {
            c = c + cpre      : update only local calcium concentration
        }
    } else if (flag == 4) {   : it is time for update
        timestep = t-lndt
        lndt = t
        FOR_NETCONS(wX, rhoX, cX) {
            : Euler method for synaptic efficacy
            rhoX = rhoX +(-bistable*rhoX*(1.-rhoX)*(rho12-rhoX)+gammap*(1.-rhoX)*heav(cX-thetap)-gammad*rhoX*heav(cX-thetad))*timestep/taurho
:                         |         STDP bistability           |          Potentiation          |          Depression        |
            : Exponential Euler for calcium concentration
            cX = cX * exp(-timestep/tauca)
        }
        net_send(STDPstep, 4) : init next update

    } else {                  : flag == 1 from INITIAL block
        WATCH (v > -20) 2     : setup post-synaptic watcher...
    }
}

