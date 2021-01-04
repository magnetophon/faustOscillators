declare author "Bart Brouns";
declare license "GPLv3";
declare name "fof";
import("stdfaust.lib");

// process =
// fof;
// adapted from:
// https://ccrma.stanford.edu/~mjolsen/pdfs/smc2016_MOlsenFOF.pdf



/////////////////////////////////////////////////////////////////////////////
//                     Hard-Syncing Wavetable Oscillator                   //
/////////////////////////////////////////////////////////////////////////////
// resettable phasor, clock val > 0 resets phase to 0

ph(f0,c)=
  inc
  :(+:d)~
  (-(_<:(_,*(_,clk))))
  :*(pl.tablesize)

with {
  clk = c>0;
  d = ma.decimal;
  inc = f0/float(ma.SR);
};
// sin lookup table with resettable phase
oscpr(f0,c) =
  rdtable(pl.tablesize
         , os.sinwaveform(pl.tablesize)
         , int(ph(f0,c)));

// sinwaveform 	= float(time)*(2.0*PI)/float(tablesize) : sin;
// sinwaveform(tablesize) = float(ba.time)*(2.0*ma.PI)/float(tablesize) : sin;
///////////////////////////////////////////////////////////////////////////////
//                            FOF Generation System                          //
///////////////////////////////////////////////////////////////////////////////

// function to generate a single Formant-Wave-Function
fof(fc,bw,a,g) =
  _ <: (_',_)
  :
  (f
   * s
  )
with {
  T
  = 1/ma.SR;
  // sampling period
  pi = ma.PI;
  u1 = exp(-a*pi*T);
  u2 = exp(-bw*pi*T);
  a1 = -1*(u1+u2);
  a2 = u1*u2;
  G0 = 1/(1+a1+a2);
  // dc magnitude response
  b0 = g/G0;
  // normalized filter gain
  s
  = oscpr(fc);
  // wavetable oscillator
  f
  = fi.tf2(b0,0,0,a1,a2); // biquad filter
};

///////////////////////////////////////////////////////////////////////////////
//                        Cyclic Impulse Train Streams                       //
///////////////////////////////////////////////////////////////////////////////

// import the oscillator library
ol = library("oscillator.lib");
// impulse train at frequency f0
clk(f0) = (1-1')+ol.lf_imptrain(f0)';
// impulse train at frequency f0 split into n cycles
clkCycle(n,f0) = clk(f0) <: par(i,n,resetCtr(n,(i+1)));
// function that lets through the mth impulse out of
// each consecutive group of n impulses
resetCtr(n,m) = _ <: (_,ctr(n)) : (_,(_==m)) : *;
                                               // function to count nonzero inputs and reset after
                                               // receiving x of them
                                               ctr(x) = (+(_)~(negSub(x)));
// function that subtracts value x from
// input stream value if input stream value >= x
negSub(x)= _<: (_>=x,_,_):((-1*_),_,_):((_*_),_):(_+_);

///////////////////////////////////////////////////////////////////////////////
//                          FOF with Sample-and-Hold                         //
///////////////////////////////////////////////////////////////////////////////

// sample and hold filter coefficients
curbw = (_,bw) : ba.latch;
cura = (_,a) : ba.latch;
// FOF sample and hold mechanism
fofSH =
  _ <: (curbw,cura,_) : (fc,_,_,g,_') : fof ;

///////////////////////////////////////////////////////////////////////////////
//                              Example Program                              //
///////////////////////////////////////////////////////////////////////////////

/************** parameters/GUI controls **************/
// fundamental freq (in Hz)
f0 = hslider("f0",220,0,2000,0.01);
// formant center freq (in Hz)
fc = hslider("Fc",800,100,6000,0.01);
// FOF filter gain (in dB)
g = ba.db2linear(hslider("Gain",0,-40,40,0.01));
// FOF bandwidth (in Hz)
bw = hslider("BW",80,1,10000,1);
// FOF attack value (in Hz)
a = hslider("A",90,1,10000,1);
// number of S&H cycling filters
n = 5;
c = hslider("c", 0, 0, 1, 0.01);


// main process
process =
  // fof
  // rdtable(pl.tablesize, (ba.time:os.sinwaveform), int(ph(f0,c)))


  clkCycle(n,f0)
  <:
  par(i,n,fofSH)
  :> _
;
