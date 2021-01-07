declare author "Bart Brouns";
declare license "GPLv3";
declare name "formantOscs";
import("stdfaust.lib");

process =
  // PAF(center,fund(freq,phase,reset),Index,1)
  // ,fuphoSlaveEvenOdd(fund(freq,phase,reset),freq,1,Index,formantFreq,evenOdd)
  // ,
  // par(i, 50,
  fofosc(fund(freq,phase,reset))
  // ):>_
  // <:(_,_)
;

FMindex = hslider("[2]FM index[tooltip: FM index of the lowest oscillator]",	1, 0, 20, 0.001):max(0.00000001);
evenOdd = hslider("evenOdd", 0.5, 0, 1, 0.001):si.smoo;

formantFreq = hslider("formant frequency", 440, 27.5, 3520, 0.1):si.smoo;
Index = hslider("[0]index[tooltip: PAF index]",	100, 0.001, 100, 0.001);
// phase = hslider("phase", 0, 0, 1, 0.001):si.smoo;
//-----------------------------------------------
// PAF oscilator
//-----------------------------------------------
// http://msp.ucsd.edu/techniques/v0.11/book-html/node96.html
//
// center = vslider("[01]bottom[style:knob][tooltip: the lowest formant frequency of the oscillators]",	1, 0.5, 7, 0.001):pow(2):si.smoo*freq;			//0.25 to 49 logarithmicly
center = formantFreq / freq;
// freq = hslider("freq", 440, 27.5, 3520, 0.1):si.smoo;
pafIndex = hslider("[0]index[tooltip: PAF index]",	100, 0.001, 100, 0.001);
// deglitch adapted from Chris Chafe:
// http://chrischafe.net/glitch-free-fm-vocal-synthesis/
// center, fund, index, volume
PAF(c,f,i,v) =
  BasicPAF(ef,f,i,ea)+
  BasicPAF(of,f,i,oa)
with {
  ci     = int(floor(c));          // integer harmonic below center freq
  ci1    = ci+1;                   // and above
  isEven = ba.if(((ci%2)<1),1,0);     // below is even?
  ef     = ba.if(isEven,ci,ci1);      // then set even harmonic to lowest
  of     = ba.if(isEven,ci1,ci);      // and odd harmonic to highest
  frac   = c-ci;                   // fractional frequency remainder
  comp   = 1-frac;                 // and it's complement
  oa     = ba.if(isEven,frac,comp)*v; // odd harmonic amplitude
  ea     = ba.if(isEven,comp,frac)*v; // even harmonic amplitude
};

sampleAndHold(sample) = select2((sample!=0):int) ~ _;

// center, fund, index, volume
BasicPAF(c,f,i,v)=
  (((cos12(c,f))*bell(f,i)) * v) * volumeCompensate
with {
  bellcurve(x) = int(x):rdtable(belltablesize + 1,curve,_)
  with 	{
  belltablesize 	= 200;
  curve 	= ba.time:(((_-100)/25)<:exp((_*_)*-1));
};
  volumeCompensate = 0.666+(i/300);
  // volumeCompensate = 1;
  wrap            = _<:(_>0,(_,1:fmod)+1,(_,1:fmod)):select2;
  centerWrap(c,f) = c:sampleAndHold(f)<:wrap;
  centerMin(c,f)  = c:sampleAndHold(f)-centerWrap(c,f);
  cos12(c,f)      = (centerMin(c,f)*f<:(_*2*ma.PI:cos)<:(_,_((_,(_+f:(_*2*ma.PI:cos))):_-_:(_*centerWrap(c,f))) )):_+_;
  bell(f,i)       = (((f*0.5)-0.25:(_*2*ma.PI:cos))*i)+100:bellcurve;
};
//-----------------------------------------------
// FM formant oscilator
//-----------------------------------------------

// declare name 		"FMVox";
// declare version		"1.0";
// declare author		"Chris Chafe";
// declare license		"BSD";
// declare copyright	"stk";

// import("stdfaust.lib");

// ------ Program A of "Glitch Free FM Vocal Synthesis" (minor corrections from article)

// ------ table lookup oscillators -----//
ts 	= 1 << 16;			// table size
fs 	= float(ts);
ts1	= ts+1;
ct 	= (+(1)~_ ) - 1;		// incrementing counter from 0
fct	= float(ct);
sc 	= fct*(2*ma.PI)/fs:sin;		// sine table for carrier
sm	= fct*(2*ma.PI)/fs:sin:/(2*ma.PI);	// sine table for modulator
dec(x)	= x-floor(x);			// fractional remainder
pha(f)	= f/float(ma.SR):(+:dec) ~ _;	// generates a index signal
tbl(t,p)= s1
          // tbl(t,p)= s1+dec(f)*(s2-s1)		// looks up linearly interpolated table value
with {
  f = p*fs;
  i = int(f):min(ts);
  s1 = rdtable(ts1,t,i);
  s2 = rdtable(ts1,t,i+1);
};

///////////////////////////////////////////////////////////////////////////////
//      formant generator using uniform (phase-synchronous) oscillators      //
///////////////////////////////////////////////////////////////////////////////

// synced to an external phaser "index"
// the parameter eo sets the proportion of even and odd harmonics.
fuphoSlaveEvenOdd(index,f0,a,b,c,eo) = (even+odd):*(a)	// outputs the sum of bracketing harmonics
with {					// from f0, amp, bandwidth, center freq
  okt = vslider("okt", 0, -2, 2, 1);
  cf 	= c/f0;
  // cf 	= c/(f0*(okt:octaveMultiplier));
  ci	= int(floor(cf));			// integer harmonic below center freq
  ci1	= ci+1;				// and above
  isEven= ba.if(((ci%2)<1),1,0);	// below is even?
  ef 	= ba.if(isEven,ci,ci1);		// then set even harmonic to lowest
  of 	= ba.if(isEven,ci1,ci);		// and odd harmonic to highest
  frac	= cf-ci;			// fractional frequency remainder
  comp	= 1-frac;
  oa 	= ba.if(isEven,frac,comp)*eo;		// odd harmonic amplitude
  ea 	= ba.if(isEven,comp,frac)*((eo*-1)+1);		// even harmonic amplitude
  // ph 	= pha(f0);			// index signal at fundamental
  m 	= tbl(sm,index):*(b);		// modulator sine signal
  even 	= ea:*(tbl(sc,(dec(ef*index+m)))); // even harmonic signal with phase modulation
  odd 	= oa:*(tbl(sc,(dec(of*index+m))));	// odd harmonic signal
};

///////////////////////////////////////////////////////////////////////////////
//                                    fof                                    //
///////////////////////////////////////////////////////////////////////////////

// declare author "Bart Brouns";
// declare license "GPLv3";
// declare name "fof";
import("CZ.lib");

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

/////////////////////////////////////////////////////////////////////////////
//                     CZ osc                   //
/////////////////////////////////////////////////////////////////////////////

oscpr2(f0,c) =
  fund(f0,0,c)
  : basicCZosc(oscType,index,res);

///////////////////////////////////////////////////////////////////////////////
//                            FOF Generation System                          //
///////////////////////////////////////////////////////////////////////////////

// function to generate a single Formant-Wave-Function
fof(fc,bw,a,g) =
  _ <: (_',_)
  : (oscil*env)
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
  oscil
  = oscpr2(fc);
  // wavetable oscillator
  env
  = fi.tf2(b0,0,0,a1,a2); // biquad filter
};

///////////////////////////////////////////////////////////////////////////////
//                        Cyclic Impulse Train Streams                       //
///////////////////////////////////////////////////////////////////////////////

// impulse train at frequency of fund
clk(fund) = (1-1')+(fund<:-(mem)<0)';
// impulse train at frequency f0 split into n cycles
clkCycle(n,fund) = clk(fund) <: par(i,n,resetCtr(n,(i+1)));
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
f0 = hslider("f0",220,0,2000,0.01):si.smoo;
// formant center freq (in Hz)
fc = hslider("Fc",800,100,6000,0.01):si.smoo;
// FOF filter gain (in dB)
g = ba.db2linear(hslider("Gain",0,-40,40,0.01)):si.smoo;
// FOF bandwidth (in Hz)
bw = hslider("BW",80,1,10000,1):si.smoo;
// FOF attack value (in Hz)
a = hslider("A",90,1,10000,1):si.smoo;
rst = button("reset");


// main process
fofosc(fund,octavation) =
  clkCycle(multi,fund)
  <: par(i,multi,fofSH*OctMuliply(i))
  :> _ * (octaviation+1) ;

OctMuliply(i) =
  (OctMuliplyPart(i,floor(octaviation+1))*    ma.decimal(octaviation)) +
  (OctMuliplyPart(i,floor(octaviation  ))* (1-ma.decimal(octaviation)));
OctMuliplyPart(i,oct) =
  select2(i>0,1,
          (((i % int(2:pow(oct)))):min(1)*-1)+1);

multi = 2:pow(maxOctavation);
maxOctavation = 4;
octaviation = hslider("octaviation",0,0,maxOctavation,0.001):si.smoo;
