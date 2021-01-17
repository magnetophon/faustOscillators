declare author "Bart Brouns";
declare license "GPLv3";
declare name "formantOscs";
import("stdfaust.lib");

process =
  // fofEnv(reset:ba.impulsify);
  // test =
  // PAF(center,fund(freq,phase,reset),Index,1)
  // fuphoSlaveEvenOdd(fund(freq,phase,reset),freq,1,Index,formantFreq,evenOdd)
  // ,
  // par(i, 50,
  fofosc(fund(freq,phase,reset))
  // nts3 (sin(2*ma.PI*fund(freq,0,0)),k1,k2,k3)*0.5
  // sigmoid
  // ,sigmoid
  // ):>_
  // fofEnv(rst:ba.impulsify)/2
  // lin_ar(attack,release,rst)
  // ,lin_ar(attack,release,rst)/2
  <:(_,_)
    // curved_ar
    // adsre(0.1,0.2,0.5,0.3,rst)
    // ade(0.1,0.2,rst)
;


// https://www.desmos.com/calculator/oufrdvzdcv
// based on:
// Adjustable Sigmoid Curve (S-Curve)
// https://math.stackexchange.com/questions/459872/adjustable-sigmoid-curve-s-curve-from-0-0-to-1-1
// https://www.desmos.com/calculator/tswgrnoosy

shaper(x,k1,k2,k3) = nts3(sin(ma.PI*(x-0.5)),k1,k2,k3);
nts3(x,k1,k2,k3) = fd(nts(fc(nts(fd(nts(fc(x),k1)),k2)),k3))
with {
  nts(x,k) = (x-x*k)/(k-abs(x)*2*k+1);
  fc(x) = x*0.5 + 0.5;
  fd(x) = 2*x-1;
};

// ntsSin(x,k1,k2,k3) =
// .99 = .956

// limiter shaper:
// https://www.desmos.com/calculator/hkqvmomfzp

k1 = hslider("k1", 0, -1, 1, 0.001):si.smoo;
k2 = hslider("k2", 0, -1, 1, 0.001):si.smoo;
k3 = hslider("k3", 0, -1, 1, 0.001):si.smoo;

// sigmoid(x) = 1.0 / (1.0 + exp(-x));
normalisedSigmoid(x,shape) = zeroSigmoid(x)/zeroSigmoid(1)
with {
  rawSigmoid(x) = 1.0 / (1.0 + exp(0-(x*2-1+n)*m));
  zeroSigmoid(x) = rawSigmoid(x)-rawSigmoid(0);
  n = l*0.8;
  m = (2+(abs(l)*4)):pow(2);
  l = shape;
};



// tanh shaper:
// https://www.desmos.com/calculator/ra0ekomhiv
// sigmoid shaper:
// https://www.desmos.com/calculator/p5lujdeecq
// improved sigmoid:
// https://www.desmos.com/calculator/nrgfzkhevw

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
//                     CZ osc for inside fof                               //
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
// f0 = hslider("f0",220,0,2000,0.01):si.smoo;
f0 = freq;
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
OLDfofosc(fund,octavation) =
  clkCycle(multi,fund)
  <: par(i,multi,fofSH*OctMuliply(i))
  :> _ * (octaviation+1) ;

fofosc(fund,octavation) =
  clkCycle(multi,fund)
  <: par(i,multi,myFOF(fund)*OctMuliply(i))
  :> _ * (octaviation:max(1)/maxOctavation) ;

myFOF(fund,fofTrig) =
  fofEnv(fund,fofTrig:ba.impulsify)*fofSlaveOsc(fofTrig);

fofEnv(fund,fofTrig) =
  // ade(attack,decay,fofTrig);
  curved_ar(fund,attack,release,ac,rc,fofTrig) ;
// attack = hslider("attack", 0, 0, 1, 0.001):pow(2);
// decay = hslider("decay", 0, 0, 1, 0.001):pow(2);
attack = hslider("attack", 2, 1, maxOctavation+1, 0.001):pow(2);
release = hslider("release", 2, 1, maxOctavation+1, 0.001):pow(2);
// release = hslider("release", 0, 0, 1, 0.001);
ac = hslider("ac", 0, -1, 1, 0.01);
rc = hslider("rc", 0, -1, 1, 0.01);

fofSlaveOsc(fofTrig) =
  fund(formantFreq,0,fofTrig)
  : basicCZosc(oscType,index,res);


lin_ar(a,r,gate) =
  FB~(_,_) // we need the previous state and value
     :!,_ // we only need the actual envelope output
with {
  FB(prevState, prevValue) =
    // states are :0 = release, 1= attack
    state
    // the actual value of the envelope
   ,envelope
  with { // so we can use prevState by name, without passing it as an argument
  // state = int(attackState * gate); // if no gate then force release
  state = (attackRamp < 1) | trigAttack;
  // attackState = prevState+trigAttack; // if no triggers, stay where you are, otherwise move to 1
  trigAttack= gate:ba.impulsify;
  attackRamp = ramp(attSamples,trigAttack);
  attSamples = int(a * ma.SR);

  LINenvelope = it.interpolate_linear(theRamp,from,to(state));
  LOGenvelope = it.interpolate_linear(fadeVal,from,to(state));
  SOFTenvelope = it.interpolate_linear(theRamp,keepDirection,to(state));
  keepDirection = ((_+oldDirection)~(_*(1-trig)))+from;
  oldDirection =
    ((prevValue-prevValue'):ba.sAndH(trig));
  envelope = LINenvelope;
  fadeVal = theRamp:LOGshaper;
  theRamp = ramp(samples,trig);
  samples = int(length(state) * ma.SR);
  length(state) = select2(state,r,a);
  trig = trigAttack, trigRelease :> _>0;
  trigRelease =
    (((prevState==1) & (attackRamp ==1)):ba.impulsify)  // doesn not auto-impulsify because prev will NOT increase
    | trigAttack & (attSamps==0);
  from = prevValue : ba.sAndH(trig);
  to(state) = select2(state,0,1); // release to 0, attack to 1, decay to sustain-level
  shaper(x) = it.interpolate_linear;
  LOGshaper(x) =
    // x;
    select2(shapeState
           ,map_log(x)
           ,x  // t == 1 should be linear, not 0
           );
  // possiby make 3rd state:  positive c mirorred
  map_log(x) = log (x * (t - 1) + 1) / log(t);  // from https://github.com/dariosanfilippo/edgeofchaos/blob/fff7e37ab80a5550421f7d4694c7de9b18b8b162/mathsEOC.lib#L881
  shapeState =
    t==1;
  // + 2*checkbox("shape"):min(3);
  // +t>1;
  t = pow((c/((c<0)+1))+1,16);
  c = curve(state);
  // adjust polarity so more is more
  curve(state) = select2(state,rc*-1,ac);
  attSamps = int(a * ma.SR);

  // ramp from 1/n to 1 in n samples.  (don't start at 0 cause when the ramp restarts, the crossfade should start right away)
  // when reset == 1, go back to 1/n.
  // ramp(n,reset) = select2(reset,_+(1/n):min(1),0)~_;
  ramp(n,reset) = (select2(reset,_+(1/n):min(1),1/n)~_)+(n<1):min(1);
};
};

curved_ar(fund,a,r,ac,rc,gate) =
  FB~(_,_) // we need the previous state and value
     :!,_ // we only need the actual envelope output
with {
  FB(prevState, prevValue) =
    // states are :0 = release, 1= attack
    state
    // the actual value of the envelope
   ,envelope
  with { // so we can use prevState by name, without passing it as an argument
  // state = int(attackState * gate); // if no gate then force release
  state = (attackRamp < 1) | trigAttack;
  // attackState = prevState+trigAttack; // if no triggers, stay where you are, otherwise move to 1
  trigAttack= gate:ba.impulsify;
  attackRamp = ramp(attSamples,trigAttack);
  // attSamples = int(a * ma.SR);
  attSamples = int( a * ma.SR / fund2freq(fund) );

  // TODO: make work with neg freq
  fund2freq(fund) = 1/delta
  with {
    absDelta = abs(fund-fund');
    delta =  select2(absDelta>0.5,absDelta,1-absDelta);
  };

  LINenvelope = it.interpolate_linear(theRamp,from,to(state));
  LOGenvelope = it.interpolate_linear(theRamp:LOGshaper,from,to(state));
  SIGenvelope = it.interpolate_linear(theRamp:SIGshaper,from,to(state));
  SOFTenvelope = it.interpolate_linear(theRamp:SIGshaper,keepDirection,to(state));
  keepDirection =
    from;
  // ((_+oldDirection)~(_*(1-trig)))+from;
  oldDirection =
    ((prevValue-prevValue'):ba.sAndH(trig));
  envelope = SOFTenvelope;
  fadeVal = theRamp:LOGshaper;
  theRamp = ramp(samples,trig);
  // samples = int(length(state) * ma.SR);
  // samples = int(length(state) * ma.SR);
  samples = int( length(state) * ma.SR / fund2freq(fund) );
  length(state) = select2(state,r,a);
  trig = trigAttack, trigRelease :> _>0;
  trigRelease =
    (((prevState==1) & (attackRamp ==1)):ba.impulsify)  // doesn not auto-impulsify because prev will NOT increase
    | trigAttack & (attSamps==0);
  from = prevValue : ba.sAndH(trig);
  to(state) = select2(state,0,1); // release to 0, attack to 1, decay to sustain-level
  // shaper(x) = it.interpolate_linear;
  LOGshaper(x) =
    select2(shapeState
           ,map_log(x)
           ,x  // t == 1 should be linear, not 0
           );
  // possiby make 3rd state:  positive c mirorred
  map_log(x) = log (x * (t - 1) + 1) / log(t);  // from https://github.com/dariosanfilippo/edgeofchaos/blob/fff7e37ab80a5550421f7d4694c7de9b18b8b162/mathsEOC.lib#L881
  shapeState =
    t==1;
  // + 2*checkbox("shape"):min(3);
  // +t>1;
  t = pow((c/((c<0)+1))+1,16);
  c = curve(state);
  // adjust polarity so more is more
  curve(state) = select2(state,rc*-1,ac);
  SIGshaper(x) = normalisedSigmoid(x,curve(state));
  // SIGshaper(x) = normalisedSigmoid(x,0);
  attSamps = int(a * ma.SR);

  // ramp from 1/n to 1 in n samples.  (don't start at 0 cause when the ramp restarts, the crossfade should start right away)
  // when reset == 1, go back to 1/n.
  // ramp(n,reset) = select2(reset,_+(1/n):min(1),0)~_;
  ramp(n,reset) = (select2(reset,_+(1/n):min(1),1/n)~_)+(n<1):min(1);
};
};

adsre(attT60,decT60,susLvl,relT60,gate) = envelope with {
  ugate = gate>0;
  samps = ugate : +~(*(ugate)); // ramp time in samples
  attSamps = int(attT60 * ma.SR);
  // the first sample of each note is alwaus the attack phase, also when attSamps==0
  attPhase = (samps<attSamps) |  (ugate:ba.impulsify);
  // attPhase = (samps<attSamps) | ((attSamps==0) & (ugate:ba.impulsify));
  target = select2(ugate, 0.0,
                   select2(attPhase, (susLvl)*float(ugate), ugate));
  t60 = select2(ugate, relT60, select2(attPhase, decT60, attT60));
  pole = ba.tau2pole(t60/6.91);
  envelope = target : si.smooth(pole);
};

ade(attT60,decT60,gate) = envelope with {
  ugate = gate>0;
  samps = ugate : +~(*(ugate)); // ramp time in samples
  attSamps = int(attT60 * ma.SR);
  // the first sample of each note is alwaus the attack phase, also when attSamps==0
  attPhase = (samps<attSamps) |  (ugate:ba.impulsify);
  // attPhase = (samps<attSamps) | ((attSamps==0) & (ugate:ba.impulsify));
  target = attPhase*ugate;
  t60 = select2(attPhase, decT60, attT60);
  pole = ba.tau2pole(t60/6.91);
  envelope = target : si.smooth(pole);
};

OctMuliply(i) =
  (
    (OctMuliplyPart(i,floor(octaviation+1))*    ma.decimal(octaviation)) +
    (OctMuliplyPart(i,floor(octaviation  ))* (1-ma.decimal(octaviation)))
  )
  // :hbargraph("octmult %2i", 0, 1)
;
OctMuliplyPart(i,oct) =
  select2(i>0,1,
          (((i % int(2:pow(oct)))):min(1)*-1)+1) ;

multi = 2:pow(maxOctavation);
maxOctavation = 4;
octaviation = hslider("octaviation",0,0,maxOctavation,0.001):si.smoo;
