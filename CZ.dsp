declare author "Bart Brouns";
declare license "GPLv3";
declare name "CZoscs";
import("stdfaust.lib");

process =
  full<:(_,_);
// basic ,medium full;


// "basic" is an almost normal CZ oscillator.
// I did a PR for it to faustlibraries, but for convenience, I pasted the code below,
// from the line marked with "Casio CZ Oscillators"
// It is documented more thoroughly there, but basically, it has the following changes:
// - it is phase alligned to a master oscillator called "fund"
// - it's anti-aliased by decreasing the index at higher frequencies.
basic =
  fund(freq,phase,reset)
  : basicCZosc(oscType,index,res);

// Like basic, but  adds a filter between the master oscillator and the CZ oscillator.
medium =
  fund(freq,phase,reset)
  : prefilterCZosc(oscType,filterType,filterFreq,filterQ,index,res);

// Like medium, but can crossfade between octaves
full =
  fund(freq*minOctMult,phase,reset)
  : prefilterOctaveCZosc(oscType,filterType,filterFreq,filterQ,index,res,oct) ;

basicCZosc(oscType,index,res,fund) =
  (fund,index,res)
  : oscillatorChooser(oscType);

prefilterCZosc(oscType,filterType,filterFreq,filterQ,index,res,fund) =
  fund
  : filterChooser(filterType,filterFreq,filterQ)
  : basicCZosc(oscType,index,res);

prefilterOctaveCZosc(oscType,filterType,filterFreq,filterQ,index,res,oct,fund) =
  fund
  :oscOctaver
   ( oct
   , prefilterCZosc(oscType,filterType,filterFreq,filterQ,index,res));

fund(freq,phase,reset) = lf_sawpos_phase_reset(freq,phase,reset);

oscillatorChooser(oscType,fund,index,res) =
  chooseOne(allOscsParallel(fund,index,res),oscType);

filterChooser(filterType,filterFreq,filterQ) =
  _<:chooseOne(allFiltersParallel(filterFreq,filterQ),filterType);

chooseOne(proc,N) = // choose the Nth processor out of a list of processors with one output each
  proc
  // : par(i, maxN, control(i==clampedN)*(i==clampedN)):>_ // only actually run one of the procs
  : par(i, maxN, _ *(i==clampedN)):>_ // run them all, but mute all but one
with {
  clampedN = int(min(N,maxN-1));
  maxN = outputs(proc);
};

allFiltersParallel(f,q) =
  _
, svf.lp(f,q)
, svf.bp(f,q)
, svf.hp(f,q)
, svf.notch(f,q)
, svf.peak(f,q)
, svf.ap(f,q);

allOscsParallel(fund,index,res) =
  sineNoise(fund,index)
, CZsawPAA(fund, index)
, CZsquarePAA(fund, index)
, CZpulsePAA(fund, index)
, CZsinePulsePAA(fund, index)
, CZhalfSinePAA(fund, index)
, CZresSawAA(fund,res)
, CZresTriangleAA(fund,res)
, CZresTrapAA(fund, res);

oscOctaver(oct,oscil,fund) =
  (
    (f0:oscil)
   ,(f1:oscil)
  ):si.interpolate(oct:abs%1)
with {
  f0 = (fund+phaseM(oct0):ma.frac):octaveSwitcher(oct0);
  f1 = (fund+phaseM(oct1):ma.frac):octaveSwitcher(oct1);
  oct0 = oct:floor+((oct<0) & (oct!=(oct:floor)));
  oct1 = oct:floor+(oct>0);
  octaveSwitcher(oct) = _*(octaveMultiplier(oct)/minOctMult)%1;
  phaseM(oct) = phase*octMult(oct);
  octMult(oct)=
    octaveMultiplier(minOct)/(1/pow(2,oct*-1));
};

octaveMultiplier =
  int<:
  (
    (_ <0) / pow(2,abs),
    (_==0),
    (_ >0) * pow(2,_)
  ):>_;

minOctMult = minOct:octaveMultiplier;

// When index is 0, you hear a sine, just like with all the other oscillators.
// When index is between 0 and 0.5, it smoothly morphs into band-pass noise at the freq of the oscillator.
// When index is between 0.5 and 1, it crossfades into low-pass noise.
sineNoise(fund,index) =
  it.interpolate_linear
  ( (index*2:max(1)-1)
  , (fund*2*ma.PI:sin*(it.interpolate_linear(index*2:min(1),1,lfnoiseF(fund))))
  , lfnoiseF(fund)
  )
with {
  index2freq(i)        = ((i-i')*ma.SR) : ba.sAndH(abs(i-i')<0.5):int :max(20):min(ma.SR*.5);
  NF = index2freq(fund) :int :max(20):min(ma.SR*.25);
  lfnoise0F(fund) = no.noise : ba.latch((fund:ma.decimal)<0.5);
  lfnoiseNF(N,fund) = lfnoise0F(fund) : fi.lowpass(N,index2freq(fund) ); // Nth-order Butterworth lowpass
  lfnoiseF(fund) = lfnoise0F(fund) : seq(i,5,fi.lowpass(1,index2freq(fund))); // non-overshooting lowpass
};

// params:
// osc type
// filter type
// filter freq
// filter q
// osc index
// osc res
// osc octave
// osc phase
//
//
// TODO:
// fix click in oct
//PR sineNoise, credit/ask?
// include subSynth?
// make brightness vs index correction, maybe PR


///////////////////////////////////////////////////////////////////////////////
//                                    GUI                                    //
///////////////////////////////////////////////////////////////////////////////

// type = hslider("type[scale:int][style:menu{'Sine-Noise':0;'Sawtooth':1;'Square':2;'Pulse':3;'Sine-Pulse':4;'Half Sine':5;'Resonant Saw':6;'Resonant Tri':7;'Resonant Trap':8}]", 2, 0, 8, 1);
oscType = hslider("[0] oscillator type[scale:int][style:menu{'Sine-Noise':0;'Sawtooth':1;'Square':2;'Pulse':3;'Sine-Pulse':4;'Half Sine':5;'Resonant Saw':6;'Resonant Tri':7;'Resonant Trap':8}]", 2, 0, 8, 1);
filterType = hslider("[1] filter type[scale:int]
[style:menu
    {
    'none':0;
    'lp':1;
    'bp':2;
    'hp':3;
	'notch':4;
	'peak':5;
	'ap':6
    }
]", 0, 0, 6, 1);

filterFreq = hslider("[2] filter freq[scale:log]", 24000, 20, 24000, 1) :si.smoo;
filterQ = hslider("[3] filter Q", 1, 0.001, 10, 0.001) :si.smoo;
// an osc has either index or res, never both, so using the same [number] is OK
index = hslider("[4] index", 0, 0, 1, stepsize) :si.smoo;
res   = hslider("[4] res", 0, 0, 64, stepsize) :si.smoo;
oct   = hslider("[5] octave", 0, minOct, maxOct, stepsize) :si.smoo;
phase = hslider("[6] phase", 0, -64, 64, stepsize) :si.smoo;
freq = hslider("[7]freq", 440, 20, 24000, 1) :si.smoo;
reset = button("[8]reset oscillator");

///////////////////////////////////////////////////////////////////////////////
//                                 constants                                 //
///////////////////////////////////////////////////////////////////////////////

minOct = -8;
maxOct = 4;

// fast
// stepsize = 0.1;
// medium
stepsize = 0.01;
// smooth
// stepsize = 0.001;

///////////////////////////////////////////////////////////////////////////////
//           not yet in my version of faust, but merged in master:           //
///////////////////////////////////////////////////////////////////////////////

lf_sawpos_phase_reset(freq,phase,reset) = lf_sawpos_reset(freq,reset) + phase : ma.frac;
lf_sawpos_reset(freq,reset) = ma.frac * (reset == 0) ~ +(freq/ma.SR);


declare svf author "Oleg Nesterov";
declare svf copyright "Copyright (C) 2020 Oleg Nesterov <oleg@redhat.com>";
declare svf license "MIT-style STK-4.3 license";

svf = environment {
	    svf(T,F,Q,G) = tick ~ (_,_) : !,!,si.dot(3, mix)
	    with {
		tick(ic1eq, ic2eq, v0) =
		  2*v1 - ic1eq,
		  2*v2 - ic2eq,
		  v0, v1, v2
		with {
		v1 = ic1eq + g *(v0-ic2eq) : /(1 + g*(g+k));
		v2 = ic2eq + g * v1;
	  };

		A = pow(10.0, G/40.0);

		g = tan(F * ma.PI/ma.SR) : case {
			  (7) => /(sqrt(A));
			  (8) => *(sqrt(A));
			  (t) => _;
		    } (T);

		k = case {
			  (6) => 1/(Q*A);
			  (t) => 1/Q;
		    } (T);

		mix = case {
			    (0) => 0, 0, 1;
			    (1) => 0, 1, 0;
			    (2) => 1, -k, -1;
			    (3) => 1, -k, 0;
			    (4) => 1, -k, -2;
			    (5) => 1, -2*k, 0;
			    (6) => 1, k*(A*A-1), 0;
			    (7) => 1, k*(A-1), A*A-1;
			    (8) => A*A, k*(1-A)*A, 1-A*A;
		      } (T);
	  };

	    // External API
	    lp(f,q)     = svf(0, f, q, 0);
	    bp(f,q)     = svf(1, f, q, 0);
	    hp(f,q)     = svf(2, f, q, 0);
	    notch(f,q)  = svf(3, f, q, 0);
	    peak(f,q)   = svf(4, f, q, 0);
	    ap(f,q)     = svf(5, f, q, 0);
	    bell(f,q,g) = svf(6, f, q, g);
	    ls(f,q,g)   = svf(7, f, q, g);
	    hs(f,q,g)   = svf(8, f, q, g);
      };

///////////////////////////////////////////////////////////////////////////////
//            https://github.com/grame-cncm/faustlibraries/pull/64           //
//            CZ Oscillators: add anti-aliased versions                      //
///////////////////////////////////////////////////////////////////////////////



//===================== Casio CZ Oscillators ==========================
// Oscillators that mimic some of the Casio CZ oscillators.
//
// There are two sets:
// - A set with an index parameter
// - A set with a res parameter
//
// The "index oscillators" output a sine wave at index=0 and gets brighter with a higher index.
// There are two versions of the "index oscillators":
// - with P appended to the name: is phase aligned with 'fund:sin'
// - without P appended to the name: has the phase of the original CZ oscillators
//
// The "res oscillators" have a resonant frequency.
// "res" is the frequency of resonance as a factor of the fundamental pitch.
//
// All oscillators also have an anti-aliased version, marked with AA
//=====================================================================

//----------`(os.)CZsaw`----------
// Oscillator that mimics the Casio CZ saw oscillator
// `CZsaw` is a standard Faust function.
//
// #### Usage
//
// ```
// CZsaw(fund,index) : _
// ```
//
// Where:
//
// * `fund`: a saw-tooth waveform between 0 and 1 that the oscillator slaves to
// * `index`: the brightness of the oscillator, 0 to 1. 0 = sine-wave, 1 = saw-wave
//------------------------------------------------------------
// Author: Bart Brouns
// License: GPLv3
// CZ oscillators by Mike Moser-Booth:
// <https://forum.pdpatchrepo.info/topic/5992/casio-cz-oscillators>
// Ported from pd to Faust by Bart Brouns

CZsaw(fund, index) = CZ.saw(fund, index);

//----------`(os.)CZsawAA`----------
// Oscillator that mimics the Casio CZ saw oscillator,
// anti-aliased by decreasing the index at higher frequencies.
// `CZsawAA` is a standard Faust function.
//
// #### Usage
//
// ```
// CZsawAA(fund,index) : _
// ```
//
// Where:
//
// * `fund`: a saw-tooth waveform between 0 and 1 that the oscillator slaves to
// * `index`: the brightness of the oscillator, 0 to 1. 0 = sine-wave, 1 = saw-wave
//------------------------------------------------------------
// Author: Bart Brouns
// License: GPLv3
// CZ oscillators by Mike Moser-Booth:
// <https://forum.pdpatchrepo.info/topic/5992/casio-cz-oscillators>
// Ported from pd to Faust by Bart Brouns

CZsawAA(fund, index) = CZ.sawAA(fund, index);

//----------`(os.)CZsawP`----------
// Oscillator that mimics the Casio CZ saw oscillator,
// with it's phase aligned to `fund:sin`.
// `CZsawP` is a standard Faust function.
//
// #### Usage
//
// ```
// CZsawP(fund,index) : _
// ```
//
// Where:
//
// * `fund`: a saw-tooth waveform between 0 and 1 that the oscillator slaves to
// * `index`: the brightness of the oscillator, 0 to 1. 0 = sine-wave, 1 = saw-wave
//------------------------------------------------------------
// Author: Bart Brouns
// License: GPLv3
// CZ oscillators by Mike Moser-Booth:
// <https://forum.pdpatchrepo.info/topic/5992/casio-cz-oscillators>
// Ported from pd to Faust by Bart Brouns

CZsawP(fund, index) = CZ.sawP(fund, index);

//----------`(os.)CZsawPAA`----------
// Oscillator that mimics the Casio CZ saw oscillator,
// with it's phase aligned to `fund:sin`,
// anti-aliased by decreasing the index at higher frequencies.
// `CZsawPAA` is a standard Faust function.
//
// #### Usage
//
// ```
// CZsawPAA(fund,index) : _
// ```
//
// Where:
//
// * `fund`: a saw-tooth waveform between 0 and 1 that the oscillator slaves to
// * `index`: the brightness of the oscillator, 0 to 1. 0 = sine-wave, 1 = saw-wave
//------------------------------------------------------------
// Author: Bart Brouns
// License: GPLv3
// CZ oscillators by Mike Moser-Booth:
// <https://forum.pdpatchrepo.info/topic/5992/casio-cz-oscillators>
// Ported from pd to Faust by Bart Brouns

CZsawPAA(fund, index) = CZ.sawPAA(fund, index);

//----------`(os.)CZsquare`----------
// Oscillator that mimics the Casio CZ square oscillator
// `CZsquare` is a standard Faust function.
//
// #### Usage
//
// ```
// CZsquare(fund,index) : _
// ```
//
// Where:
//
// * `fund`: a saw-tooth waveform between 0 and 1 that the oscillator slaves to
// * `index`: the brightness of the oscillator, 0 to 1. 0 = sine-wave, 1 = square-wave
//------------------------------------------------------------
// Author: Bart Brouns
// License: GPLv3
// CZ oscillators by Mike Moser-Booth:
// <https://forum.pdpatchrepo.info/topic/5992/casio-cz-oscillators>
// Ported from pd to Faust by Bart Brouns

CZsquare(fund, index) = CZ.square(fund, index);

//----------`(os.)CZsquareAA`----------
// Oscillator that mimics the Casio CZ square oscillator,
// anti-aliased by decreasing the index at higher frequencies.
// `CZsquareAA` is a standard Faust function.
//
// #### Usage
//
// ```
// CZsquareAA(fund,index) : _
// ```
//
// Where:
//
// * `fund`: a saw-tooth waveform between 0 and 1 that the oscillator slaves to
// * `index`: the brightness of the oscillator, 0 to 1. 0 = sine-wave, 1 = square-wave
//------------------------------------------------------------
// Author: Bart Brouns
// License: GPLv3
// CZ oscillators by Mike Moser-Booth:
// <https://forum.pdpatchrepo.info/topic/5992/casio-cz-oscillators>
// Ported from pd to Faust by Bart Brouns

CZsquareAA(fund, index) = CZ.squareAA(fund, index);

//----------`(os.)CZsquareP`----------
// Oscillator that mimics the Casio CZ square oscillator,
// with it's phase aligned to `fund:sin`.
// `CZsquareP` is a standard Faust function.
//
// #### Usage
//
// ```
// CZsquareP(fund,index) : _
// ```
//
// Where:
//
// * `fund`: a saw-tooth waveform between 0 and 1 that the oscillator slaves to
// * `index`: the brightness of the oscillator, 0 to 1. 0 = sine-wave, 1 = square-wave
//------------------------------------------------------------
// Author: Bart Brouns
// License: GPLv3
// CZ oscillators by Mike Moser-Booth:
// <https://forum.pdpatchrepo.info/topic/5992/casio-cz-oscillators>
// Ported from pd to Faust by Bart Brouns

CZsquareP(fund, index) = CZ.squareP(fund, index);

//----------`(os.)CZsquarePAA`----------
// Oscillator that mimics the Casio CZ square oscillator,
// with it's phase aligned to `fund:sin`,
// anti-aliased by decreasing the index at higher frequencies.
// `CZsquarePAA` is a standard Faust function.
//
// #### Usage
//
// ```
// CZsquarePAA(fund,index) : _
// ```
//
// Where:
//
// * `fund`: a saw-tooth waveform between 0 and 1 that the oscillator slaves to
// * `index`: the brightness of the oscillator, 0 to 1. 0 = sine-wave, 1 = square-wave
//------------------------------------------------------------
// Author: Bart Brouns
// License: GPLv3
// CZ oscillators by Mike Moser-Booth:
// <https://forum.pdpatchrepo.info/topic/5992/casio-cz-oscillators>
// Ported from pd to Faust by Bart Brouns

CZsquarePAA(fund, index) = CZ.squarePAA(fund, index);

//----------`(os.)CZpulse`----------
// Oscillator that mimics the Casio CZ pulse oscillator
// `CZpulse` is a standard Faust function.
//
// #### Usage
//
// ```
// CZpulse(fund,index) : _
// ```
//
// Where:
//
// * `fund`: a saw-tooth waveform between 0 and 1 that the oscillator slaves to
// * `index`: the brightness of the oscillator, 0 gives a sine-wave, 1 is closer to a pulse
//------------------------------------------------------------
// Author: Bart Brouns
// License: GPLv3
// CZ oscillators by Mike Moser-Booth:
// <https://forum.pdpatchrepo.info/topic/5992/casio-cz-oscillators>
// Ported from pd to Faust by Bart Brouns

CZpulse(fund, index) = CZ.pulse(fund, index);

//----------`(os.)CZpulseAA`----------
// Oscillator that mimics the Casio CZ pulse oscillator,
// anti-aliased by decreasing the index at higher frequencies.
// `CZpulseAA` is a standard Faust function.
//
// #### Usage
//
// ```
// CZpulseAA(fund,index) : _
// ```
//
// Where:
//
// * `fund`: a saw-tooth waveform between 0 and 1 that the oscillator slaves to
// * `index`: the brightness of the oscillator, 0 gives a sine-wave, 1 is closer to a pulse
//------------------------------------------------------------
// Author: Bart Brouns
// License: GPLv3
// CZ oscillators by Mike Moser-Booth:
// <https://forum.pdpatchrepo.info/topic/5992/casio-cz-oscillators>
// Ported from pd to Faust by Bart Brouns

CZpulseAA(fund, index) = CZ.pulseAA(fund, index);

//----------`(os.)CZpulseP`----------
// Oscillator that mimics the Casio CZ pulse oscillator,
// with it's phase aligned to `fund:sin`.
// `CZpulseP` is a standard Faust function.
//
// #### Usage
//
// ```
// CZpulseP(fund,index) : _
// ```
//
// Where:
//
// * `fund`: a saw-tooth waveform between 0 and 1 that the oscillator slaves to
// * `index`: the brightness of the oscillator, 0 gives a sine-wave, 1 is closer to a pulse
//------------------------------------------------------------
// Author: Bart Brouns
// License: GPLv3
// CZ oscillators by Mike Moser-Booth:
// <https://forum.pdpatchrepo.info/topic/5992/casio-cz-oscillators>
// Ported from pd to Faust by Bart Brouns

CZpulseP(fund, index) = CZ.pulseP(fund, index);

//----------`(os.)CZpulsePAA`----------
// Oscillator that mimics the Casio CZ pulse oscillator,
// with it's phase aligned to `fund:sin`,
// anti-aliased by decreasing the index at higher frequencies.
// `CZpulsePAA` is a standard Faust function.
//
// #### Usage
//
// ```
// CZpulsePAA(fund,index) : _
// ```
//
// Where:
//
// * `fund`: a saw-tooth waveform between 0 and 1 that the oscillator slaves to
// * `index`: the brightness of the oscillator, 0 gives a sine-wave, 1 is closer to a pulse
//------------------------------------------------------------
// Author: Bart Brouns
// License: GPLv3
// CZ oscillators by Mike Moser-Booth:
// <https://forum.pdpatchrepo.info/topic/5992/casio-cz-oscillators>
// Ported from pd to Faust by Bart Brouns

CZpulsePAA(fund, index) = CZ.pulsePAA(fund, index);

//----------`(os.)CZsinePulse`----------
// Oscillator that mimics the Casio CZ sine/pulse oscillator
// `CZsinePulse` is a standard Faust function.
//
// #### Usage
//
// ```
// CZsinePulse(fund,index) : _
// ```
//
// Where:
//
// * `fund`: a saw-tooth waveform between 0 and 1 that the oscillator slaves to
// * `index`: the brightness of the oscillator, 0 gives a sine-wave, 1 is a sine minus a pulse
//------------------------------------------------------------
// Author: Bart Brouns
// License: GPLv3
// CZ oscillators by Mike Moser-Booth:
// <https://forum.pdpatchrepo.info/topic/5992/casio-cz-oscillators>
// Ported from pd to Faust by Bart Brouns

CZsinePulse(fund, index) = CZ.sinePulse(fund, index);

//----------`(os.)CZsinePulseAA`----------
// Oscillator that mimics the Casio CZ sine/pulse oscillator,
// anti-aliased by decreasing the index at higher frequencies.
// `CZsinePulseAA` is a standard Faust function.
//
// #### Usage
//
// ```
// CZsinePulseAA(fund,index) : _
// ```
//
// Where:
//
// * `fund`: a saw-tooth waveform between 0 and 1 that the oscillator slaves to
// * `index`: the brightness of the oscillator, 0 gives a sine-wave, 1 is a sine minus a pulse
//------------------------------------------------------------
// Author: Bart Brouns
// License: GPLv3
// CZ oscillators by Mike Moser-Booth:
// <https://forum.pdpatchrepo.info/topic/5992/casio-cz-oscillators>
// Ported from pd to Faust by Bart Brouns

CZsinePulseAA(fund, index) = CZ.sinePulseAA(fund, index);

//----------`(os.)CZsinePulseP`----------
// Oscillator that mimics the Casio CZ sine/pulse oscillator,
// with it's phase aligned to `fund:sin`.
// `CZsinePulseP` is a standard Faust function.
//
// #### Usage
//
// ```
// CZsinePulseP(fund,index) : _
// ```
//
// Where:
//
// * `fund`: a saw-tooth waveform between 0 and 1 that the oscillator slaves to
// * `index`: the brightness of the oscillator, 0 gives a sine-wave, 1 is a sine minus a pulse
//------------------------------------------------------------
// Author: Bart Brouns
// License: GPLv3
// CZ oscillators by Mike Moser-Booth:
// <https://forum.pdpatchrepo.info/topic/5992/casio-cz-oscillators>
// Ported from pd to Faust by Bart Brouns

CZsinePulseP(fund, index) = CZ.sinePulseP(fund, index);

//----------`(os.)CZsinePulsePAA`----------
// Oscillator that mimics the Casio CZ sine/pulse oscillator,
// with it's phase aligned to `fund:sin`,
// anti-aliased by decreasing the index at higher frequencies.
// `CZsinePulsePAA` is a standard Faust function.
//
// #### Usage
//
// ```
// CZsinePulsePAA(fund,index) : _
// ```
//
// Where:
//
// * `fund`: a saw-tooth waveform between 0 and 1 that the oscillator slaves to
// * `index`: the brightness of the oscillator, 0 gives a sine-wave, 1 is a sine minus a pulse
//------------------------------------------------------------
// Author: Bart Brouns
// License: GPLv3
// CZ oscillators by Mike Moser-Booth:
// <https://forum.pdpatchrepo.info/topic/5992/casio-cz-oscillators>
// Ported from pd to Faust by Bart Brouns

CZsinePulsePAA(fund, index) = CZ.sinePulsePAA(fund, index);

//----------`(os.)CZhalfSine`----------
// Oscillator that mimics the Casio CZ half sine oscillator
// `CZhalfSine` is a standard Faust function.
//
// #### Usage
//
// ```
// CZhalfSine(fund,index) : _
// ```
//
// Where:
//
// * `fund`: a saw-tooth waveform between 0 and 1 that the oscillator slaves to
// * `index`: the brightness of the oscillator, 0 gives a sine-wave, 1 is somewhere between a saw and a square
//------------------------------------------------------------
// Author: Bart Brouns
// License: GPLv3
// CZ oscillators by Mike Moser-Booth:
// <https://forum.pdpatchrepo.info/topic/5992/casio-cz-oscillators>
// Ported from pd to Faust by Bart Brouns

CZhalfSine(fund, index) = CZ.halfSine(fund, index);

//----------`(os.)CZhalfSineAA`----------
// Oscillator that mimics the Casio CZ half sine oscillator,
// anti-aliased by decreasing the index at higher frequencies.
// `CZhalfSineAA` is a standard Faust function.
//
// #### Usage
//
// ```
// CZhalfSineAA(fund,index) : _
// ```
//
// Where:
//
// * `fund`: a saw-tooth waveform between 0 and 1 that the oscillator slaves to
// * `index`: the brightness of the oscillator, 0 gives a sine-wave, 1 is somewhere between a saw and a square
//------------------------------------------------------------
// Author: Bart Brouns
// License: GPLv3
// CZ oscillators by Mike Moser-Booth:
// <https://forum.pdpatchrepo.info/topic/5992/casio-cz-oscillators>
// Ported from pd to Faust by Bart Brouns

CZhalfSineAA(fund, index) = CZ.halfSineAA(fund, index);

//----------`(os.)CZhalfSineP`----------
// Oscillator that mimics the Casio CZ half sine oscillator,
// with it's phase aligned to `fund:sin`.
// `CZhalfSineP` is a standard Faust function.
//
// #### Usage
//
// ```
// CZhalfSineP(fund,index) : _
// ```
//
// Where:
//
// * `fund`: a saw-tooth waveform between 0 and 1 that the oscillator slaves to
// * `index`: the brightness of the oscillator, 0 gives a sine-wave, 1 is somewhere between a saw and a square
//------------------------------------------------------------
// Author: Bart Brouns
// License: GPLv3
// CZ oscillators by Mike Moser-Booth:
// <https://forum.pdpatchrepo.info/topic/5992/casio-cz-oscillators>
// Ported from pd to Faust by Bart Brouns

CZhalfSineP(fund, index) = CZ.halfSineP(fund, index);

//----------`(os.)CZhalfSinePAA`----------
// Oscillator that mimics the Casio CZ half sine oscillator,
// with it's phase aligned to `fund:sin`,
// anti-aliased by decreasing the index at higher frequencies.
// `CZhalfSinePAA` is a standard Faust function.
//
// #### Usage
//
// ```
// CZhalfSinePAA(fund,index) : _
// ```
//
// Where:
//
// * `fund`: a saw-tooth waveform between 0 and 1 that the oscillator slaves to
// * `index`: the brightness of the oscillator, 0 gives a sine-wave, 1 is somewhere between a saw and a square
//------------------------------------------------------------
// Author: Bart Brouns
// License: GPLv3
// CZ oscillators by Mike Moser-Booth:
// <https://forum.pdpatchrepo.info/topic/5992/casio-cz-oscillators>
// Ported from pd to Faust by Bart Brouns

CZhalfSinePAA(fund, index) = CZ.halfSinePAA(fund, index);

//----------`(os.)CZresSaw`----------
// Oscillator that mimics the Casio CZ resonant saw-tooth oscillator
// `CZresSaw` is a standard Faust function.
//
// #### Usage
//
// ```
// CZresSaw(fund,res) : _
// ```
//
// Where:
//
// * `fund`: a saw-tooth waveform between 0 and 1 that the oscillator slaves to
// * `res`: the frequency of resonance as a factor of the fundamental pitch.
//------------------------------------------------------------
// Author: Bart Brouns
// License: GPLv3
// CZ oscillators by Mike Moser-Booth:
// <https://forum.pdpatchrepo.info/topic/5992/casio-cz-oscillators>
// Ported from pd to Faust by Bart Brouns

CZresSaw(fund,res) = CZ.resSaw(fund,res);

//----------`(os.)CZresSawAA`----------
// Oscillator that mimics the Casio CZ resonant saw-tooth oscillator,
// anti-aliased by decreasing the resonance at higher frequencies.
// `CZresSawAA` is a standard Faust function.
//
// #### Usage
//
// ```
// CZresSawAA(fund,res) : _
// ```
//
// Where:
//
// * `fund`: a saw-tooth waveform between 0 and 1 that the oscillator slaves to
// * `res`: the frequency of resonance as a factor of the fundamental pitch.
//------------------------------------------------------------
// Author: Bart Brouns
// License: GPLv3
// CZ oscillators by Mike Moser-Booth:
// <https://forum.pdpatchrepo.info/topic/5992/casio-cz-oscillators>
// Ported from pd to Faust by Bart Brouns

CZresSawAA(fund,res) = CZ.resSawAA(fund,res);

//----------`(os.)CZresTriangle`----------
// Oscillator that mimics the Casio CZ resonant triangle oscillator
// `CZresTriangle` is a standard Faust function.
//
// #### Usage
//
// ```
// CZresTriangle(fund,res) : _
// ```
//
// Where:
//
// * `fund`: a saw-tooth waveform between 0 and 1 that the oscillator slaves to
// * `res`: the frequency of resonance as a factor of the fundamental pitch.
//------------------------------------------------------------
// Author: Bart Brouns
// License: GPLv3
// CZ oscillators by Mike Moser-Booth:
// <https://forum.pdpatchrepo.info/topic/5992/casio-cz-oscillators>
// Ported from pd to Faust by Bart Brouns

CZresTriangle(fund,res) = CZ.resTriangle(fund,res);

//----------`(os.)CZresTriangleAA`----------
// Oscillator that mimics the Casio CZ resonant triangle oscillator,
// anti-aliased by decreasing the resonance at higher frequencies.
// `CZresTriangleAA` is a standard Faust function.
//
// #### Usage
//
// ```
// CZresTriangleAA(fund,res) : _
// ```
//
// Where:
//
// * `fund`: a saw-tooth waveform between 0 and 1 that the oscillator slaves to
// * `res`: the frequency of resonance as a factor of the fundamental pitch.
//------------------------------------------------------------
// Author: Bart Brouns
// License: GPLv3
// CZ oscillators by Mike Moser-Booth:
// <https://forum.pdpatchrepo.info/topic/5992/casio-cz-oscillators>
// Ported from pd to Faust by Bart Brouns

CZresTriangleAA(fund,res) = CZ.resTriangleAA(fund,res);

//----------`(os.)CZresTrap`----------
// Oscillator that mimics the Casio CZ resonant trapeze oscillator
// `CZresTrap` is a standard Faust function.
//
// #### Usage
//
// ```
// CZresTrap(fund,res) : _
// ```
//
// Where:
//
// * `fund`: a saw-tooth waveform between 0 and 1 that the oscillator slaves to
// * `res`: the frequency of resonance as a factor of the fundamental pitch.
//------------------------------------------------------------
// Author: Bart Brouns
// License: GPLv3
// CZ oscillators by Mike Moser-Booth:
// <https://forum.pdpatchrepo.info/topic/5992/casio-cz-oscillators>
// Ported from pd to Faust by Bart Brouns

CZresTrap(fund, res) = CZ.resTrap(fund, res);

//----------`(os.)CZresTrapAA`----------
// Oscillator that mimics the Casio CZ resonant trapeze oscillator,
// anti-aliased by decreasing the resonance at higher frequencies.
// `CZresTrapAA` is a standard Faust function.
//
// #### Usage
//
// ```
// CZresTrapAA(fund,res) : _
// ```
//
// Where:
//
// * `fund`: a saw-tooth waveform between 0 and 1 that the oscillator slaves to
// * `res`: the frequency of resonance as a factor of the fundamental pitch.
//------------------------------------------------------------
// Author: Bart Brouns
// License: GPLv3
// CZ oscillators by Mike Moser-Booth:
// <https://forum.pdpatchrepo.info/topic/5992/casio-cz-oscillators>
// Ported from pd to Faust by Bart Brouns

CZresTrapAA(fund, res) = CZ.resTrapAA(fund, res);

CZ = environment {

       saw(fund, index) = sawChooseP(fund, index, 0);
       sawP(fund, index) = sawChooseP(fund, index, 1);
       sawAA(fund, index) = saw(fund, indexAA(fund, index));
       sawPAA(fund, index) = sawP(fund, indexAA(fund, index));
       sawChooseP(fund, index, p) =
         (((FUND(fund,align,p)*((.5-INDEX)/INDEX)),(-1*FUND(fund,align,p)+1)*((.5-INDEX)/(1-INDEX))):min+FUND(fund,align,p))*2*ma.PI:cos
       with {
         INDEX = (.5-(index*.5)):max(0.01):min(0.5);
         align = si.interpolate(index, 0.75, 0.5);
       };

       square(fund, index) = squareChooseP(fund, index, 0);
       squareP(fund, index) = squareChooseP(fund, index, 1);
       squareAA(fund, index) = square(fund, indexAA(fund, index));
       squarePAA(fund, index) = squareP(fund, indexAA(fund, index));
       squareChooseP(fund, index, p) = (FUND(fund,align,p)>=0.5), (ma.decimal((FUND(fund,align,p)*2)+1)<:_-min(_,(-1*_+1)*((INDEX)/(1-INDEX)))) :+ *ma.PI:cos
       with {
         INDEX = (index:pow(0.25)):max(0):min(1);
         align = si.interpolate(INDEX, -0.25, 0);
       };

       pulse(fund, index) = pulseChooseP(fund, index, 0);
       pulseP(fund, index) = pulseChooseP(fund, index, 1);
       pulseAA(fund, index) = pulse(fund, indexAA(fund, index));
       pulsePAA(fund, index) = pulseP(fund, indexAA(fund, index));
       pulseChooseP(fund, index, p) = ((FUND(fund,align,p)-min(FUND(fund,align,p),((-1*FUND(fund,align,p)+1)*(INDEX/(1-INDEX)))))*2*ma.PI):cos
       with {
         INDEX = index:min(0.99):max(0);
         align = si.interpolate(index, -0.25, 0.0);
       };

       sinePulse(fund, index) = sinePulseChooseP(fund, index, 0);
       sinePulseP(fund, index) = sinePulseChooseP(fund, index, 1);
       sinePulseAA(fund, index) = sinePulse(fund, indexAA(fund, index));
       sinePulsePAA(fund, index) = sinePulseP(fund, indexAA(fund, index));
       sinePulseChooseP(fund, index, p) = (min(FUND(fund,align,p)*((0.5-INDEX)/INDEX),(-1*FUND(fund,align,p)+1)*((.5-INDEX)/(1-INDEX)))+FUND(fund,align,p))*4*ma.PI:cos
       with {
         INDEX = ((index*-0.49)+0.5);
         align = si.interpolate(index, -0.125, -0.25);
       };

       halfSine(fund, index) = halfSineChooseP(fund, index, 0);
       halfSineP(fund, index) = halfSineChooseP(fund, index, 1);
       halfSineAA(fund, index) = halfSine(fund, indexAA(fund, index));
       halfSinePAA(fund, index) = halfSineP(fund, indexAA(fund, index));
       halfSineChooseP(fund, index, p) = (select2(FUND(fund,align,p)<.5, .5*(FUND(fund,align,p)-.5)/INDEX+.5, FUND(fund,align,p)):min(1))*2*ma.PI:cos
       with {
         INDEX = (.5-(index*0.5)):min(.5):max(.01);
         align = si.interpolate(index:min(0.975), -0.25, -0.5);
       };

       FUND =
         case {
           (fund,align,0) => fund;
           (fund,align,1) => (fund+align) : ma.frac; // align phase with fund
         };
       resSaw(fund,res) = (((-1*(1-fund))*((cos((ma.decimal((max(1,res)*fund)+1))*2*ma.PI)*-.5)+.5))*2)+1;
       resSawAA(fund, res) = resSaw(fund, resAA(fund, res));
       resTriangle(fund, res) = select2(fund<.5, 2-(fund*2), fund*2)*RES*2-1
       with {
         RES = ((fund*(res:max(1)))+1:ma.decimal)*2*ma.PI:cos*.5+.5;
       };
       resTriangleAA(fund, res) = resTriangle(fund, resAA(fund, res));
       resTrap(fund, res) = (((1-fund)*2):min(1)*sin(ma.decimal(fund*(res:max(1)))*2*ma.PI));
       resTrapAA(fund,  res) = resTrap(fund,resAA(fund,  res));

       fund2freq(fund)        = ((fund-fund')*ma.SR) : ba.sAndH(abs(fund-fund')<0.5);
       indexAA(fund,index) =  // Anti Alias => lower the index for higher freqs
         index*(1-
                (( (fund2freq(fund)-(ma.SR/256))
                   / (ma.SR/8))
                 :max(0):min(1)
                ));
       resAA(fund,res) = res*fund2freq(fund):max(0):min(ma.SR/4)/fund2freq(fund);
     };
