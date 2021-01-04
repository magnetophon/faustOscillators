declare author "Bart Brouns";
declare license "GPLv3";
declare name "formantOscs";
import("stdfaust.lib");

process =
  // PAF(center,fund,pafIndex,1);
  fof(phase,fofCenter,multiK,freq,fofSkirt,fofDecay,octaviation);

fund = os.lf_sawpos(freq);

fofCenter = hslider("center", 440, 27.5, 3520, 0.1):si.smoo;
phase = hslider("phase", 0, 0, 1, 0.001):si.smoo;
fofSkirt= hslider("[2]skirt[style:knob]", 2,0.01,10,0.001)*t(4)*(ma.SR/freq*ma.SR):si.smoo;
fofDecay = vslider("[3]decay[style:knob]", 400,2,2000,1)*multi:si.smoo;
octaviation = hslider("octaviation",0,0,maxOctavation,0.001):si.smoo;

multiK = os.lf_rawsaw(ma.SR/freq*multi);

// FOF
// the maximum number of octaves of the fof oscilators
// carefull: the CPU usage will go up with pow(2,maxOctavation) !
// and be prepared for extremely long compile times.
maxOctavation             = 4;
//-----------------------------------------------
// PAF oscilator
//-----------------------------------------------
// http://msp.ucsd.edu/techniques/v0.11/book-html/node96.html
//
// center = vslider("[01]bottom[style:knob][tooltip: the lowest formant frequency of the oscillators]",	1, 0.5, 7, 0.001):pow(2):si.smoo*freq;			//0.25 to 49 logarithmicly
center = hslider("center", 440, 27.5, 3520, 0.1)/freq :si.smoo;
freq = hslider("freq", 440, 27.5, 3520, 0.1):si.smoo;
pafIndex= hslider("[0]index[style:knob][tooltip: PAF index]",	100, 0.001, 100, 0.001);
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

FMfund(freq,index,phase) = fundPhase(freq,index,FMoctave,phase);

//-------------------------------------
//   A formant wave function generator
//-------------------------------------
import("stdfaust.lib");

// Based on
// https://ccrma.stanford.edu/~mjolsen/220a/fp/Foflet.dsp
// by Michael Jørgen Olsen

// modified for dynamic parameters and octavation by Bart Brouns

// Sample period
T = 1.0/ma.SR;
t(4) = 0.000001;
// functions used in foflet calculation
multi = 2:pow(maxOctavation);
expy(fund,BW) = exp(-BW * ma.PI * T)^(fund/multi); // exponential env (BW = BW)
// bug in original: we don't want * ma.SR.
envAttack(fund,beta) = 0.5 * (1.0 - cos(beta * (fund) )); // attack discontinuity smoother (beta=beta)
// envAttack(beta) = 0.5 * (1.0 - cos(beta * k * ma.SR)); // attack discontinuity smoother (y=beta)
// most interesting sounds with sine at phase 0 and envelope phased
sinus(fund,fc,phase) = sin((2.0 * ma.PI * fc * (fund) * T)+(phase*0*ma.PI)); // sinusoid (z=fc)

// functions to calculate fof attack and decay sections
// fofStop(BW) = k < sigLen(BW); // gate
fofAttack(fund,phase,BW,beta,fc) = expy(fund,BW) * envAttack(fund,beta) * sinus(fund,fc,phase); // first part of fof calculation
fofRemainder(fund,phase,BW,fc) =   expy(fund,BW) * sinus(fund,fc,phase); // 2nd part of fof calculation
// function to generate single fof
// most interesting sounds with sine at phase 0 and envelope phased:
// fofPart(fund,k1,phase,BW,fc) = (((fund) < int(k1)) * fofAttack(fund,phase,BW,beta,fc)) + (((fund) >= int(k1)) * fofRemainder(fund,phase,BW,fc)) with {
fofPart(fund,freq,k1,phase,BW,fc) = (((phasedFund) < int(k1)) * fofAttack(phasedFund,phase,BW,beta,fc)) + (((phasedFund) >= int(k1)) * fofRemainder(phasedFund,phase,BW,fc)) with {
  beta = ma.PI / (float(k1));
  phasedFund = (((fund/(multi*f0Period))+(phase)):ma.decimal)*f0Period*multi;
  f0Period = ma.SR/freq;
}; // v = k1
// combines multiple fofParts to create octavation:
// octaviation index, normally zero. If greater than fi.zero, lowers the effective xfund frequency by attenuating odd-numbered sinebursts. Whole numbers are full octaves, fractions transitional.
// adapted from: https://csound.github.io/docs/manual/fof2.html
fof(phase,fc,fund,freq,k1,BW,octaviation) =
  // fofPart(fundI(0),k1,phase,BW,fc)
  par(i,multi,fofPart(fundI(i),freq,k1,phase,BW,fc)*(OctMuliply(i):hbargraph("mult%i", 0, 1))):>_
with {
  fundI(i) = (((fund/(multi*f0Period))+(i/multi)):ma.decimal)*f0Period*multi;
  OctMuliply(i) =
    (OctMuliplyPart(i,floor(octaviation+1))*    ma.decimal(octaviation)) +
    (OctMuliplyPart(i,floor(octaviation  ))* (1-ma.decimal(octaviation)));
  OctMuliplyPart(i,oct) = select2(i>0,1,
                                  (((i % int(2:pow(oct)))):min(1)*-1)+1
                                 );

  f0Period = ma.SR/freq;
};
