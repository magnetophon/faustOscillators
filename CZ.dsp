declare author "Bart Brouns";
declare license "GPLv3";
declare name "CZoscs";
import("stdfaust.lib");

process(fund,index,res) =
 CZsawPAA(fund, index)
, CZsquarePAA(fund, index)
, CZpulsePAA(fund, index)
, CZsinePulsePAA(fund, index)
, CZhalfSinePAA(fund, index)
, CZresSawAA(fund,res)
, CZresTriangleAA(fund,res)
, CZresTrapAA(fund, res);

// params:
// osc type
// filter type
// filter freq
// filter q
// osc index
// osc octave
// osc phase


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

CZsawAA(fund, index) = CZ.sawPAA(fund, index);

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

CZsquareAA(fund, index) = CZ.square(fund, index);

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

CZpulseAA(fund, index) = CZ.pulse(fund, index);

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

CZsinePulseAA(fund, index) = CZ.sinePulse(fund, index);

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
    indexAA(index,fund) =  // Anti Alias => lower the index for higher freqs
      index*(1-
             (( (fund2freq(fund)-(ma.SR/256))
                / (ma.SR/8))
              :max(0):min(1)
             ));
    resAA(fund,res) = res*fund2freq(fund):max(0):min(ma.SR/4)/fund2freq(fund);
};
