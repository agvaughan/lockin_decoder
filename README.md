# lockin_decoder

This is a software implementation of lockin amplifier, which estimates the power of a ~sinusoidally modulated signal in the presence of massive noise. [See eg. http://www.thinksrs.com/downloads/PDFs/ApplicationNotes/AboutLIAs.pdf]

One trick here is that we're trying to implement this without the reference signal for the modulation, and thus need to simultaneously estimate the demodulation frequency/phase as well as its amplitude modulation.

One annoyance is that, given our long recordings, 0.001% error in estimating the modulation frequency led to catastrophic phase precession. In the end we gave up trying to estimate the modulation frequency precisely, and simply used a phase demodulator to compensate for that phase precession.
 
This code is experimental, with a few unused code paths for deprecated approaches. It works well enough, though - probably better than the $4k instrument from SRS for our purposes.
