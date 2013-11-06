(
// hop = 40 msec == 0.04 seconds == 1 / 25 fps == 25 Hz analysis frequency
SynthDef(\power_mfcc24BusAnalyzerNRT, { |inbus=20, savebufNum=0, hop=0.04|
	var in, chain, power, mfccs, driver;
	in = In.ar(inbus, 1);
	chain = FFT(LocalBuf(8192,1), in);

	power =			FFTPower.kr(chain);          // empirical multiplier
	mfccs =			MFCC.kr(chain, numcoeff:24); // or 13|24|42...

	// log the metadata into a buffer and signal sclang to read from the buffer
	driver = Impulse.kr( 1 / hop );
	Logger.kr([power, mfccs].flatten, driver, savebufNum);
	Out.ar(0, in);
}).load(s);

SynthDef(\sigmaSynth, { |outbus=0, srcbufNum, start=0.0, dur=1.0, transp=1.0, c_delay=0.0, c_decay=0.0, d_mult=1.0, d_amp=0.7, ms_bins=0, alpha=1, beta=1, gamma=1|

	var env, chain;
	env = EnvGen.kr(Env.linen(0.01, ((dur / transp) - 0.02), 0.01, alpha), gate: 1, doneAction: 2);
	//chain = PlayBuf.ar(1, srcbufNum, BufRateScale.kr(srcbufNum) * transp, startPos: (start * BufSampleRate.kr(srcbufNum))) * env;
	/*chain = CombC.ar(chain, 1.0, c_delay, c_decay) * beta;
	chain = (chain * d_mult).clip2(d_amp) * gamma;
	chain = FFT(LocalBuf(8192), chain);
	chain = PV_MagSmear(chain, ms_bins);
	Out.ar(outbus, Pan2.ar(IFFT(chain), 0.0));
	Out.ar(outbus, Pan2.ar(chain));
	*/
	Out.ar(0, Pan2.ar(SinOsc.ar(0.1)*env));
}).load(s);

SynthDef(\sigmaSynthAnalyzer, { |outbus=0, srcbufNum=0, savebufNum=1, hop=0.04, start=0.0, dur=1.0, transp=1.0, c_delay=0.0, c_decay=0.0, d_mult=1.0, d_amp=0.7, ms_bins=0, alpha=1, beta=1, gamma=1|

	var env, chain, power, mfccs, driver;
	env = EnvGen.kr(Env.linen(0.01, ((BufDur.kr(srcbufNum) / transp) - 0.02), 0.01, alpha), gate: 1, doneAction: 2);
	chain = PlayBuf.ar(1, srcbufNum, BufRateScale.kr(srcbufNum) * transp, startPos: (start * BufSampleRate.kr(srcbufNum))) * env;
	chain = CombC.ar(chain, 1.0, c_delay, c_decay) * beta;
	chain = (chain * d_mult).clip2(d_amp) * gamma;
	chain = FFT(LocalBuf(8192), chain);
	chain = PV_MagSmear(chain, ms_bins);

	power =			FFTPower.kr(chain);          // empirical multiplier
	mfccs =			MFCC.kr(chain, numcoeff:24); // or 13|24|42...

	// log the metadata into a buffer and signal sclang to read from the buffer
	driver = Impulse.kr( 1 / hop );
	Logger.kr([power, mfccs].flatten, driver, savebufNum);

	Out.ar(outbus, IFFT(chain));
}).load(s);
)