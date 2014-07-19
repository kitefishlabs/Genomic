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
	chain = PlayBuf.ar(1, srcbufNum, BufRateScale.kr(srcbufNum) * transp, startPos: (start * BufSampleRate.kr(srcbufNum))) * env;
	chain = CombC.ar(chain, 1.0, c_delay, c_decay) * beta;
	chain = (chain * d_mult).clip2(d_amp) * gamma;
	chain = FFT(LocalBuf(8192), chain);
	chain = PV_MagSmear(chain, ms_bins);
	Out.ar(outbus, Pan2.ar(IFFT(chain), 0.0));
	// Out.ar(outbus, Pan2.ar(chain));

	// Out.ar(0, Pan2.ar(SinOsc.ar(0.1)*env));
}).load(s);


SynthDef(\sigmaSynth2, { |outbus=0, srcbufNum, start=0.0, dur=1.0, transp=1.0, c_delay=0.0, c_decay=0.0, d_mult=1.0, d_amp=0.7, ms_bins=0, alpha=1, beta=1, gamma=1|

	var env, chain, chainA, chainB, chainC;
	env = EnvGen.kr(Env.linen(0.01, ((dur / transp) - 0.02), 0.01, alpha), gate: 1, doneAction: 2);
	chain = PlayBuf.ar(1, srcbufNum, BufRateScale.kr(srcbufNum) * transp, startPos: (start * BufSampleRate.kr(srcbufNum))) * env;

	chainA = CombC.ar(chain, 1.0, c_delay, c_decay) * beta;

	chainB = (chain * d_mult).clip2(d_amp) * gamma;

	chainC = FFT(LocalBuf(8192), chain);
	chainC = PV_MagSmear(chainC, ms_bins);
	chainC = IFFT(chainC);

	Out.ar(outbus, Pan2.ar(Mix([chain, chainA, chainB, chainC]), 0.0));
}).load(s);

//

b = Buffer.read(s, "/Users/kfl/dev/git/Genomic/snd/viola_07.wav");

Synth(\sigmaSynth2, ["srcbufNum", b.bufnum, "start", 0.0, "dur", b.duration, "transp", 0.7, "c_delay", 0.1, "c_decay", 1.0, "d_mult", 13, "d_amp", 0.01, "ms_bins", 6, "alpha", 1, "beta", 1, "gamma", 1]);

synth

SynthDef(\sigmaSynthAnalyzer, { |outbus=0, srcbufNum=0, savebufNum=1,
	hop=0.04, start=0.0, dur=1.0, transp=1.0, c_delay=0.0, c_decay=0.0,
	d_mult=1.0, d_amp=1.0, ms_bins=0, alpha=1, beta=1, gamma=1|

	var env, chain, power, mfccs, driver;
	env = EnvGen.kr(
		Env.linen(
			0.01,
			((BufDur.kr(srcbufNum) / transp) - 0.02),
			0.01, alpha),
		gate: 1, doneAction: 2);
	chain = PlayBuf.ar(
		1,
		srcbufNum,
		BufRateScale.kr(srcbufNum) * transp,
		startPos: (start * BufSampleRate.kr(srcbufNum))) * env;

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

SynthDef(\sigmaSynthAnalyzer, { |outbus=0, srcbufNum=0, savebufNum=1,
	hop=0.04, start=0.0, dur=1.0, transp=1.0, c_delay=0.0, c_decay=0.0,
	d_mult=1.0, d_amp=1.0, ms_bins=0, alpha=1, beta=1, gamma=1|

	var env, chain, power, mfccs, driver;
	env = EnvGen.kr(
		Env.linen(
			0.01,
			((BufDur.kr(srcbufNum) / transp) - 0.02),
			0.01, alpha),
		gate: 1, doneAction: 2);
	chain = PlayBuf.ar(
		1,
		srcbufNum,
		BufRateScale.kr(srcbufNum) * transp,
		startPos: (start * BufSampleRate.kr(srcbufNum))) * env;

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


SynthDef(\sigmaSynthAnalyzer2, { |outbus=0, srcbufNum=0, savebufNum=1,
	hop=0.04, start=0.0, dur=1.0, transp=1.0, c_delay=0.0, c_decay=0.0,
	d_mult=1.0, d_amp=1.0, ms_bins=0, alpha=1, beta=1, gamma=1, delta=1|

	var env, chain, chainA, chainB, chainC, chainSigma, power, mfccs, driver;
	env = EnvGen.kr(
		Env.linen(
			0.01,
			((BufDur.kr(srcbufNum) / transp) - 0.02),
			0.01, alpha),
		gate: 1, doneAction: 2);
	chain = PlayBuf.ar(
		1,
		srcbufNum,
		BufRateScale.kr(srcbufNum) * transp,
		startPos: (start * BufSampleRate.kr(srcbufNum))) * env;

	chainA = CombC.ar(chain, 1.0, c_delay, c_decay) * beta;

	chainB = (chain * d_mult).clip2(d_amp) * gamma;

	chainC = FFT(LocalBuf(8192), chain);
	chainC = PV_MagSmear(chainC, ms_bins);
	chainC = IFFT(chainC) * delta;

	chainSigma = Mix([chainA, chainB, chainC]);

	power =			FFTPower.kr(chainSigma);
	mfccs =			MFCC.kr(chainSigma, numcoeff:24); // or 13|24|42...

	// log the metadata into a buffer and signal sclang to read from the buffer
	driver = Impulse.kr( 1 / hop );
	Logger.kr([power, mfccs].flatten, driver, savebufNum);

	Out.ar(outbus, IFFT(chain));
}).load(s);


////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////
SynthDef(\monoSamplerNRT, { |outbus=20, srcbufNum, start=0, dur=1, transp=1|
	var env, chain;
	env = EnvGen.ar(Env.linen(0.01, ((BufDur.kr(srcbufNum) / transp) - 0.02), 0.01, 1), gate: 1, doneAction: 2);
	Out.ar(outbus, PlayBuf.ar(1, srcbufNum, BufRateScale.kr(srcbufNum) * transp, startPos: (start * BufSampleRate.kr(srcbufNum))) * env);
}).load(s);


SynthDef(\sigmaMod, { |outbus=21, inbus=20, dur=1, transp=1, c_delay=0.0, c_decay=0.0, d_mult=1.0, d_amp=1.0, ms_bins=0, alpha=1, beta=1, gamma=1|
	var chain, env = EnvGen.ar(Env.linen(0.01, ((dur / transp) - 0.02), 0.01, alpha), gate: 1, doneAction: 2);
	chain = In.ar(inbus, 1);
	chain = CombC.ar(chain, 2.0, c_delay, c_decay) * beta;
	chain = (chain * d_mult).clip2(d_amp) * gamma;
	chain = FFT(LocalBuf(8192), chain);
	chain = PV_MagSmear(chain, ms_bins);
	chain = IFFT(chain);
	Out.ar(outbus, chain * env);
}).load(s);

SynthDef(\power_mfcc24BusAnalyzerNRT, { |inbus=21, savebufNum=0, hop=0.04|
	var in, chain, power, mfccs, driver, array;
	in = In.ar(inbus, 1);
	chain = FFT(LocalBuf(8192,1), in);

	power =			FFTPower.kr(chain);          // empirical multiplier
	mfccs =			MFCC.kr(chain, numcoeff:24); // or 13|24|42...

	// log the metadata into a buffer and signal sclang to read from the buffer
	driver = Impulse.kr( 1 / hop );
	Logger.kr([power, mfccs].flatten, driver, savebufNum);
	Out.ar(0, in);
}).load(s);

)

b = Buffer.readChannel(s, "/Users/kfl/dev/git/Genomic/snd/test0.wav", channels:[0]);
a = Buffer.alloc(s, b.numFrames, 1);

~sigma = Synth(\sigmaSynthAnalyzer, [\srcbufNum, b, \savebufNum, a, \dur, b.duration, \c_delay, 0.0, \c_decay, 0.0, \d_mult, 1.0, \d_amp, 1.0, \ms_bins, 0]);
