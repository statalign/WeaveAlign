package wvalign;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayDeque;
import java.util.Arrays;
import java.util.List;

import wvalign.io.FastaReader;
import wvalign.io.RawSequences;
import wvalign.io.SampleReader;


public class WeaveInterface {
	
	private MinRiskGlobalFast mpd;
	
	private int maxNoSamples;
	private int sampleRate = 1;

	private SampleReader sReader;
	private ArrayDeque<String> fastaList;
	private RawSequences lastSample;

	private boolean timeStats = true;
	
	/**
	 * Logging (partial MPD results) is off by default.
	 * Call {@link #enableLogging(File)} to change this.
	 */
	public WeaveInterface(double gValue, boolean optGi, boolean outGi) throws FileNotFoundException, IOException {
		mpd = new MinRiskGlobalFast(gValue, optGi, outGi);
	}

	public MinRiskGlobalFast getMpd() {
		return mpd;
	}

	public void setMaxNoSamples(int value) {
		maxNoSamples = value;
	}

	public void setSampleRate(int sampleRate) {
		this.sampleRate = sampleRate;
	}
	
	//>test
	//	public void setClassSizeFiles(String fwdFile, String bwdFile) throws IOException {
	//		mpd.fwdClassFile = new FileWriter(fwdFile);
	//		mpd.bwdClassFile = new FileWriter(bwdFile);
	//	}
	//<test
		
	public void setTimeStats(boolean timeStats) {
		this.timeStats = timeStats;
	}
	void activateTwoState(boolean t) {
		mpd.network.activateTwoState(t);
	}
	void computeEquivalenceClassFreqs() {
		mpd.network.computeEquivalenceClassFreqs();
	}
	double logNPaths() {
		return mpd.network.logNPaths();
	}
	private void initMpd(List<String> inFiles, String outputFile, String scoreFile, int scoreSamples) throws FileNotFoundException, IOException {
		if(inFiles.size() == 1)
			sReader = new SampleReader(new FileReader(inFiles.get(0)));
		else
			fastaList = new ArrayDeque<String>(inFiles);
		getNextSample();		// read first sample to get alignment size
		
		if(lastSample == null || lastSample.sequences.size() == 0)
			throw new Error("No samples found.");

		if(scoreSamples < 2) {
			System.out.println("Using g value "+mpd.getGValue());
			System.out.println("Number of sequences = "+lastSample.sequences.size());

			mpd.outputFile = outputFile;
			mpd.scoreFile = scoreFile;
		}
	}
	
	//>test
	//	public void setClassSizeFiles(String fwdFile, String bwdFile) throws IOException {
	//		mpd.fwdClassFile = new FileWriter(fwdFile);
	//		mpd.bwdClassFile = new FileWriter(bwdFile);
	//	}
	//<test

	public void doMpd(String logFile, String outputFile, String scoreFile, int scoreSamples, boolean computePosterior) throws FileNotFoundException, IOException {
		doMpd(Arrays.asList(logFile), outputFile, scoreFile, scoreSamples,computePosterior);
	}
	
	public void doMpd(List<String> inFiles, String outputFile, String scoreFile, int scoreSamples, boolean computePosterior) throws FileNotFoundException, IOException {
		long totalTime = -System.currentTimeMillis();
		long ioTime = totalTime;
		initMpd(inFiles, outputFile, scoreFile, scoreSamples);

		FileWriter writer = null;
		if(scoreSamples == 2) {
			writer = new FileWriter(outputFile);
		}
		
		int no = 0;
		
		if (computePosterior) {
			mpd.computeEquivalenceClassFreqs();
		}

		while(lastSample != null && (maxNoSamples == 0 || no < maxNoSamples)) {
			//			System.out.println("Sample no. "+(no+1));
			String[] align = getAlign();
			ioTime += System.currentTimeMillis();
			if(scoreSamples < 2)
				mpd.addAlignment(align);
			else {
				mpd.scoreSample(no, align, writer,computePosterior);
			}
			ioTime -= System.currentTimeMillis();
			for(int i = 0; i < sampleRate; i++)
				getNextSample();
			no++;
		}
		ioTime += System.currentTimeMillis();
		if(scoreSamples == 0) {
			mpd.finalise();

			totalTime += System.currentTimeMillis();
			if(timeStats) {
				System.out.println("Time spent in:");
				System.out.println(" * Input and output  : "+ioTime+" ms");
				System.out.println(" * Building network  : "+mpd.getBuildTime()+" ms");
				if(mpd.annotator == null)
					System.out.println(" * Viterbi algorithm : "+mpd.getViterbiTime()+" ms");
				else
					System.out.println(" * Annotation total  : "+mpd.getAnnotTime()+" ms");
				System.out.println("Total time: "+totalTime+" ms");
			}
		} else if(scoreSamples == 2) {
			writer.close();
		}

	}

	private void getNextSample() {
		try {
			if(sReader != null) {
				if(!sReader.isEof()) {
					FastaReader fReader = new FastaReader();
					lastSample = fReader.read(sReader);
					if(fReader.getErrors() > 0)
						throw new Error("Error reading samples");
					sReader.nextSample();
					return;
				}
			} else if(fastaList.size() > 0) {
				FastaReader fReader = new FastaReader();
				String file = fastaList.pollFirst();
//				System.out.println("Reading alignment from "+file);
				lastSample = fReader.read(file);
				if(fReader.getErrors() > 0)
					throw new Error("Error reading alignment "+file);
				return;
			}
		} catch (IOException e) {
		}
		lastSample = null;
	}
	
	//>test
	//	public void setClassSizeFiles(String fwdFile, String bwdFile) throws IOException {
	//		mpd.fwdClassFile = new FileWriter(fwdFile);
	//		mpd.bwdClassFile = new FileWriter(bwdFile);
	//	}
	//<test
		
	private String[] getAlign() {
		int size = lastSample.sequences.size();

		String[] align = new String[size];
		for(int i = 0; i < size; i++) {
			align[i] = lastSample.seqNames.get(i)+"\t"+lastSample.sequences.get(i);
		}
		return align;
	}

	public void setAnnotator(MinRiskAnnotator annotator) {
		mpd.setAnnotator(annotator);
	}
	
	public MinRiskAnnotator getAnnotator() {
		return mpd.annotator;
	}
	
}