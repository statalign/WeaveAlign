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


public class DagInterface {
	
	private AlignmentDAG dag;
	
	private int maxNoSamples;
	private int sampleRate = 1;
	private int firstSample;
	private int nSamples = 0;
	
	public int getNSamples() { return nSamples; }	

	private SampleReader sReader;
	private ArrayDeque<String> fastaList;
	private RawSequences lastSample;

	private boolean timeStats = true;
	
	/**
	 * Logging (partial MinRisk results) is off by default.
	 * Call {@link #enableLogging(File)} to change this.
	 */
	public DagInterface(double gValue, boolean optGi, boolean outGi) throws FileNotFoundException, IOException {
		dag = new AlignmentDAG(gValue, optGi, outGi);
	}

	public AlignmentDAG getDag() {
		return dag;
	}

	public void setMaxNoSamples(int value) {
		maxNoSamples = value;
	}

	public void setSampleRate(int sampleRate) {
		this.sampleRate = sampleRate;
	}
	
	public void setFirstSample(int firstSample) {
		this.firstSample = firstSample;
	}
	
	//>test
	//	public void setClassSizeFiles(String fwdFile, String bwdFile) throws IOException {
	//		dag.fwdClassFile = new FileWriter(fwdFile);
	//		dag.bwdClassFile = new FileWriter(bwdFile);
	//	}
	//<test
		
	public void setTimeStats(boolean timeStats) {
		this.timeStats = timeStats;
	}
	void activateTwoState(boolean t) {
		dag.columnNetwork.activateTwoState(t);
	}
	void computeEquivalenceClassFreqs() {
		dag.columnNetwork.computeEquivalenceClassFreqs();
	}
	double logNPaths() {
		return dag.columnNetwork.logNPaths();
	}
	private void initDag(List<String> inFiles, String outputFile, String scoreFile, 
			boolean scoreSamples) throws FileNotFoundException, IOException {
		if(inFiles.size() == 1)
			sReader = new SampleReader(new FileReader(inFiles.get(0)));
		else
			fastaList = new ArrayDeque<String>(inFiles);
		getNextSample();		// read first sample to get alignment size
		
		if(lastSample == null || lastSample.sequences.size() == 0)
			throw new Error("No samples found.");

		if(!scoreSamples) {
			System.out.println("Using g value "+dag.getGValue());
			System.out.println("Number of sequences = "+lastSample.sequences.size());
			//System.out.println("Number of samples = "+fastaList.size());
			
			dag.outputFile = outputFile;
			dag.scoreFile = scoreFile;
		}
	}
	
	//>test
	//	public void setClassSizeFiles(String fwdFile, String bwdFile) throws IOException {
	//		dag.fwdClassFile = new FileWriter(fwdFile);
	//		dag.bwdClassFile = new FileWriter(bwdFile);
	//	}
	//<test

	public void scoreSamples(String logFile, String outputFile, String scoreFile,
			boolean computePosterior) 
			throws FileNotFoundException, IOException {
		networkOperation(logFile,outputFile,scoreFile,true,false,computePosterior);
	}
	public void scoreSamples(List<String> logFile, String outputFile, String scoreFile,
			boolean computePosterior) 
			throws FileNotFoundException, IOException {
		networkOperation(logFile,outputFile,scoreFile,true,false,computePosterior);
	}
	public void setupNetwork(String logFile, String outputFile, String scoreFile) 
			throws FileNotFoundException, IOException {
		networkOperation(logFile,outputFile,scoreFile,false,false,false);
	}
	public void setupNetwork(List<String> logFile, String outputFile, String scoreFile) 
			throws FileNotFoundException, IOException {
		networkOperation(logFile,outputFile,scoreFile,false,false,false);
	}
	public void computeMinRisk(String logFile, String outputFile, String scoreFile) 
			throws FileNotFoundException, IOException {
		networkOperation(logFile,outputFile,scoreFile,false,true,false);
	}
	public void computeMinRisk(List<String> logFile, String outputFile, String scoreFile) 
			throws FileNotFoundException, IOException {
		networkOperation(logFile,outputFile,scoreFile,false,true,false);
	}
	public void networkOperation(String logFile, String outputFile, String scoreFile, 
			boolean scoreSamples, boolean computeMinRisk, boolean computePosterior) 
					throws FileNotFoundException, IOException {
		networkOperation(Arrays.asList(logFile), outputFile, scoreFile, scoreSamples,computeMinRisk,computePosterior);
	}
	
	public void networkOperation(List<String> inFiles, String outputFile, String scoreFile, 
			boolean scoreSamples, boolean computeMinRisk, boolean computePosterior) 
					throws FileNotFoundException, IOException {
		long totalTime = -System.currentTimeMillis();
		long ioTime = totalTime;
		initDag(inFiles, outputFile, scoreFile, scoreSamples);

		FileWriter writer = null;
		if(scoreSamples) {
			writer = new FileWriter(outputFile);
		}
		
		
		if (computePosterior) {
			dag.computeEquivalenceClassFreqs();
		}

		int no = 0;
		int sampleIndex = 0;
		
		while(lastSample != null && (maxNoSamples == 0 || no < maxNoSamples)) {
			//			System.out.println("Sample no. "+(no+1));
			if (sampleIndex >= firstSample) {
				String[] align = getAlign();
				ioTime += System.currentTimeMillis();									
				if (scoreSamples){
					dag.scoreSample(no, align, writer,computePosterior);
				}
				else dag.addAlignment(align);
				ioTime -= System.currentTimeMillis();
				no++;
			}
			for(int i = 0; i < sampleRate; i++) {
				getNextSample();
				nSamples++;
				sampleIndex++;
			}	
		}
		ioTime += System.currentTimeMillis();
		if(computeMinRisk) {
			dag.finalise();

			totalTime += System.currentTimeMillis();
			if(timeStats) {
				System.out.println("Time spent in:");
				System.out.println(" * Input and output  : "+ioTime+" ms");
				System.out.println(" * Building network  : "+dag.getBuildTime()+" ms");
				if(dag.annotator == null)
					System.out.println(" * Viterbi algorithm : "+dag.getViterbiTime()+" ms");
				else
					System.out.println(" * Annotation total  : "+dag.getAnnotTime()+" ms");
				System.out.println("Total time: "+totalTime+" ms");
			}
		} else if(scoreSamples) {
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
	//		dag.fwdClassFile = new FileWriter(fwdFile);
	//		dag.bwdClassFile = new FileWriter(bwdFile);
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
		dag.setAnnotator(annotator);
	}
	
	public MinRiskAnnotator getAnnotator() {
		return dag.annotator;
	}
	
}
