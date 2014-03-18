package wvalign.old;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

import wvalign.io.FastaReader;
import wvalign.io.RawSequences;
import wvalign.io.SampleReader;


public class MpdInterfaceTree {
	
	private static final String DEFAULT_ALIGN_TYPE = "StatAlign";
	
	private CurrentAlignment alig;
	private MpdGlobalTree mpd;
	
	private List<MpdGlobalTree> splitMpds = new ArrayList<MpdGlobalTree>();
	private List<List<Integer>> splits = new ArrayList<List<Integer>>();
	private String fwdIngroupFile;
	private String bwdIngroupFile;
	
	private int maxNoSamples;
	private int sampleRate = 1;
	
	private SampleReader sReader;
	private RawSequences lastSample;
	
	/**
	 * Logging (partial MPD results) is off by default.
	 * Call {@link #enableLogging(File)} to change this.
	 */
	public MpdInterfaceTree(double gValue) {
		alig = new CurrentAlignment();
		mpd = new MpdGlobalTree();
		mpd.setGValue(gValue);
		System.out.println("Using g value "+gValue);
	}
	
	public void addSplits(List<List<Integer>> splits, String fwdIngroupFile, String bwdIngroupFile) {
		this.splits = splits;
		this.fwdIngroupFile = fwdIngroupFile;
		this.bwdIngroupFile = bwdIngroupFile;
		for(int i = 0; i < splits.size(); i++) {
			MpdGlobalTree splitMpd = new MpdGlobalTree();
			splitMpd.setGValue(mpd.getGValue());
			splitMpds.add(splitMpd);
		}
	}
	
	public void setMaxNoSamples(int value) {
		maxNoSamples = value;
	}

	public void setSampleRate(int sampleRate) {
		this.sampleRate = sampleRate;
	}
	
	private void initMpd(String logFile, String outputFile) throws FileNotFoundException, IOException {
		sReader = new SampleReader(new FileReader(logFile));
		getNextSample();		// read first sample to get alignment size
		
		if(lastSample == null)
			throw new Error("No samples found.");
		
		System.out.println("Alignment size is "+lastSample.sequences.size());
		
		mpd.beforeFirstSample(lastSample.sequences.size(), new FileWriter(outputFile));
		
		for(int i = 0; i < splits.size(); i++) {
			MpdGlobalTree splitMpd = splitMpds.get(i);
			splitMpd.beforeFirstSample(splits.get(i).size(), new FileWriter("NUL"));
		}
	}
	
	private void getNextSample() {
		if(!sReader.isEof()) {
			try {
				FastaReader fReader = new FastaReader();
				lastSample = fReader.read(sReader);
				sReader.nextSample();
				return;
			} catch (IOException e) {
			}
		}
		lastSample = null;
	}
	
	public MpdGlobalTree getMpd() {
		return mpd;
	}
	
	/**
	 * Call this before {@link #doMpd()}.
	 * 
	 * @param logOutput file to write MPD log into
	 */
	public void enableLogging(String logFile) throws IOException {
		mpd.file = new FileWriter(logFile);
		mpd.sampling = true;
	}

	public void setCondMargFiles(String fwdFile, String bwdFile) throws IOException {
		mpd.fwdCondMargFile = new FileWriter(fwdFile);
		mpd.bwdCondMargFile = new FileWriter(bwdFile);
	}
	
	public void setWriteScores(boolean value) {
		mpd.writeScores = value;
	}
	
//	public void setAlignType(String type) {
//		mpd.alignmentType = type;
//	}
	
	public void doMpd(String logFile, String outputFile) throws FileNotFoundException, IOException {
		initMpd(logFile, outputFile);
		
		int no = 0;
		while(lastSample != null && (maxNoSamples == 0 || no < maxNoSamples)) {
//			System.out.println("Sample no. "+(no+1));
			injectSampleIntoAlig();
			mpd.newSample(no, no+1);
			String[] leafAligSorted = alig.leafAlignment;
			Arrays.sort(leafAligSorted);
			for(int i = 0; i < splits.size(); i++) {
//				System.out.println("dompd split "+i);
				injectSampleIntoAlig(leafAligSorted, splits.get(i));
				splitMpds.get(i).newSample(no, no+1);
			}
			for(int i = 0; i < sampleRate; i++)
				getNextSample();
			no++;
		}
		mpd.afterLastSample();
		if(fwdIngroupFile != null && bwdIngroupFile != null)
			try {
				mpd.findIngroups(splits, splitMpds, new FileWriter(fwdIngroupFile), new FileWriter(bwdIngroupFile));
			} catch (IOException e) {
			}
	}

	private void injectSampleIntoAlig() {
		int size = lastSample.sequences.size();
		
		alig.leafAlignment = new String[size];
		for(int i = 0; i < size; i++) {
			alig.leafAlignment[i] = lastSample.seqNames.get(i)+"\t"+lastSample.sequences.get(i);
		}
	}
	
	private void injectSampleIntoAlig(String[] leafAligSorted, List<Integer> split) {
		int size = split.size();

		alig.leafAlignment = new String[size];
		for(int i = 0; i < size; i++) {
			int ind = split == null ? i : split.get(i);
			alig.leafAlignment[i] = leafAligSorted[ind];
		}
	}
}
