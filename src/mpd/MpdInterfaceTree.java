package mpd;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

import mpd.io.SampleReader;

import statalign.base.Mcmc;
import statalign.base.Tree;
import statalign.base.Vertex;
import statalign.io.RawSequences;
import statalign.io.input.plugins.FastaReader;
import statalign.postprocess.Postprocess;

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
		
		Tree tree = new Tree();
		tree.vertex = new Vertex[2*lastSample.sequences.size()-1];
		mpd.mcmc = new Mcmc(tree, null, null);
		
		mpd.refToDependences(new Postprocess[] { alig });
		mpd.show = false;
		mpd.sampling = false;
		mpd.outputFile = new FileWriter(outputFile);
		mpd.frequency = 0;
		mpd.beforeFirstSample();
		
		for(int i = 0; i < splits.size(); i++) {
			tree = new Tree();
			tree.vertex = new Vertex[2*splits.get(i).size()-1];
			MpdGlobalTree splitMpd = splitMpds.get(i);
			splitMpd.mcmc = new Mcmc(tree, null, null);
			
			splitMpd.refToDependences(new Postprocess[] { alig });
			splitMpd.show = false;
			splitMpd.sampling = false;
			splitMpd.frequency = 0;
			splitMpd.alignmentType = DEFAULT_ALIGN_TYPE;
			splitMpd.outputFile = new FileWriter("NUL");	// null writer
//			System.out.println("beforefirst mpd "+i);
			splitMpd.beforeFirstSample();
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
	
	public Postprocess getMpd() {
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
	
	public void setAlignType(String type) {
		mpd.alignmentType = type;
	}
	
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
