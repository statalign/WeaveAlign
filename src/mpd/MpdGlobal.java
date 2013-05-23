package mpd;

import java.awt.BorderLayout;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Comparator;
import java.util.HashMap;
import java.util.List;

import javax.swing.Icon;
import javax.swing.ImageIcon;
import javax.swing.JPanel;
import javax.swing.JScrollPane;

import statalign.base.Utils;
import statalign.postprocess.Postprocess;
import statalign.postprocess.gui.AlignmentGUI;

public class MpdGlobal extends statalign.postprocess.Postprocess {

	public String title;
	AlignmentGUI gui;
	JPanel pan;
	
	public int frequency = 100;
	double gValue;
	
	CurrentAlignment curAlig;
	
	ColumnNetwork network;
	Column firstVector, lastVector;
	int sizeOfAlignments;
	
	int[] firstDescriptor; 
	String t[][];
	String[] sequences;
	StringBuilder[] alignBuilder;

	int totalSamples;
	
	List<Double> decoding;
	String[] alignment;
	
	boolean writeScores;
	
	long buildTime;		// time spent in building the network (creating columns, hashing etc.)
	long viterbiTime;	// time spent in Viterbi algorithm

	public MpdGlobal(){
		screenable = true;
		outputable = true;
		postprocessable = true;
		postprocessWrite = true;
		sampling = true;
	}
	
	@Override
	public JPanel getJPanel() {
		pan = new JPanel(new BorderLayout());
		return pan;
	}

	@Override
	public Icon getIcon() {
		return new ImageIcon(ClassLoader.getSystemResource("icons/MPD.gif"));
	}

	@Override
	public String getTabName() {
		return "MPD";
	}

	@Override
	public String getTip() {
		return "Maximum Posterior Decoding (consensus) alignment";
	}

	@Override
	public String getFileExtension() {
		return "mpd";
	}

	@Override
	public String[] getDependences() {
		return new String[] { "statalign.postprocess.plugins.CurrentAlignment" };
	}


	@Override
	public void refToDependences(Postprocess[] plugins) {
		curAlig = (CurrentAlignment) plugins[0];
	}
	
	static Comparator<String[]> compStringArr = new Comparator<String[]>() {
		public int compare(String[] a1, String[] a2) {
			return a1[0].compareTo(a2[0]);
		}};
		
		
	public void setGValue(double gValue) {
		this.gValue = gValue;
	}
	
	public double getGValue() {
		return gValue;
	}
		
	@Override
	public void setSampling(boolean enabled) {
		sampling = enabled;
	}
	
	public long getBuildTime() {
		return buildTime;
	}
	
	public long getViterbiTime() {
		return viterbiTime;
	}

	@Override
	public void beforeFirstSample() {
		if(show) {
			pan.removeAll();
			title = mcmc.tree.title;
			JScrollPane scroll = new JScrollPane();
			scroll.setViewportView(gui = new AlignmentGUI(title,mcmc.tree.substitutionModel));//, mcmc.tree.printedAlignment()));
			pan.add(scroll, BorderLayout.CENTER);
			pan.getParent().validate();
		}
		
		sizeOfAlignments = (mcmc.tree.vertex.length+1)/2;
		alignment = new String[sizeOfAlignments];
		if(show)
			gui.alignment = alignment;
		t = new String[sizeOfAlignments][];
		alignBuilder = new StringBuilder[sizeOfAlignments];
		for(int i = 0; i < sizeOfAlignments; i++)
			alignBuilder[i] = new StringBuilder();
		sequences = null;

		network = new ColumnNetwork(gValue);
		
		firstDescriptor = new int[sizeOfAlignments];
		Arrays.fill(firstDescriptor, -1);
		firstVector = network.add(firstDescriptor);
		lastVector = null;
	}

	@Override
	public void newSample(int no, int total) {
		//System.out.println(curAlig);
		//System.out.println(curAlig.leafAlignment);
		totalSamples++;
		
		for(int i = 0; i < t.length; i++){
			t[i] = curAlig.leafAlignment[i].split("\t");
		}
		Arrays.sort(t, compStringArr);

		int[] previousDescriptor = firstDescriptor;
		
		buildTime -= System.currentTimeMillis();
		
		int len = t[0][1].length();
		for(int j = 0; j < len; j++){
			int[] nextDescriptor =  new int[sizeOfAlignments];
			boolean allGap = true;
			for(int k = 0; k < sizeOfAlignments; k++){
				if(t[k][1].charAt(j) == '-')
					nextDescriptor[k] = ColumnKey.colNext(previousDescriptor[k]);
				else {
					nextDescriptor[k] = ColumnKey.colNext(previousDescriptor[k])+1;
					allGap = false;
				}
			}
			if(!allGap)
				network.add(nextDescriptor);//[j]);
			
			previousDescriptor = nextDescriptor;
		}//j (length of alignments)
		
		if(no == 0) {		// add last vector once only
			int[] lastDescriptor =  new int[sizeOfAlignments];
			for(int j = 0; j < sizeOfAlignments; j++){
				lastDescriptor[j] = ColumnKey.colNext(previousDescriptor[j])+1;
			}
			lastVector = network.add(lastDescriptor);
		}
		
		buildTime += System.currentTimeMillis();
		
//		if(no == 0 || (total-1-no) % frequency == 0) {
		if(frequency != 0 && no % frequency == 0) {
			updateAll();
			
			if(show) {
				gui.decoding = new double[decoding.size()];
				for(int i = 0; i < decoding.size(); i++)
					gui.decoding[i] = decoding.get(i);
				gui.alignment = alignment;
				gui.repaint();
			}
		}
		if(sampling) {
			try {
				String[] aln = Utils.alignmentTransformation(alignment, alignmentType, mcmc.tree);
				for(int i = 0; i < aln.length; i++){
					file.write("Sample "+no+"\tMPD alignment:\t"+aln[i]+"\n");
				}
				if(decoding != null){
					for(int i = 0; i < decoding.size(); i++){
						file.write("Sample "+no+"\tMPD alignment probabilities:\t"+decoding.get(i)+"\n");
					}
				}
				else{
					file.write("Sample "+no+"\tMPD alignment:\tNo posterior values so far\n");
				}
			} catch (IOException e) {
				// TODO Auto-generated catch block
				e.printStackTrace();
			}

		}
	}

	private void updateAll() {
		viterbiTime -= System.currentTimeMillis();
		network.updateViterbi(totalSamples);
		viterbiTime += System.currentTimeMillis();
//		System.out.println("viterbi (slow): "+time);
		
		if(sequences == null) {
			sequences = new String[sizeOfAlignments];
			StringBuilder b = new StringBuilder();
			int len = t[0][1].length();
			for(int i = 0; i < sizeOfAlignments; i++){
				b.setLength(0);
				for(int j = 0; j < len; j++){
					if(t[i][1].charAt(j) != '-'){
						b.append(t[i][1].charAt(j));
					}
				}
				sequences[i] = b.toString();
			}
		}
		for(int i = 0; i < sizeOfAlignments; i++)
			alignBuilder[i].setLength(0);
		
		Column actualVector = firstVector.viterbi;
		decoding = new ArrayList<Double>();
		while(!actualVector.equals(lastVector)) {
			int[] desc = actualVector.key.desc;
			
			decoding.add((double)actualVector.count/totalSamples);
			for(int i = 0; i < desc.length; i++) {
				alignBuilder[i].append((desc[i] & 1) == 0 ? '-' : sequences[i].charAt(desc[i] >> 1));
			}
			actualVector = actualVector.viterbi;
		}

		for(int i = 0; i < sizeOfAlignments; i++){
			alignment[i] = t[i][0]+"\t"+alignBuilder[i].toString();
		}
	}
	
	@Override
	public void afterLastSample() {
		updateAll();
		try{
			String[] aln = Utils.alignmentTransformation(alignment, alignmentType, mcmc.tree);
			for(int i = 0; i < aln.length; i++){
				outputFile.write(aln[i]+"\n");
			}
			if(writeScores) {
				outputFile.write("\n#scores\n\n");
				if(decoding != null){
					for(int i = 0; i < decoding.size(); i++){
						outputFile.write(decoding.get(i)+"\n");
					}
				} else {
					outputFile.write("No posterior values so far\n");
				}
			}
			outputFile.close();
		}
		catch(IOException e){
		}
		
	}

	private static class ColumnNetwork {
		HashMap<ColumnKey, Column> contMap = new HashMap<ColumnKey, Column>();
		HashMap<ColumnKey, ArrayList<Column>> succList = new HashMap<ColumnKey, ArrayList<Column>>();

		int numberOfNodes = 0;
		Column first;

		double gValue;

		public ColumnNetwork(double gValue) {
			this.gValue = gValue;
		}

		/**
		 * Adds a new alignment column into the network. If already in the network, column count is incremented.
		 * @param descriptor Alignment column represented by an array of signed integers
		 * @return alignment column or null if it was in network before
		 */
		Column add(int[] descriptor) {
			Column column;

			ColumnKey key = new ColumnKey(descriptor);
			if((column=contMap.get(key)) != null) {
				++column.count;
				return null;
			}

			column = new Column(key);
			contMap.put(key, column);
			if(numberOfNodes == 0)
				first = column;
			numberOfNodes++;

			ArrayList<Column> arr;

			ColumnKey predKey = key.pred();
			if((arr = succList.get(predKey)) == null)
				succList.put(predKey, arr = new ArrayList<Column>());
			arr.add(column);

			return column;
		}

		void updateViterbi(int n) {
			for(Column cols : contMap.values()) {
				cols.score = Double.NEGATIVE_INFINITY;
				cols.viterbi = null;
			}
//			double score = 
			updateViterbi(first, n);
//			System.out.println("viterbi score: "+score);
		}

		double updateViterbi(Column col, int n) {
			if(col.score != Double.NEGATIVE_INFINITY)
				return col.score;
			List<Column> list = succList.get(col.key.succ());
			if(list == null) {
				return col.score = 0;
			} else {
				for(Column succ : list) {
					double sc = updateViterbi(succ, n);
					if(sc > col.score) {
						col.score = sc;
						col.viterbi = succ;
					}
				}
				return col.score += (double)col.count/n - gValue;
				//			return col.score += Math.log((double)col.count/n);
			}
		}

	}

	private static class Column {
		ColumnKey key;

		int count = 1;
		double score;
		Column viterbi;

		Column(ColumnKey _key) {
			key = _key;
		}

		@Override
		public String toString() {
			return Arrays.toString(key.desc);
		}
	}

	private static class ColumnKey {
		public int[] desc;

		ColumnKey(int[] arr) {
			desc = arr;
		}

		@Override
		public boolean equals(Object o)	{
			return (o instanceof ColumnKey) && Arrays.equals(desc, ((ColumnKey)o).desc);
		}

		@Override
		public int hashCode()	{
			return Arrays.hashCode(desc);
		}

		/**
		 * Calculates key of predecessor column class for this column.
		 * (Column x can directly precede column y iff succ(x)=pred(y) ) 
		 */
		ColumnKey pred() {
			int[] desc = this.desc, ret = new int[desc.length];
			for(int i = 0; i < desc.length; i++)
				ret[i] = desc[i] >> 1;
			return new ColumnKey(ret);
		}

		/**
		 * Calculates key of successor column class for this column.
		 */
		ColumnKey succ() {
			int[] desc = this.desc, ret = new int[desc.length];
			for(int i = 0; i < desc.length; i++)
				ret[i] = (desc[i]+1) >> 1;
			return new ColumnKey(ret);
		}

		static int colNext(int n) {
			return n + (n & 1);
		}

	}

}