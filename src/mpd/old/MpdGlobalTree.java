package mpd.old;

import java.io.FileWriter;
import java.io.IOException;
import java.text.NumberFormat;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Comparator;
import java.util.HashMap;
import java.util.List;
import java.util.Locale;

import mpd.utils.Utils;


public class MpdGlobalTree {
	
	public String title;
	
	boolean sampling;
	FileWriter file;
	FileWriter outputFile;
	
	public int frequency = 100;
	double gValue;
	
	CurrentAlignment curAlig;
	
	ColumnNetwork3 network;
	Column3 firstVector, lastVector;
	int sizeOfAlignments;
	
	int[] firstDescriptor; 
	String t[][];
	String[] sequences;
	StringBuilder[] alignBuilder;
	
	int totalSamples;

	List<Double> decoding;
	String[] alignment;
	
	List<Double> fwdCondMarg;		// forward conditional marginals of the MPD columns
	List<Double> bwdCondMarg;		// backward ~
	public FileWriter fwdCondMargFile;
	public FileWriter bwdCondMargFile;
	
	public boolean writeScores = true;		// writes scores into MPD file after alignment
	
	static Comparator<String[]> compStringArr = new Comparator<String[]>() {
		@Override
		public int compare(String[] a1, String[] a2) {
			return a1[0].compareTo(a2[0]);
		}};
		
		
	public void setGValue(double gValue) {
		this.gValue = gValue;
	}
	
	public double getGValue() {
		return gValue;
	}
	
	public void setFrequency(int frequency) {
		this.frequency = frequency;
	}
	
	public int getFrequency() {
		return frequency;
	}
		
	public void beforeFirstSample(int sizeOfAlignments, FileWriter outputFile) {
		this.outputFile = outputFile;
		this.sizeOfAlignments = sizeOfAlignments;
		alignment = new String[sizeOfAlignments];
		t = new String[sizeOfAlignments][];
		alignBuilder = new StringBuilder[sizeOfAlignments];
		for(int i = 0; i < sizeOfAlignments; i++)
			alignBuilder[i] = new StringBuilder();
		sequences = null;

		network = new ColumnNetwork3(gValue);
		
		firstDescriptor = new int[sizeOfAlignments];
		Arrays.fill(firstDescriptor, -1);
		firstVector = network.add(firstDescriptor);
		lastVector = null;
	}

	public void newSample(int no, int total) {
		//System.out.println(curAlig);
		//System.out.println(curAlig.leafAlignment);
		totalSamples++;

		for(int i = 0; i < t.length; i++){
			t[i] = curAlig.leafAlignment[i].split("\t");
		}
		Arrays.sort(t, compStringArr);

		int[] previousDescriptor = firstDescriptor;
		
		int len = t[0][1].length();
		for(int j = 0; j < len; j++){
			int[] nextDescriptor =  new int[sizeOfAlignments];
			boolean allGap = true;
			for(int k = 0; k < sizeOfAlignments; k++){
				if(t[k][1].charAt(j) == '-')
					nextDescriptor[k] = ColumnKey3.colNext(previousDescriptor[k]);
				else {
					nextDescriptor[k] = ColumnKey3.colNext(previousDescriptor[k])+1;
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
				lastDescriptor[j] = ColumnKey3.colNext(previousDescriptor[j])+1;
			}
			lastVector = network.add(lastDescriptor);
		}
		
//		if(no == 0 || (total-1-no) % frequency == 0) {
		if(frequency != 0 && no % frequency == 0) {
			updateAll();
			
		}
		if(sampling) {
			try {
				String[] aln = alignment;
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
		network.updateViterbi(totalSamples);
		
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
		
		Column3 actualVector = firstVector.viterbi;
		Column3 prevVector = firstVector;
		fwdCondMarg = new ArrayList<Double>();
		bwdCondMarg = new ArrayList<Double>();
		decoding = new ArrayList<Double>();
		while(!actualVector.equals(lastVector)) {
			int[] desc = actualVector.key.desc;
			
			decoding.add((double)actualVector.count/totalSamples);
			for(int i = 0; i < desc.length; i++) {
				alignBuilder[i].append((desc[i] & 1) == 0 ? '-' : sequences[i].charAt(desc[i] >> 1));
			}
			fwdCondMarg.add(network.getFwdCondMarg(prevVector, actualVector));
			bwdCondMarg.add(network.getBwdCondMarg(actualVector, actualVector.viterbi));
			prevVector = actualVector;
			actualVector = actualVector.viterbi;
		}

		for(int i = 0; i < sizeOfAlignments; i++){
			alignment[i] = t[i][0]+"\t"+alignBuilder[i].toString();
		}
	}
	
	public void afterLastSample() {
		updateAll();
		try{
			String[] aln = alignment;
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
			
			if(fwdCondMargFile != null && fwdCondMarg != null) {
				for(double value : fwdCondMarg)
					fwdCondMargFile.write(value+"\n");
				fwdCondMargFile.close();
			}
			
			if(bwdCondMargFile != null && bwdCondMarg != null) {
				for(double value : bwdCondMarg)
					bwdCondMargFile.write(value+"\n");
				bwdCondMargFile.close();
			}
			
		}
		catch(IOException e){
		}
		
	}
	
	public void findIngroups(List<List<Integer>> splits, List<MpdGlobalTree> splitMpds, FileWriter fwdIngroupFile, FileWriter bwdIngroupFile) {
		List<Integer> fullSplit = new ArrayList<Integer>();
		for(int i = 0; i < firstVector.key.desc.length; i++)
			fullSplit.add(i);
		
		Column3 prevVector = firstVector;
		Column3 actualVector = firstVector.viterbi;
		while(!actualVector.equals(lastVector)) {
			List<Integer> fwdBestSplit = fullSplit, bwdBestSplit = fullSplit;
			double fwdBestCM = network.getFwdCondMarg(prevVector, actualVector);
			double bwdBestCM = network.getBwdCondMarg(actualVector, actualVector.viterbi);
			for(int i = 0; i < splits.size(); i++) {
				List<Integer> split = splits.get(i);
				ColumnNetwork3 splitNet = splitMpds.get(i).network;
				
				ColumnKey3 key = prevVector.key.split(split);
				Column3 prev = splitNet.contMap.get(key);
				if(prev == null && key != null) throw new Error("Key error!");
				
				key = actualVector.key.split(split);
				Column3 col = splitNet.contMap.get(key);
				if(col == null && key != null) throw new Error("Key error!");
				
//				System.out.println(""+actualVector.viterbi+Arrays.toString(actualVector.viterbi.key.split(split).desc));
				key = actualVector.viterbi.key.split(split);
				Column3 next = splitNet.contMap.get(key);
				if(next == null && key != null) throw new Error("Key error!");
//				{
//					int[] desc = actualVector.viterbi.key.desc;
//					StringBuilder s = new StringBuilder();
//
//					for(int j = 0; j < desc.length; j++) {
//						s.append((desc[j] & 1) == 0 ? '-' : sequences[j].charAt(desc[j] >> 1));
//					}
//					String x = s.toString();
//					desc = actualVector.viterbi.key.split(split).desc;
//					s = new StringBuilder();
//
//					for(int j = 0; j < desc.length; j++) {
//						s.append((desc[j] & 1) == 0 ? '-' : sequences[split.get(j)].charAt(desc[j] >> 1));
//					}
//					System.out.println(x+split+s);
//				}
//				System.out.println(""+prev+col+next);
				if(prev != null && col != null) {
					double cm = splitNet.getFwdCondMarg(prev, col);
					if(cm > fwdBestCM) {
						fwdBestCM = cm;
						fwdBestSplit = split;
					}
				}
				if(col != null && next != null) {
					double cm = splitNet.getBwdCondMarg(col, next);
					if(cm > bwdBestCM) {
						bwdBestCM = cm;
						bwdBestSplit = split;
					}
				}
			}
			try {
				fwdIngroupFile.write(Utils.joinStrings(fwdBestSplit.toArray()," ")+"\n");
				bwdIngroupFile.write(Utils.joinStrings(bwdBestSplit.toArray()," ")+"\n");
			} catch (IOException e) {
			}
			prevVector = actualVector;
			actualVector = actualVector.viterbi;
		}
		try {
			fwdIngroupFile.close();
			bwdIngroupFile.close();
		} catch (Exception e) {
		}
	}

	public void setSampling(boolean enabled, FileWriter file) {
		sampling = enabled;
		this.file = file;
	}
	
	private static class ColumnNetwork3 {
		HashMap<ColumnKey3, Column3> contMap = new HashMap<ColumnKey3, Column3>();
		HashMap<ColumnKey3, ArrayList<Column3>> succList = new HashMap<ColumnKey3, ArrayList<Column3>>();
		HashMap<ColumnKey3, ArrayList<Column3>> predList = new HashMap<ColumnKey3, ArrayList<Column3>>();

		int numberOfNodes = 0;
		Column3 first;

		double gValue;

		static NumberFormat usFormat = NumberFormat.getInstance(Locale.US);

		public ColumnNetwork3(double gValue) {
			this.gValue = gValue;
		}

		public double getFwdCondMarg(Column3 col1, Column3 col2) {
			int sum = 0;
			for(Column3 col : succList.get(col1.key.succ())) {
				sum += col.count;
			}
			return (double)col2.count/sum;
		}

		public double getBwdCondMarg(Column3 col1, Column3 col2) {
			int sum = 0;
			for(Column3 col : predList.get(col2.key.pred())) {
				sum += col.count;
			}
			return (double)col1.count/sum;
		}

		/**
		 * Adds a new alignment column into the network. If already in the network, column count is incremented.
		 * @param descriptor Alignment column represented by an array of signed integers
		 * @return alignment column or null if it was in network before
		 */
		Column3 add(int[] descriptor) {
			Column3 column;

			ColumnKey3 key = new ColumnKey3(descriptor);
			if((column=contMap.get(key)) != null) {
				++column.count;
				return null;
			}
			column = new Column3(key);
			contMap.put(key, column);
			if(numberOfNodes == 0)
				first = column;
			numberOfNodes++;

			ArrayList<Column3> arr;

			ColumnKey3 predKey = key.pred();
			if((arr = succList.get(predKey)) == null)
				succList.put(predKey, arr = new ArrayList<Column3>());
			arr.add(column);

			ColumnKey3 succKey = key.succ();
			if((arr = predList.get(succKey)) == null)
				predList.put(succKey, arr = new ArrayList<Column3>());
			arr.add(column);

			return column;
		}


		void updateViterbi(int n) {
			for(Column3 cols : contMap.values()) {
				cols.score = Double.NEGATIVE_INFINITY;
				cols.viterbi = null;
			}
			usFormat.setMaximumFractionDigits(3);
			System.out.println("Viterbi score: "+usFormat.format(updateViterbi(first, n)));
		}

		double updateViterbi(Column3 col, int n) {
			if(col.score != Double.NEGATIVE_INFINITY)
				return col.score;
			List<Column3> list = succList.get(col.key.succ());
			if(list == null) {
				return col.score = 0;
			} else {
				for(Column3 succ : list) {
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

	private static class Column3 {
		ColumnKey3 key;

		int count = 1;
		double score;
		Column3 viterbi;

		Column3(ColumnKey3 _key) {
			key = _key;
		}

		@Override
		public String toString() {
			return Arrays.toString(key.desc);
		}
	}

	private static class ColumnKey3 {
		public int[] desc;

		ColumnKey3(int[] arr) {
			desc = arr;
		}

		@Override
		public boolean equals(Object o)	{
			return (o instanceof ColumnKey3) && Arrays.equals(desc, ((ColumnKey3)o).desc);
		}

		@Override
		public int hashCode()	{
			return Arrays.hashCode(desc);
		}

		/**
		 * Calculates key of predecessor column class for this column.
		 * (Column x can directly precede column y iff succ(x)=pred(y) ) 
		 */
		ColumnKey3 pred() {
			int[] desc = this.desc, ret = new int[desc.length];
			for(int i = 0; i < desc.length; i++)
				ret[i] = desc[i] >> 1;
			return new ColumnKey3(ret);
		}

		/**
		 * Calculates key of successor column class for this column.
		 */
		ColumnKey3 succ() {
			int[] desc = this.desc, ret = new int[desc.length];
			for(int i = 0; i < desc.length; i++)
				ret[i] = (desc[i]+1) >> 1;
			return new ColumnKey3(ret);
		}

		ColumnKey3 split(List<Integer> split) {
			boolean allGap = true;
			int[] desc = this.desc, sdesc = new int[split.size()];
			for(int i = 0; i < split.size(); i++) {
				sdesc[i] = desc[split.get(i)];
				if((sdesc[i] & 1) == 1)
					allGap = false;
			}
			if(allGap)
				return null;
			return new ColumnKey3(sdesc);
		}

		static int colNext(int n) {
			return n + (n & 1);
		}

	}

}
