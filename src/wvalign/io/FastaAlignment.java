package wvalign.io;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;

import org.apache.commons.lang3.tuple.ImmutablePair;
import org.apache.commons.lang3.tuple.Pair;

import static wvalign.io.GappedSequence.GAP;

public class FastaAlignment  implements Comparable<FastaAlignment>{

	private List<GappedSequence> alignment = new ArrayList<GappedSequence>();
	private int width = -1;
	private int homologPairNum = -1;
	private int nonHomologPairNum = -1;
	private int id; //matrix id when using randomly perturbated matrices
	private Map<Pair<Integer, Integer>,
				Set<Pair<Integer, Integer> > > homologyPairs = null;
	
	public int getNumOfSequences() {
		return alignment.size();
	}
	
	public int getWidth() {
		if (alignment.isEmpty()) {
			System.out.println("Warning: empty seq list in alignment!");
			return 0;
		}
		if (width < 0) {
			width = alignment.get(0).getLen(); 
		}	
		return width;
	}
	
	public int getId() {
		return id;
	}

	public void setId(int id) {
		this.id = id;
	}

	// TODO: check alignment length 
	public void addSeq(String header, String seq) {
		GappedSequence gs = new GappedSequence();
		gs.id = header;
		gs.idxs = createIndexes(seq);
		alignment.add(gs);
		sortByName();
	}

	public void sortByName() {
		java.util.Collections.sort(alignment);
	}
	
	private int[] createIndexes(String seq) {
		int[] idxs = new int[seq.length()];
		int cnt = 0;
		for (int i = 0; i < seq.length(); ++i) {
			if (seq.charAt(i) == '-') {
				idxs[i] = GappedSequence.GAP;
			} else {
				idxs[i] = cnt;
				++cnt;
			}
		}
		return idxs;
	}

	public GappedSequence get(int index) {
		return alignment.get(index);
	}

	public double getScolScore(FastaAlignment testAlignment) throws Exception {
		this.sortByName();
		testAlignment.sortByName();
		int goodCols = 0;
		int allCols = getWidth();
		for (int i = 0; i < allCols; ++i) {
			int[] refCol = getColumn(i);
			int[] otherCol = testAlignment.getCandidateCol(refCol);
			
			if (Arrays.equals(refCol, otherCol)) {
				goodCols ++;
			}
		}
		return 1.0 * goodCols / allCols;
	}

	
	private int[] getCandidateCol(int[] refCol) throws Exception {
		// find first non-gap element of refcol
		int firstNonGapPos = -1;
		for (int i = 0; i < refCol.length ; ++i) {
			if (refCol[i] != GAP) {
				firstNonGapPos = i;
			}
		}
		if (firstNonGapPos == -1) {
			throw new Exception("ERROR, all-gap column in alignment!");
		}
		
		for (int i = 0; i < this.getWidth() ; ++i) {
			int[] candidateCol = getColumn(i);
			if (candidateCol[firstNonGapPos] == refCol[firstNonGapPos]) {
				return candidateCol;
			}
		}
		throw new Exception("ERROR, alignments can not be matched!");
	}

	private int[] getColumn(int idx) {
		int colSize = this.getNumOfSequences();
		int[] col = new int[colSize];
		for (int i = 0; i < colSize; ++i) {
			col[i] = this.get(i).idxs[idx];
		}
		return col;
	}

	public double getFsaScore(FastaAlignment testAlignment) {
		this.sortByName();
		testAlignment.sortByName();
		long t = -System.currentTimeMillis();
		int hom = getHomologNum();
		int nonHom = getNonHomologNum();
		//int totalPairNum = getTotalPairsNum();		
		Pair<Integer, Integer> matchNums = getMatchedNums(testAlignment);
		int goodHom = matchNums.getLeft();
		int goodNonHom = matchNums.getRight();
		return (2.0 * goodHom + goodNonHom) / (2 * hom  + nonHom);
	}

	// first of pair is matched homology num; second is matched non-homology num
	private Pair<Integer, Integer> getMatchedNums(FastaAlignment testAlignment) {
		if (homologyPairs == null) {
			buildHomologyPairs();
		}
		int matchedHom = 0;
		int matchedNonHom = 0;
		for (int i = 0; i < alignment.size(); ++i) {
			int[] first = testAlignment.alignment.get(i).idxs;
			for(int j = i + 1; j < alignment.size(); ++j) {
				int[] second = testAlignment.alignment.get(j).idxs;
				Pair<Integer, Integer> key =
						new ImmutablePair<Integer, Integer>(i, j);
				Set<Pair<Integer, Integer> > refPairs = homologyPairs.get(key);

				for (int k = 0; k < first.length; ++k) {
					// skip if both gaps
					if (first[k] != GAP || second[k] != GAP) {
						Pair<Integer, Integer> seqCol = 
								new ImmutablePair<Integer, Integer>(first[k], second[k]);
						
						if (refPairs.contains(seqCol)) {
							if (first[k] != GAP && second[k] != GAP) {
								matchedHom ++;
							} else {
								matchedNonHom ++;
							}
						}						
					}  
				}
			}
		}
		return new ImmutablePair<Integer, Integer>(matchedHom, matchedNonHom);
	}

	private void buildHomologyPairs() {
		homologyPairs = new HashMap<Pair<Integer, Integer>,
				Set<Pair<Integer, Integer> > >();
		for (int i = 0; i < alignment.size(); ++i) {
			int[] first = alignment.get(i).idxs;
			for(int j = i + 1; j < alignment.size(); ++j) {
				Pair<Integer, Integer> key =
						new ImmutablePair<Integer, Integer>(i, j);
				int[] second = alignment.get(j).idxs;
				Set<Pair<Integer, Integer> > value = collectPairs(first, second);
				homologyPairs.put(key, value);
			}
		}
	}

	private Set<Pair<Integer, Integer> > collectPairs(int[] first, int[] second) {
		Set<Pair<Integer, Integer> > pairs =
				new HashSet<Pair<Integer, Integer> >();
		for (int i = 0; i < first.length; ++i) {
			pairs.add(new ImmutablePair<Integer, Integer>(first[i], second[i]));
		}
		return pairs;
	}

	private int getNonGapNum(int[] col) {
		int ret = 0;
		for (int i = 0; i < col.length; ++i) {
			if (col[i] != GAP) {
				++ret;
			}
		}
		return ret;
	}

	public int getHomologNum() {
		if (homologPairNum < 0) {
			int homologies = 0;
			int allCols = getWidth();
			for (int i = 0; i < allCols; ++i) {
				int[] col = getColumn(i);
				int hom = (getNonGapNum(col) * (getNonGapNum(col) - 1)) / 2;
				homologies += hom;
			}
			homologPairNum = homologies;
		}
		return homologPairNum;
	}
	
	public int getNonHomologNum() {
		if (nonHomologPairNum < 0) {
			int nonhomologies = 0;
			int allCols = getWidth();
			for (int i = 0; i < allCols; ++i) {
				int[] col = getColumn(i);
				int nongaps = getNonGapNum(col);
				int nhom = (nongaps * (col.length - nongaps));
				nonhomologies += nhom;
			}
			nonHomologPairNum = nonhomologies;
		}
		return nonHomologPairNum;
	}
	

	@Override
	public String toString() {
		StringBuilder builder = new StringBuilder();
		for (int i = 0; i < alignment.size(); i++) {
			GappedSequence gs = alignment.get(i); 
			builder.append(gs.id + "\n");
			int[] indexes = gs.idxs;
			for (int j = 0; j < indexes.length; j++) {
				builder.append(indexes[j] + " ");
			}
			builder.append("\n");
		}
		return builder.toString();
	}

	@Override
	public int compareTo(FastaAlignment other) {
		String thisStr = this.toString();
		String otherStr = other.toString();
		return thisStr.compareTo(otherStr);
	}
	
}
