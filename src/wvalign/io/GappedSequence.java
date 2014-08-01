package wvalign.io;

public class GappedSequence implements Comparable<GappedSequence> {

	final static int GAP = -1;
	public String id;
	public int[] idxs;
	
	@Override
	public int compareTo(GappedSequence other) {
		return id.compareTo(other.id);		
	}
	
	public int getLen() {
		return idxs.length;
	}
	
}
