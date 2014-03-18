package wvalign;

import java.util.ArrayList;
import java.util.List;

class ColClass {
	List<Column> succList = new ArrayList<Column>();
	double score;
	Column viterbi;
	// annotation
	List<Column> predList;
	double[] fwd, bwd;
	int vitState;
	int predFreq = 0;
	int succFreq = 0;
}