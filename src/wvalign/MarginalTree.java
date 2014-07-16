package wvalign;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Set;

import wvalign.io.TreeReader.ParseException;
import wvalign.model.SubstitutionModel;
import wvalign.tree.Tree;
import wvalign.tree.TreeSplits;
import wvalign.utils.Utils;


public class MarginalTree {
	
	private static final char GAPCHAR = '-';
	
	private Tree tree;
	String outFile;
	AlignmentDAG m;

	public Tree getTree() { return tree; }
	void setTree(Tree _tree) { 
		tree = _tree;
		//tree.indexNodes(tree.getRoot());
		tree.sortNames(m.getSeqNames());
	}
	
	private boolean gapIsChar;		// gap is part of the alphabet (with subst rates etc.)
	
	SubstitutionModel model;
	
	HashMap<Set<List<String> >,ArrayList<Double> > sampledTrees;
	HashMap<String,Boolean> sampledNewicks;
	//HashMap<Double,String> sampledLikelihoods;
	
	long time = 0;
	
	public MarginalTree(AlignmentDAG _m, 
			SubstitutionModel _model,
			String _outfile) {
		
		m = _m;
		model = _model;
		outFile = _outfile;
		checkGap(model);
		System.out.println("Gaps are treated as "+(gapIsChar?"characters":"missing data"));
		time -= System.currentTimeMillis();
	}

	private void checkGap(SubstitutionModel model) {
		gapIsChar = new String(model.alphabet).contains(""+GAPCHAR);
	}

	public double computeLogLikelihood() {

		// annotate
		calcEmissions();
		double fwd = forward(m.columnNetwork.lastCol);
		//System.out.println(fwd);
		
		time += System.currentTimeMillis();
		
		return fwd;
	}
	
	public void sampleTrees(int iterations) throws IOException {
		sampledTrees = new HashMap<Set<List<String> >,ArrayList<Double>>();
		sampledNewicks = new HashMap<String,Boolean>();
	//	sampledTrees = new HashMap<String,Double>();
		tree.indexNodes();
		doubleLink();
		tree.setSubstModel(model);
		//System.out.println("Initial log likelihood = "+computeLogLikelihood());
		double logLikelihood = computeLogLikelihood();
		Set<List<String> > splits = TreeSplits.getNamedSplits(tree);
		String treeString = tree.newickString();
		sampledTrees.put(splits,new ArrayList<Double>(Arrays.asList(logLikelihood)));
		sampledNewicks.put(treeString,true);
		//sampledTrees.put(treeString,logLikelihood);
		try {
			BufferedWriter writer = new BufferedWriter(new FileWriter(outFile));
	
			for (int i=0; i<iterations; i++) {
				// Update tree via random NNI move,
				// keeping branch lengths fixed
				treeString = tree.newickString();
//				System.out.println(treeString);
				tree.NNI();		
				splits = TreeSplits.getNamedSplits(tree);
				if (sampledNewicks.containsKey(treeString)) continue;
				sampledNewicks.put(treeString,true);
				resetForward(m.columnNetwork.lastCol.pred);	
				tree.setSubstModel(model);
				logLikelihood = computeLogLikelihood();
				if (!sampledTrees.containsKey(splits)) {
					sampledTrees.put(splits,new ArrayList<Double>(Arrays.asList(logLikelihood)));
				}
				else {
					sampledTrees.get(splits).add(logLikelihood);
					//System.out.println(sampledTrees.get(splits));
				}
			}
			double marginalLikelihood = Utils.log0;
			for (Set<List<String> > split : sampledTrees.keySet()) {				
				for (double ll : sampledTrees.get(split)) marginalLikelihood = Utils.logAdd(marginalLikelihood,ll);
			}			
			for (Set<List<String> > split : sampledTrees.keySet()) {
				double treeLikelihood = Utils.log0;
				for (double ll : sampledTrees.get(split)) treeLikelihood = Utils.logAdd(treeLikelihood,ll);
				writer.write(split+"\t"+(treeLikelihood-marginalLikelihood)+"\n");
			}
			writer.close();
		} catch (IOException e) {
			System.out.println("Cannot write to file "+outFile);
		}
	}
	
	/**
	 * Calculates emission probabilities for columns and marginals for
	 * column equivalence classes.
	 */
	private void calcEmissions() {
		System.err.print("Calculating emissions...");
		for(ColClass cl : m.columnNetwork.succMap.values()) {
			int n = m.columnNetwork.n;
			double classProb = 0;
			for(Column c : cl.succList)
				classProb += (double)c.count/n;
			cl.score = classProb;
			
			for(Column c : cl.succList) {
				c.scores = new double[1];
				if(c.succ == null) {
					Arrays.fill(c.scores, 0);
					continue;
				}
//				double condMarg = ((double)c.count/n)/classProb;
				char[] obs = getObserv(c.key.desc);
				//System.out.println(obs);
//					c.scores[i] = Math.log(tree.calcSubstLike(obs, i)*condMarg);
				c.scores[0] = Math.log(tree.calcSubstLike(obs, 0));
//					System.out.println("emission "+c+": "+Math.exp(c.scores[0]));
				
			}
		}
		m.columnNetwork.firstCol.scores = new double[1];
		Arrays.fill(m.columnNetwork.firstCol.scores, 0);
		System.err.println("done.");
	}
	
	private void resetForward(ColClass cl) {
		if (cl==null || cl.fwd == null) return;
		cl.fwd = null;
		for(Column c : cl.predList) resetForward(c.pred);
	}
	void resetForward() {
		resetForward(m.columnNetwork.lastCol.pred);
	}
	private double forward(Column lastCol) {
		System.err.print("Computing marginal likelihood...");
		double fwd = forward(lastCol.pred);
		System.err.println("done.");
		return fwd;		
	}

	private double forward(ColClass cl) {
		if(cl == null) {//cl.predList.size() == 1 && cl.predList.get(0).pred == null)
			return 0;
		}
		if(cl.fwd != null)
			return cl.fwd[0];
		double fwd = Utils.log0;
		double fcs;
		for(Column c : cl.predList) {
			double pfwd = forward(c.pred);
//			double jmp = 0;
//			System.out.println(c.count+" "+m.columnNetwork.n+" "+cl.score+" "+c.scores[0]);
			double jmp = Math.log(((double)c.count/m.columnNetwork.n)/cl.score);
			fcs = pfwd+c.scores[0];	// real (log) forward score of column c, state s
			fwd = Utils.logAdd(fwd, fcs+jmp);	
//			System.out.println(c + ": " + jmp + " " + " " + pfwd + " " + fcs + " " +fwd);

		}
		if (cl.fwd == null) cl.fwd = new double[1];
		cl.fwd[0] = fwd;
		return fwd;
	}

	private char[] getObserv(int[] desc) {
		char gap = gapIsChar ? GAPCHAR : 0;
		int i, len = desc.length;
		char[] observ = new char[len];
		for(i = 0; i < len; i++) {
			int x = desc[i];
			observ[i] = (x & 1) == 0 ? gap : Character.toUpperCase(m.sequences[i].charAt(x >> 1));
		}
		return observ;
	}

	/**
	 * Creates predecessor links in the network.
	 */
	void doubleLink() {
		for(ColClass cl : m.columnNetwork.succMap.values()) {
			cl.predList = new ArrayList<Column>();
			for(Column col : cl.succList)
				col.pred = cl;
		}
		for(Column col : m.columnNetwork.contMap.values()) {
			if(col.succ != null)
				col.succ.predList.add(col);
		}
	}

}
