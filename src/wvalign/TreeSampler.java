package wvalign;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.List;

import wvalign.model.SubstitutionModel;
import wvalign.tree.Tree;
import wvalign.utils.Utils;


public class TreeSampler {
	
	private static final char GAPCHAR = '-';
	
	private Tree tree;
	String outFile;
	MinRiskGlobalFast m;

	public Tree getTree() { return tree; }
	void setTree(Tree _tree) { tree = _tree; }
	
	private boolean gapIsChar;		// gap is part of the alphabet (with subst rates etc.)
	
	SubstitutionModel model;
	
	HashMap<String,Double> sampledTrees;
	
	long time = 0;
	
	public TreeSampler(MinRiskGlobalFast m, 
			SubstitutionModel model,
			String outfile) {
		
		this.m = m;
		this.model = model;
		this.outFile = outfile;
		checkGap(model);
		System.out.println("Gaps are treated as "+(gapIsChar?"characters":"missing data"));
		time -= System.currentTimeMillis();
	}

	private void checkGap(SubstitutionModel model) {
		gapIsChar = new String(model.alphabet).contains(""+GAPCHAR);
	}

	public double computeLogLikelihood() {

		System.out.println("Calculating emissions...");
		// annotate
		calcEmissions();
		
		System.out.println("Computing marginal likelihood...");
		double fwd = forward(m.network.lastCol);
		System.out.println(fwd);
		
		time += System.currentTimeMillis();
		
		return fwd;
	}
	
	public void sampleTrees(int iterations) {
		sampledTrees = new HashMap<String,Double>();
		tree.indexNodes();
		doubleLink();
		try {
			BufferedWriter writer = new BufferedWriter(new FileWriter(outFile));
	
			for (int i=0; i<iterations; i++) {
				// Update tree via random NNI move,
				// keeping branch lengths fixed
				tree.NNI();		
				String treeString = tree.newickString();
				if (!sampledTrees.containsKey(treeString)) {
					tree.setSubstModel(model);
					sampledTrees.put(treeString, computeLogLikelihood());
				}
			}
			for (String t : sampledTrees.keySet()) {
				writer.write(t+"\t"+sampledTrees.get(t)+"\n");
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
		for(ColClass cl : m.network.succMap.values()) {
			int n = m.network.n;
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
//					c.scores[i] = Math.log(tree.calcSubstLike(obs, i)*condMarg);
				c.scores[0] = Math.log(tree.calcSubstLike(obs, 0));
//					System.out.println("emission "+c+" state "+i+": "+Math.exp(c.scores[i]));
				
			}
		}
		m.network.firstCol.scores = new double[1];
		Arrays.fill(m.network.firstCol.scores, 0);
	}
	
	private double forward(Column lastCol) {
		double fwd = forward(lastCol.pred);
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
			double jmp = Math.log(((double)c.count/m.network.n)/cl.score);
			fcs = pfwd+c.scores[0];	// real (log) forward score of column c, state s
			fwd = Utils.logAdd(fwd, fcs+jmp);		
		}
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
	private void doubleLink() {
		for(ColClass cl : m.network.succMap.values()) {
			cl.predList = new ArrayList<Column>();
			for(Column col : cl.succList)
				col.pred = cl;
		}
		for(Column col : m.network.contMap.values()) {
			if(col.succ != null)
				col.succ.predList.add(col);
		}
	}

}
