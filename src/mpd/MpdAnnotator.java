package mpd;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

import mpd.tree.Tree;
import statalign.base.Utils;
import statalign.model.subst.SubstitutionModel;

public class MpdAnnotator {
	
	private static final char GAPCHAR = '-';
	
	private Tree tree;
	private double[][] transMat;	// transition matrix of the HMM
	private double[] initState;		// initial state vector of the HMM
	private boolean gapIsChar;		// gap is part of the alphabet (with subst rates etc.)
	
	private ColumnNetwork network;
	
	private int mpdMode = 1;		// 0: P(col), 1: P(col, state)
	private File outFile;
	private String projectSeq;
	private File pseqFile;
	
	long annotTime = 0;

	private double dataProb;
	
	public MpdAnnotator(Tree tree, double[][] transMat, double[] initState, List<? extends SubstitutionModel> models) {
		this.tree = tree;
		this.transMat = Utils.toLog(transMat);
		this.initState = Utils.toLog(initState);
		tree.setSubstModel(models.toArray(new SubstitutionModel[models.size()]));
		checkGap(models.get(0));
	}
	
	public MpdAnnotator(Tree tree, double[][] transMat, double[] initState, SubstitutionModel model, List<Double> rhos) {
		this.tree = tree;
		this.transMat = Utils.toLog(transMat);
		this.initState = Utils.toLog(initState);
		double[] lenMultips = new double[rhos.size()];
		for(int i = 0; i < rhos.size(); i++)
			lenMultips[i] = rhos.get(i);
		tree.setSubstModel(model, lenMultips);
		checkGap(model);
	}

	private void checkGap(SubstitutionModel model) {
		gapIsChar = new String(model.alphabet).contains(""+GAPCHAR);
	}

	public Column annotate(ColumnNetwork network, String[] sequences, String[] seqNames) {
		this.network = network;
		
		System.out.println("Gaps are treated as "+(gapIsChar?"characters":"missing data"));
		System.out.println("MPD mode is "+(mpdMode+1));
		
		annotTime -= System.currentTimeMillis();
		
		// init
		doubleLink();
		tree.sortNames(seqNames);
		
		System.out.println("Calculating emissions...");
		// annotate
		calcEmissions(sequences);
		
		System.out.println("Forward-Backward-Viterbi...");
		double fwd = forward(network.lastCol);
		System.out.println("Forward score: "+fwd);
		dataProb = fwd;
		
		double bwd = backward(network.firstCol);
		System.out.println("Backward score: "+bwd);
		
		double vit = viterbi(network.firstCol, fwd);
		System.out.println("Viterbi score: "+vit);

		annotTime += System.currentTimeMillis();

		try {
			BufferedWriter writer = new BufferedWriter(new FileWriter(outFile));
			
			Column col = network.firstCol.succ.viterbi;
			while(col.succ != null) {
				double total = Utils.log0;
				for(int i = 0; i < initState.length; i++)
					total = Utils.logAdd(total, col.pred.fwd[i]+col.scores[i]+col.succ.bwd[i]);
				for(int state = 0; state < initState.length; state++) {
					double score = Math.exp(col.pred.fwd[state]+col.scores[state]+col.succ.bwd[state]-total);///((double)col.count/network.n);
					if(state > 0)
						writer.write('\t');
					writer.write(""+score);
				}
				writer.write('\n');
				col = col.succ.viterbi;
			}
			writer.close();
			
			
			if(projectSeq != null) {
				int id;
				for(id = 0; id < seqNames.length; id++)
					if(seqNames[id].equals(projectSeq))
						break;
				if(id == seqNames.length)
					throw new IOException("MpdAnnotator: sequence id "+projectSeq+" not found - projected annotation probabilites omitted");
				
				writer = new BufferedWriter(new FileWriter(pseqFile));

				col = network.firstCol.succ.viterbi;
				while(col.succ != null) {
					if((col.key.desc[id] & 1) == 1) {
						double total = Utils.log0;
						for(int i = 0; i < initState.length; i++)
							total = Utils.logAdd(total, col.pred.fwd[i]+col.scores[i]+col.succ.bwd[i]);
						for(int state = 1; state < initState.length; state++) {
							double score = Math.exp(col.pred.fwd[state]+col.scores[state]+col.succ.bwd[state]-total);///((double)col.count/network.n);
							if(state > 1)
								writer.write('\t');
							writer.write(""+score);
						}
						writer.write('\n');
					}
					col = col.succ.viterbi;
				}
				writer.close();
			}

		} catch (IOException e) {
			e.printStackTrace();
		}
		
		return network.firstCol;
	}
	
	/**
	 * Calculates emission probabilities for columns and marginals for
	 * column equivalence classes.
	 */
	private void calcEmissions(String[] sequences) {
		int states = initState.length;
		for(ColClass cl : network.succMap.values()) {
			int n = network.n;
			double classProb = 0;
			for(Column c : cl.succList)
				classProb += (double)c.count/n;
			cl.score = classProb;
			
			for(Column c : cl.succList) {
				c.scores = new double[states];
				if(c.succ == null) {
					Arrays.fill(c.scores, 0);
					continue;
				}
//				double condMarg = ((double)c.count/n)/classProb;
				char[] obs = getObserv(c.key.desc, sequences);
				for(int i = 0; i < states; i++) {
//					c.scores[i] = Math.log(tree.calcSubstLike(obs, i)*condMarg);
					c.scores[i] = Math.log(tree.calcSubstLike(obs, i));
//					System.out.println("emission "+c+" state "+i+": "+Math.exp(c.scores[i]));
				}
			}
		}
		network.firstCol.scores = new double[states];
		Arrays.fill(network.firstCol.scores, 0);
	}
	
	private double forward(Column lastCol) {
		double[] fwd = forward(lastCol.pred);
		double total = Utils.log0;
		for(int i = 0; i < fwd.length; i++)
			total = Utils.logAdd(total, fwd[i]);
		return total;
	}

	private double[] forward(ColClass cl) {
		if(cl == null) {//cl.predList.size() == 1 && cl.predList.get(0).pred == null)
			return initState;
		}
		if(cl.fwd != null)
			return cl.fwd;
		int states = initState.length, s, t;
		double[] fwd = new double[states];
		Arrays.fill(fwd, Utils.log0);
		double fcs;
		for(Column c : cl.predList) {
			double[] pfwd = forward(c.pred);
//			double jmp = 0;
			double jmp = Math.log(((double)c.count/network.n)/cl.score);
			for(s = 0; s < states; s++) {
				fcs = pfwd[s]+c.scores[s];	// real (log) forward score of column c, state s
//				System.out.println("forward "+c+" state "+s+": "+Math.exp(fcs));
				for(t = 0; t < states; t++)
					fwd[t] = Utils.logAdd(fwd[t], fcs+transMat[s][t]+jmp);
			}
		}
		cl.fwd = fwd;
		return fwd;
	}

	private double backward(Column firstCol) {
		double[] bwd = backward(firstCol.succ);
		double total = Utils.log0;
		for(int i = 0; i < bwd.length; i++)
			total = Utils.logAdd(total, initState[i]+bwd[i]);
		return total;
	}

	private double[] backward(ColClass cl) {
		if(cl == null) {
			double[] ret = new double[initState.length];
			Arrays.fill(ret, 0);
			return ret;
		}
		if(cl.bwd != null)
			return cl.bwd;
		int states = initState.length, s, t;
		double[] bwd = new double[states];
		Arrays.fill(bwd, Utils.log0);
		double bcs;
		for(Column c : cl.succList) {
			double[] sbwd = backward(c.succ);
//			double jmp = 0;
			double jmp = Math.log(((double)c.count/network.n)/cl.score);
			for(t = 0; t < states; t++) {
				bcs = sbwd[t]+c.scores[t];
//				System.out.println("backward "+c+" state "+t+": "+Math.exp(sbwd[t]));
				for(s = 0; s < states; s++)
					bwd[s] = Utils.logAdd(bwd[s], bcs+transMat[s][t]+jmp);
			}
		}
		cl.bwd = bwd;
		return bwd;
	}

	private double viterbi(Column firstCol, double dataProb) {
		for(ColClass succClass : network.succMap.values()) {
			succClass.score = Double.NEGATIVE_INFINITY;
			succClass.viterbi = null;
		}
		return viterbi(firstCol.succ, dataProb);
	}

	private double viterbi(ColClass colClass, double dataProb) {
		if(colClass == null)
			return 0;
		if(colClass.score != Double.NEGATIVE_INFINITY)
			return colClass.score;
		double[] fwd = colClass.fwd, bwd;
		for(Column col : colClass.succList) {
			if(col.succ == null) {
				colClass.viterbi = col;
				return 0;
			}
			bwd = col.succ.bwd;
			double total = Utils.log0, score;
			if(mpdMode == 1)
				for(int i = 0; i < initState.length; i++)
					total = Utils.logAdd(total, fwd[i]+col.scores[i]+bwd[i]);
			for(int i = 0; i < initState.length; i++) {
				if(mpdMode == 0)
					score = ((double)col.count/network.n) - network.gValue + viterbi(col.succ, dataProb);
				else
					score = Math.exp(fwd[i]+col.scores[i]+bwd[i]-total) * ((double)col.count/network.n) - network.gValue + viterbi(col.succ, dataProb);
				if(score > colClass.score) {
					colClass.score = score;
					colClass.viterbi = col;
					colClass.vitState = i;
				}
			}
		}
		return colClass.score;
	}

	private char[] getObserv(int[] desc, String[] sequences) {
		char gap = gapIsChar ? GAPCHAR : 0;
		int i, len = desc.length;
		char[] observ = new char[len];
		for(i = 0; i < len; i++) {
			int x = desc[i];
			observ[i] = (x & 1) == 0 ? gap : Character.toUpperCase(sequences[i].charAt(x >> 1));
		}
		return observ;
	}

	/**
	 * Creates predecessor links in the network.
	 */
	private void doubleLink() {
		for(ColClass cl : network.succMap.values()) {
			cl.predList = new ArrayList<Column>();
			for(Column col : cl.succList)
				col.pred = cl;
		}
		for(Column col : network.contMap.values()) {
			if(col.succ != null)
				col.succ.predList.add(col);
		}
	}

	/**
	 * Sets output file for annotation.
	 */
	public void setOutFile(File outFile) {
		this.outFile = outFile;
	}
	
	/**
	 * Sets MPD mode.
	 * @param mpdMode 0 for P(col) and 1 for P(col,state)
	 */
	public void setMpdMode(int mpdMode) {
		this.mpdMode = mpdMode;
	}
	
	/**
	 * Selects sequence to project annotation probabilites onto and
	 * file where this is printed.
	 * @param projectSeq a sequence id
	 * @param pseqFile output file
	 */
	public void setProjectSeq(String projectSeq, File pseqFile) {
		this.projectSeq = projectSeq;
		this.pseqFile = pseqFile;
	}
	
	public double getDataProb() {
		return dataProb;
	}
}
