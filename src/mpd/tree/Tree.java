package mpd.tree;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.LinkedHashMap;
import java.util.List;

import mpd.model.SubstitutionModel;
import mpd.utils.Utils;


/**
 * This class represents a guide tree for the multiple sequence alignment.
 */
public class Tree {

	TreeNode root;
	
	SubstitutionModel substModel;
	boolean unknownCharWarn = true;
	
	public Tree(TreeNode root) {
		this.root = root;
	}

	public int numberOfLeaves() {
		return root.numberOfLeaves();
	}

	/**
	 * Finds a node in the tree by a given name.
	 * 
	 * @param name
	 *            node name to search for
	 * @return <code>TreeNode</code> or null if not found
	 */
	public TreeNode findNodeByName(String name) {
		return root.traverseFindDeep(name);
	}

	/** @return the root */
	public TreeNode getRoot() {
		return root;
	}

	/** @param root the root to set */
	public void setRoot(TreeNode root) {
		this.root = root;
	}

	public String newickString() {
		return root.newickString() + ";";
	}

	public List<String> getNames() {
		List<String> names = new ArrayList<String>();
		root.recursiveCollectNodeNames(names);
		return names;
	}
	
	/**
	 * Provides the reference ordering of the sequences.
	 * Also checks whether tree is compatible with sequences and throws {@link Error} if not.
	 * @param seqNames names of sequences in the same order as character observations in subsequent
	 * 		calls to {@link #calcSubstLike(char[], int)}
	 */
	public void sortNames(String[] seqNames) {
		int leaves = numberOfLeaves();
		if(leaves != seqNames.length)
			throw new Error("Tree is not compatible with sequences (#leaves="+leaves+" , #seqs="+seqNames.length+")!");
		HashMap<String, Integer> nameLookup = new HashMap<String, Integer>();
		for(int i = 0; i < seqNames.length; i++) {
			String name = seqNames[i];
			if(nameLookup.get(name) != null)
				throw new Error("Duplicate name "+name+" in tree!");
			nameLookup.put(name, i);
		}
		root.setNodeIds(nameLookup);
	}
	
	double[] gapFelsen;
	
	/**
	 * Sets the substitution model for the tree and precalculates exponentials along edges
	 * using a list of edge length multipliers (model indices).
	 */
	public void setSubstModel(SubstitutionModel substModel, double[] lenMultip) {
		this.substModel = substModel;
		root.precalcSubstMats(substModel, lenMultip);
		init(lenMultip.length);
	}
	
	/**
	 * Sets the default substitution model for the tree that is used to translate characters
	 * as the first one provided and precalculates exponentials along edges from the list
	 * of rate matrices.
	 */
	public void setSubstModel(SubstitutionModel[] substModels) {
		substModel = substModels[0];
		root.precalcSubstMats(substModels);
		init(substModels.length);
	}
	
	@SuppressWarnings("unchecked")
	private void init(int states) {
		gapFelsen = new double[substModel.e.length];
		Arrays.fill(gapFelsen, 1.0);

		caches = new LinkedHashMap[states];
		for(int i = 0; i < caches.length; i++) {
			caches[i] = new LinkedHashMap<String, Double>() {
				private static final long serialVersionUID = 1L;
				
				@Override
				protected boolean removeEldestEntry(java.util.Map.Entry<String,Double> eldest) {
					return size() > 100000;
				};
			};
		}
	}

	LinkedHashMap<String, Double>[] caches;
	
	/**
	 * Calculates likelihood of observation vector given the tree and the model index.
	 * The method {@link #sortNames(ArrayList)} must be called before this method.
	 * @param observation character vector (zero for gaps)
	 * @param modelInd index of substitution model relative to the list of models specified in
	 * 		{@link #setSubstModel(SubstitutionModel, double[])}) or
	 * 		{@link #setSubstModel(SubstitutionModel[])}
	 * @return
	 */
	public double calcSubstLike(char[] observation, int modelInd) {
		String obs = new String(observation);
		Double d = caches[modelInd].get(obs);
		if(d != null)
			return d;
		int[][] which = substModel.attachedScoringScheme.which;
		double[][] felsenObserv = new double[observation.length][substModel.e.length];
		for(int i = 0; i < observation.length; i++) {
			char ch = observation[i];
			if(ch == 0) {
				felsenObserv[i] = gapFelsen;
			} else {
				int[] subWh = which[ch];
				if(subWh == null) {
					if(unknownCharWarn)
						System.err.println("Warning: unknown character '"+ch+"' - treating as missing data");
					felsenObserv[i] = gapFelsen;
					continue;
				}
				double[] subFel = felsenObserv[i];
				for(int j = 0; j < subWh.length; j++)
					subFel[j] = subWh[j];
			}
		}
		double[] felsen = root.calcSubstLike(felsenObserv, modelInd);
		double v = Utils.calcEmProb(felsen, substModel.e);
		caches[modelInd].put(obs, v);
		return v;
	}

}
