package wvalign.tree;

import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import wvalign.model.SubstitutionModel;
import wvalign.utils.Utils;




/**
 * Class that represents a single node of the guide tree with reference
 * to ancestor and descendant nodes.
 */
public class TreeNode {

	private String name;
	private TreeNode leftChild;
	private TreeNode rightChild;
	private TreeNode parent;
	
	/** Evolutionary distance from its parent node. Read from the Newick file
	 * of the guide tree or computed (with NJ/UPGMA) */
	private double evolDist; // distance from its parent

	int nodeId = -1;		// id of this node (based on rank in the alphabetical order of names)
	double[][][] substMat;	// precalculated character transition likelihoods for edge length (a matrix per edge multiplier)

	/** Number of leaves under this node */
	private int numOfLeaves = -1;
	
	private List<Integer> nodesBelow;
	
	public TreeNode() {
		name = "";		
	}

	public TreeNode(String name) {
		this.name = name;
	}

	public TreeNode(TreeNode left, TreeNode right) {
		leftChild = left;
		rightChild = right;
		name = left.getName() + " - " + right.getName();
	}

	/** @param leftChild the left child to set */
	public void setLeftChild(TreeNode leftChild) {
		this.leftChild = leftChild;
	}

	/** @param rightChild the right child to set */
	public void setRightChild(TreeNode rightChild) {
		this.rightChild = rightChild;
	}

	/** @return the left child */
	public TreeNode getLeftChild() {
		return leftChild;
	}

	/** @return the right child */
	public TreeNode getRightChild() {
		return rightChild;
	}

	/** @return the parent */
	public TreeNode getParent() {
		return parent;
	}

	/** @param parent the parent to set */
	public void setParent(TreeNode parent) {
		this.parent = parent;
	}

	public int numberOfLeaves() {
		return leftChild == null ? 1 : leftChild.numberOfLeaves()
				+ rightChild.numberOfLeaves();
	}

	/** @return the name of the sequence */
	public String getName() {
		return name;
	}

	/**
	 * @param name
	 *            the name to set
	 */
	public void setName(String name) {
		this.name = name;
	}

	/** @return the evolDist */
	public double getEvolDist() {
		return evolDist;
	}

	/**
	 * @param evolDist
	 *            the evolutionary distance (read from the Newick file of the
	 *            guide tree) to set
	 */
	public void setEvolDist(double evolDist) {
		this.evolDist = evolDist;
	}

	/** @return <tt>true</tt> if this node is a leaf in the tree (has no children). */
	public boolean isLeaf() {
		return leftChild == null && rightChild == null;
	}

	public boolean hasLeftChild() {
		return leftChild != null;
	}

	public boolean hasRightChild() {
		return rightChild != null;
	}

	public String newickString() {
		if (rightChild != null && leftChild != null)
			return "(" + leftChild.newickString() + ","
					+ rightChild.newickString() + ")"
					+ (parent == null ? "" : (":" + getEvolDist()));
		return name + ":" + getEvolDist();
	}

	/**
	 * Finds and returns node in the subtree rooted in this node with the given
	 * name.
	 * 
	 * @param name
	 *            node name to search for
	 * @return <code>TreeNode</code> or null if not found
	 */
	public TreeNode traverseFindDeep(String name) {
		if (this.name.equals(name))
			return this;
		if (isLeaf())
			return null;
		TreeNode l = leftChild.traverseFindDeep(name);
		if (l != null)
			return l;
		return rightChild.traverseFindDeep(name);
	}

	public void recursiveCollectNodeNames(List<String> names) {
		if (leftChild == null && rightChild == null) {
			names.add(getName());
			return;
		}
		if (leftChild != null)
			leftChild.recursiveCollectNodeNames(names);
		if (rightChild != null)
			rightChild.recursiveCollectNodeNames(names);
	}

	public int resetNumOfLeaves() {
		numOfLeaves = -1;
		return getNumOfLeaves(); 
	}
	
	public int getNumOfLeaves() {
		if (numOfLeaves < 0) {
			if (!hasLeftChild() && !hasRightChild()) {
				numOfLeaves = 1;
			} else if (!hasLeftChild()) {
				System.err.println("TreeNode.java: unbalanced tree!");
				numOfLeaves = rightChild.getNumOfLeaves();
			} else { // has 2 children
				numOfLeaves = leftChild.getNumOfLeaves() + rightChild.getNumOfLeaves();
			}
		}
		return numOfLeaves;
	}
	
	public List<Integer> getNodesBelow() {
		return nodesBelow;
	}
	
	public List<Integer> findNodesBelow(Map<String,Integer> nameMap) {
		nodesBelow = new ArrayList<Integer>();
		if(isLeaf()) {
			nodesBelow.add(nameMap.get(name));
		} else {
			nodesBelow.addAll(leftChild.findNodesBelow(nameMap));
			nodesBelow.addAll(rightChild.findNodesBelow(nameMap));
		}
		Collections.sort(nodesBelow);
		return nodesBelow;
	}

	public void collectSplits(List<Integer> nodesAbove, List<List<Integer>> splits) {
		if(isLeaf())
			return;
		
		splits.add(nodesBelow);
		
		List<Integer> split = new ArrayList<Integer>();
		split.addAll(nodesAbove);
		split.addAll(rightChild.nodesBelow);
		Collections.sort(split);
		splits.add(split);
		leftChild.collectSplits(split, splits);
		
		split = new ArrayList<Integer>();
		split.addAll(nodesAbove);
		split.addAll(leftChild.nodesBelow);
		Collections.sort(split);
		splits.add(split);
		rightChild.collectSplits(split, splits);
	}

	/**
	 * Recursively precalculates exponentials of the substitution rate matrix along edges
	 * using different edge length multipliers.
	 */
	public void precalcSubstMats(SubstitutionModel substModel, double[] lenMultip) {
		substMat = new double[lenMultip.length][][];
		for(int i = 0; i < lenMultip.length; i++)
			substMat[i] = substModel.updateTransitionMatrix(null, evolDist*lenMultip[i]);
		if(leftChild != null)
			leftChild.precalcSubstMats(substModel, lenMultip);
		if(rightChild != null)
			rightChild.precalcSubstMats(substModel, lenMultip);
	}

	/**
	 * Recursively precalculates exponentials of a list of substitution rate matrices along edges.
	 */
	public void precalcSubstMats(SubstitutionModel[] substModels) {
		substMat = new double[substModels.length][][];
		for(int i = 0; i < substModels.length; i++)
			substMat[i] = substModels[i].updateTransitionMatrix(null, evolDist);
		if(leftChild != null)
			leftChild.precalcSubstMats(substModels);
		if(rightChild != null)
			rightChild.precalcSubstMats(substModels);
	}
	public void precalcSubstMats(SubstitutionModel substModel) {
		substMat = new double[1][][];		
		substMat[0] = substModel.updateTransitionMatrix(null, evolDist);
		if(leftChild != null)
			leftChild.precalcSubstMats(substModel);
		if(rightChild != null)
			rightChild.precalcSubstMats(substModel);
	}
	
	public double[] calcSubstLike(double[][] observation, int modelInd) {
		if(leftChild == null || rightChild == null)
			return observation[nodeId];
		double[] left = leftChild.calcSubstLike(observation, modelInd);
		double[] right = rightChild.calcSubstLike(observation, modelInd);
//		if(leftChild.leftChild == null && rightChild.leftChild == null && leftChild.name.equals("dp4")) {
//			String lc = left[4] == 0 ? "ng" : "g";
//			String rc = right[4] == 0 ? "ng" : "g";
//			System.out.println(lc+" to "+rc+" in "+(leftChild.evolDist+rightChild.evolDist));
//		}
		double[] res = new double[left.length];
		Utils.calcFelsen(res, left, leftChild.substMat[modelInd], right, rightChild.substMat[modelInd]);
		return res;
	}

	public void setNodeIds(HashMap<String, Integer> nameLookup) {
		if(leftChild == null || rightChild == null) {
			Integer id = nameLookup.get(name);
			if(id == null)
				throw new Error("Sequence "+name+" not found among sequences");
			nodeId = id;
		} else {
			leftChild.setNodeIds(nameLookup);
			rightChild.setNodeIds(nameLookup);
		}
	}
}