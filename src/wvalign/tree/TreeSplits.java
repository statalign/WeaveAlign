package wvalign.tree;

import java.io.IOException;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;

import wvalign.io.TreeReader;
import wvalign.io.TreeReader.ParseException;


public class TreeSplits {
	
	public static Set<List<Integer>> getSplits(String treeFile) throws IOException, ParseException {
		Tree tree = new TreeReader(treeFile).getTree();
		//TreeNode root = new TreeReader(treeFile).readTree().getRoot();
		return getSplits(tree);
	}
	public static Set<List<Integer>> getSplits(Tree tree) throws IOException {
		TreeNode root = tree.getRoot();
		List<String> leafNames = new ArrayList<String>();
		root.recursiveCollectNodeNames(leafNames);
		Collections.sort(leafNames);
		//System.out.println(leafNames);
		Map<String, Integer> nameMap = new HashMap<String, Integer>();
		for(int i = 0; i < leafNames.size(); i++)
			nameMap.put(leafNames.get(i), i);
		
		root.findNodesBelow(nameMap);
		
		Set<List<Integer>> splits = new HashSet<List<Integer>>();
		if(!root.isLeaf()) {
			root.getLeftChild().collectSplits(root.getRightChild().getNodesBelow(), splits);
			root.getRightChild().collectSplits(root.getLeftChild().getNodesBelow(), splits);
		}
//		for (List<Integer> split : splits) {
//			for (int i : split) System.out.print(leafNames.get(i)+" ");
//			System.out.println();
//		}
		return splits;
	}
	public static Set<List<String>> getNamedSplits(Tree tree) throws IOException {
		Set<List<Integer>> splits = getSplits(tree);
		TreeNode root = tree.getRoot();
		List<String> leafNames = new ArrayList<String>();
		root.recursiveCollectNodeNames(leafNames);
		Collections.sort(leafNames);
		
		Set<List<String>> namedSplits = new HashSet<List<String>>();
		for (List<Integer> split : splits) {
			List<String> names = new ArrayList<String>();
			for (int node : split) names.add(leafNames.get(node));
			namedSplits.add(names);
		}
		return namedSplits;
	}
	
	public static void main(String[] args) {
		try {
			//System.out.println(getSplits("eval/uncertain/1aboA.fasta.dnd"));
			System.out.println(getSplits("testdata/tree.nwk"));
		} catch (Exception e) {
			e.printStackTrace();
		}
	}
}
