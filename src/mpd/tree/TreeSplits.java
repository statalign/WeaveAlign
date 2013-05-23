package mpd.tree;

import java.io.IOException;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import mpd.io.TreeReader;
import mpd.io.TreeReader.ParseException;

public class TreeSplits {
	
	public static List<List<Integer>> getSplits(String treeFile) throws IOException, ParseException {
		TreeNode root = new TreeReader(treeFile).readTree().getRoot();
		
		List<String> leafNames = new ArrayList<String>();
		root.recursiveCollectNodeNames(leafNames);
		Collections.sort(leafNames);
		System.out.println(leafNames);
		Map<String, Integer> nameMap = new HashMap<String, Integer>();
		for(int i = 0; i < leafNames.size(); i++)
			nameMap.put(leafNames.get(i), i);
		
		root.findNodesBelow(nameMap);
		
		List<List<Integer>> splits = new ArrayList<List<Integer>>();
		if(!root.isLeaf()) {
			root.getLeftChild().collectSplits(root.getRightChild().getNodesBelow(), splits);
			root.getRightChild().collectSplits(root.getLeftChild().getNodesBelow(), splits);
		}
		return splits;
	}
	
	public static void main(String[] args) {
		try {
			System.out.println(getSplits("eval/uncertain/1aboA.dnd"));
		} catch (Exception e) {
			e.printStackTrace();
		}
	}
}
