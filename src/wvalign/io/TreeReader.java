package wvalign.io;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Stack;

import wvalign.io.NewickReader.FormatException;
import wvalign.tree.Tree;
import wvalign.tree.TreeNode;


/**
 * Reads a guide tree in Newick format.
 * 
 * @author miklos, novak
 */
public class TreeReader {

	private BufferedReader reader;
	private Tree tree;
	private ArrayList<Tree> trees;
	public Tree getTree() { return tree; }
	public ArrayList<Tree> getTrees() { return trees; }
	
	private double minEdgeLen = 0.01;
	private String newickString;
	
	public TreeReader(String fileName) throws IOException, FileNotFoundException {
		this(new File(fileName));
	}
	
	public TreeReader(File file) throws FileNotFoundException, IOException {
		this(new BufferedReader(new FileReader(file)));
	}

	/**
	 * Constructs default {@link TreeReader} with a minimum edge length of 0.01.
	 */
	public TreeReader(BufferedReader _reader) throws IOException {
		reader = _reader;
		String line;		
		boolean nexus = false;
		while((line = reader.readLine()) != null) {
			if (line.isEmpty()) continue;
			if (line.contains("#NEXUS") || line.contains("#nexus")) {
				nexus = true;
				trees = new ArrayList<Tree>();
				continue;
			}
			String token = line;
			if (nexus) {
				if (line.contains(" tree ")) {
					String[] s = line.split("=");
					token = s[1];
					//token = token.substring(0, token.length()-1); // remove semicolon
				}
				else continue;
			}	
			token = token.replaceAll("\\[\\d+\\]", "");
			//System.out.println(token);
			NewickReader tr = new NewickReader(token, 0);
			try {
				tree = tr.parseTree();				
				if (nexus) trees.add(tree);
			} catch (FormatException e) {
				throw new Error("TreeReader: tree parse exception: "+e);			
			}
		}
	}
	/**
	 * Constructs {@link TreeReader} with a specified minimum edge length.
	 */
	public TreeReader(BufferedReader reader, double minEdgeLen) {
		this.reader = reader;
		this.minEdgeLen = minEdgeLen;
	}
	
	public TreeReader(String newickString, double minEdgeLen) {
		this.minEdgeLen = minEdgeLen;
		this.newickString = newickString;
	}
	
	/**
	 * Reads tree in Newick format. Should be called only once.
	 * 
	 * @return {@link Tree} read
	 * @throws IOException  when I/O error occurs
	 * @throws ParseException  when parsing error occurs
	 */
	public Tree readTree() throws IOException, ParseException {
		TreeNode root = null;
		// we will collect the tree nodes in a stack in order to build a tree
		// input is in newick format (first line of file), example:
		// ((a:0.1,b:0.1):0.1,(c:0.1,d:0.1):0.1);
		// here a,b,c,d are leaf names and every distance (from its parent) is
		// set to 0.1
		Stack<TreeNode> stack = new Stack<TreeNode>();
		String treeStr = makeRootedNewickString();
		boolean isInNumberState = false;
		boolean isLeaf = false;
		String name = "";
		String number = "-1";
		StringBuffer collector = new StringBuffer();
		TreeNode leftNode, rightNode;

		// parsing the tree
		for (int i = 0; i < treeStr.length(); i++) {
			char c = treeStr.charAt(i);
			if (c == '(') { // new non-leaf node
				TreeNode tn = new TreeNode();
				if (stack.empty()) { // the first node
					root = tn; // setting this tree's root
				} else { // setting parent-child links
					TreeNode parent = stack.peek();
					tn.setParent(parent);
					if (!parent.hasLeftChild()) {
						parent.setLeftChild(tn);
					} else {
						parent.setRightChild(tn);
					}
				}
				stack.push(tn);
				collector = new StringBuffer();
				isLeaf = false; // we don' really know it yet but will set this
								// to true if a normal char arrives
			} else if (c == ',') { // end of left child data
				number = collector.toString();
				collector = new StringBuffer();
				isInNumberState = false;
				if (isLeaf) { // end of a leaf
					isLeaf = false; // restore default
					// last seen name is belonging to this leaf
					leftNode = new TreeNode(name); 
					// last seen number
					leftNode.setEvolDist(Math.max(minEdgeLen, new Double(number)));
					TreeNode parent = stack.peek();
					parent.setLeftChild(leftNode);
					// setting parent link
					leftNode.setParent(parent);

				} else { 
					// end of a left-side inner node. It has its children already
					TreeNode tnow = stack.pop();
					tnow.setEvolDist(Math.max(minEdgeLen, new Double(number)));
				}
			} else if (c == ':') { 
				// end of name, or inner node, beginning of numberstate
				isInNumberState = true;
				if (isLeaf) { // end of name (of a leaf node)
					name = collector.toString();
					collector = new StringBuffer();
				} else { // end of children of an inner node
					// for collecting the distance
					collector = new StringBuffer(); 
				}
			} else if (c == ')') { 
				// end of an inner node (except for its evolutionary dist.)
				number = collector.toString();
				collector = new StringBuffer();
				isInNumberState = false;
				if (isLeaf) { // end of a leaf
					// restore default
					isLeaf = false; 
					// last seen name is belonging to this leaf
					rightNode = new TreeNode(name); 
					// last seen number
					rightNode.setEvolDist(Math.max(minEdgeLen, new Double(number))); 
					TreeNode parent = stack.peek();
					parent.setRightChild(rightNode);
					// setting parent link
					rightNode.setParent(parent);
				} else { 
					// End of a right-side inner node. It has its children already.
					TreeNode tnow = stack.pop();
					tnow.setEvolDist(Math.max(minEdgeLen, new Double(number)));
					if (!stack.empty()) {
						TreeNode parent = stack.peek();
						tnow.setParent(parent);
						parent.setRightChild(tnow);
					} else {
						throw new ParseException("Empty stack! ff char: " + treeStr.charAt(i - 1));
					}

				}
			} else if (isInNumberState) { 
				// collecting chars of the evolution distance number (captures negatives)
				if (c >= '0' && c <= '9' || c == '.' || c == '-') {
					collector.append(c);
				} else {
					throw new ParseException("Error while parsing newick tree: found invalid character in \"length\" field!");
				}
			} else {
				// collecting chars of the name of the sequence
				collector.append(c);
				isLeaf = true;
			}
		}
		
		Tree gt = new Tree(root);
		
		return gt;

		// check tree:
//		System.out.println("Printing the tree:");
//		printTree();
	}

	private String makeRootedNewickString() throws IOException {
		if(newickString == null) {
			StringBuilder b = new StringBuilder();
			String s;
			while((s = reader.readLine()) != null)
				b.append(s);
			newickString = b.toString();
		}
		String s = newickString;

		// check if it is rooted;
		// System.out.println(s);
		int par = 0;
		int first;
		// find first ,
		for (first = 0; s.charAt(first) != ',' || par != 1; first++) {
			if (s.charAt(first) == '(') {
				par++;
			} else if (s.charAt(first) == ')') {
				par--;
			}
		}
		// check if we can find a second
		int sec;
		for (sec = first + 1; sec < s.length() && (s.charAt(sec) != ',' || par != 1); sec++) {
			// System.out.println(s.charAt(sec)+" "+par);
			if (s.charAt(sec) == '(') {
				par++;
			} else if (s.charAt(sec) == ')') {
				par--;
			}
		}
		if (sec < s.length()) {
			// root the tree
			s = s.substring(0, first + 1) + "("
					+ s.substring(first + 1, s.length() - 2) + "):0.001);";
		}
		return s;
	}
	
	public static class ParseException extends Exception {
		
		private static final long serialVersionUID = 1L;

		public ParseException(String message) {
			super(message);
		}
		
	}

}