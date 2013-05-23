package mpd.io;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;

import mpd.tree.Tree;
import mpd.tree.TreeNode;
import statalign.base.CircularArray;

/**
 * Class to parse Newick format tree files or strings.
 *
 * The tree must either be an rooted or unrooted binary tree, the latter will be rooted after
 * parsing at the middle of the third branch of the top node as the Tree class does not support
 * unrooted trees.
 *
 * The vertex[] array of the tree are filled such that the leaf nodes come first in the order
 * they are visited by a left-first preorder traversal then the internal nodes in the reversed
 * order they are visited by the same traversal (root will be the last). The latter ones get
 * labelled by their index in the array. Vertex.vertexNum fields are filled.
 *
 * You can specify a lower limit on edge lengths using the appropriate constructors. Default
 * is a lower limit of 0.
 *
 * Any number of white space (including new lines) may be present.
 *
 * @author novak
 *
 */
public class NewickReader {

	private BufferedReader reader;

	// internal buffer holding the String representation of the complete tree
	private String buf;
	// current parsing position within buf.
	private int pos;

	CircularArray<TreeNode> leaves;
	CircularArray<TreeNode> inNodes;

	private double minEdgeLen = 0.0;

	/**
	 * Constructs a NewickReader from a file with default minimum edge length of 0.
	 * @throws FileNotFoundException if file not found
	 */
	public NewickReader(File file) throws FileNotFoundException {
		reader = new BufferedReader(new FileReader(file));
	}

	/**
	 * Constructs a NewickReader from a file with specified minimum edge length.
	 * @throws FileNotFoundException if file not found
	 */
	public NewickReader(File file, double minEdgeLen) throws FileNotFoundException {
		this(file);
		this.minEdgeLen = minEdgeLen;
	}

	/**
	 * Constructs a NewickReader from a String with default minimum edge length of 0.
	 */
	public NewickReader(String tree) {
		buf = tree;
	}

	/**
	 * Constructs a NewickReader from a String with specified minimum edge length.
	 */
	public NewickReader(String tree, double minEdgeLen) {
		buf = tree;
		this.minEdgeLen = minEdgeLen;
	}

	/**
	 * Parses tree. Should be called once for each NewickReader object.
	 *
	 * @return tree built
	 * @throws IOException only if parsing from file and I/O error occurs
	 * @throws FormatException
	 */
	public Tree parseTree() throws IOException, FormatException {
		if(reader != null) {
			readBuf();
		}

		leaves = new CircularArray<TreeNode>();
		inNodes = new CircularArray<TreeNode>();

		TreeNode root = new TreeNode();

		inNodes.push(root);

		expectNext("(");
		parseSubtree(root, 0, true);
		if(expectNext(",)") == 1) {			// handle unrooted tree
			parseSubtree(root, 1, false);
			root.getLeftChild().setEvolDist(Math.max(minEdgeLen, root.getLeftChild().getEvolDist()/2));
			root.getRightChild().setEvolDist(root.getLeftChild().getEvolDist());
		} else {
			parseSubtree(root, 1, true);
			if(expectNext(",)") == 0) {			// handle unrooted tree and transform to rooted one
				root.setParent(new TreeNode());
				root.getParent().setLeftChild(root);
				root = root.getParent();
				inNodes.unshift(root);		// put new root to the beginning of inNodes
				parseSubtree(root, 1, true);
				if(root.getRightChild().getEvolDist() < 2*minEdgeLen) {
					root.getLeftChild().setEvolDist(minEdgeLen);
					root.getRightChild().setEvolDist(minEdgeLen);
				} else {
					root.getLeftChild().setEvolDist(root.getRightChild().getEvolDist()/2);
					root.getRightChild().setEvolDist(root.getLeftChild().getEvolDist());
				}
				expectNext(")");
			}
		}
		expectNext(";");

		// name internal nodes by numbers starting from the number of leaves
		int id = leaves.length();
		TreeNode vertex;
		while((vertex = inNodes.pop()) != null) {
			if(vertex.getLeftChild() == null)		// has been deleted
				continue;
			vertex.setName(Integer.toString(id++));
			leaves.push(vertex);
		}

		Tree tree = new Tree(root);
		return tree;
	}

	/**
	 * Recursively parses a subtree with edge length to parent from buf starting at pos.
	 * Moves pos to first character after the subtree's representation.
	 * Subtree can also be a single named node. Root of the built subtree will be set
	 * as the child of parent.
	 * @param parent parent of the subtree
	 * @param child 0 if subtree should be left child 1 if right
	 * @throws FormatException
	 */
	protected void parseSubtree(TreeNode parent, int child, boolean withLen) throws FormatException {
		TreeNode v = new TreeNode();
		v.setParent(parent);
		if(child == 0) {
			parent.setLeftChild(v);
		} else {
			parent.setRightChild(v);
		}

		char x = next();
		if(x != '(') {		// leaf node
			leaves.push(v);
			// get name
			int found = expectLater("():,;");
			v.setName(buf.substring(pos, found).trim());
			pos = found;		// skip name
		} else {		// internal node
			pos++;				// skip '('
			inNodes.push(v);
			parseSubtree(v, 0, withLen);
			if(expectNext(",)") == 1) {		// no right subtree: delete node
				TreeNode w = v.getLeftChild();
				w.setParent(parent);
				if(child == 0)
					parent.setLeftChild(w);
				else
					parent.setRightChild(w);
				v.setLeftChild(null);		// delete from inNodes
				v = w;				// add edge length to existing node
			} else {
				parseSubtree(v, 1, withLen);
				expectNext(")");
			}
		}

		if(withLen) {
			// get edge length
			expectNext(":");
			int found = expectLater("():,");
			try {
				v.setEvolDist(Math.max(v.getEvolDist()+Double.parseDouble(buf.substring(pos, found)), minEdgeLen));
			} catch (NumberFormatException e) {
				throw new FormatException(FormatExceptType.EDGELEN_ERROR);
			}
			pos = found;
		}
	}

	/**
	 * Reads tree file and stores it in buf.
	 * @throws IOException
	 */
	protected void readBuf() throws IOException {
		StringBuilder build = new StringBuilder();
		String str;
		while((str=reader.readLine()) != null) {
			build.append(str);
		}
		reader.close();
		buf = build.toString();
	}

	/**
	 * Returns the index of the character in chars that is the next non-whitespace at pos
	 * (pos is then incremented)
	 * @throws FormatException if next non-whitespace at pos isn't one of chars
	 */
	protected int expectNext(String chars) throws FormatException {
		char x = next();
		int found = chars.indexOf(x);
		if(found == -1) {
			throw new FormatException(FormatExceptType.PARSE_ERROR, chars);
		}
		pos++;
		return found;
	}

	/**
	 * Returns the first index (starting from pos) where any of the characters in chars is found,
	 * pos is not moved.
	 */
	protected int expectLater(String chars) throws FormatException {
		for(int i = pos; i < buf.length(); i++) {
			if(chars.indexOf(buf.charAt(i)) != -1) {
				if(i == pos) {
					throw new FormatException(FormatExceptType.EMPTY_FIELD);
				}
				return i;
			}
		}
		throw new FormatException(FormatExceptType.PREMATURE_END, chars);
	}

	/**
	 * Returns next non-whitespace character and moves pos to its position.
	 */
	protected char next() throws FormatException {
		for(;;) {
			if(pos >= buf.length()) {
				throw new FormatException(FormatExceptType.PREMATURE_END);
			}
			char x = buf.charAt(pos);
			if(!Character.isWhitespace(x)) {
				return x;
			}
			pos++;
		}
	}

	public class FormatException extends Exception {
		private static final long serialVersionUID = 1L;

		FormatExceptType type;
		String expected;

		public FormatException(FormatExceptType type) {
			this.type = type;
		}

		public FormatException(FormatExceptType type,
				String expected) {
			this.type = type;
			this.expected = expected;
		}

		@Override
		public String getMessage() {
			return type.getMessage()+(expected != null ? " - expected "+
					(expected.length()>1?" one of `":" `")+expected+"'":"")+
					" at "+(pos+1)+" near `"+buf.substring(pos)+"'";
		}

	}

	public static enum FormatExceptType {

		PARSE_ERROR("Parse error"),
		PREMATURE_END("Premature end of input"),
		EDGELEN_ERROR("Ill-formatted edge length"),
		EMPTY_FIELD("Edge length or node name field has zero length");

		private String message;

		private FormatExceptType(String message) {
			this.message = message;
		}

		public String getMessage() {
			return message;
		}

	}

	public static void main(String[] args) throws IOException, FormatException {
		Tree tree = new NewickReader("  (x:1e1,(t:3,u:3):3,(x:9,y:9):1);").parseTree();
		System.out.println(tree.newickString());
	}
}
