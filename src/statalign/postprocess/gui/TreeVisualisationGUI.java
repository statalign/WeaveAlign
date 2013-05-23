package statalign.postprocess.gui;

import java.awt.Color;
import java.awt.Dimension;
import java.awt.Font;
import java.awt.Graphics;
import java.awt.Graphics2D;

import javax.swing.JPanel;

import statalign.base.Vertex;
import statalign.postprocess.plugins.TreeVisualisation;

/**
 * The graphical interface for showing the current tree.
 *
 * @author miklos, novak
 *
 */
public class TreeVisualisationGUI extends JPanel{
	/**
	 *
	 */
	private static final long serialVersionUID = 1L;
	JPanel panel;
	TreeVisualisation owner;

	/**
	 *
	 * @param panel The main panel.
	 * @param owner The TreeVisualisation postprocess handler.
	 */
	public TreeVisualisationGUI(JPanel panel, TreeVisualisation owner){
		this.panel = panel;
		this.owner = owner;
	}

	/**
	 * It repaints the graphics
	 */
	@Override
	public void paintComponent(Graphics gr){
		//System.out.println("Updating tree");

		//String newick = owner.mcmc.tree.printedTree();
		Graphics2D g = (Graphics2D)gr.create();

		g.setBackground(Color.WHITE);

	    g.setFont(new Font("MONOSPACED", Font.PLAIN, 12));
	    g.setClip(0,0, panel.getWidth(), panel.getHeight());

		g.clearRect(0, 0, this.getWidth(), this.getHeight());
	    g.setColor(Color.BLACK);
		//g.drawString(newick,100,100);
		Vertex root = owner.mcmc.tree.root;
		double maxDepth = root.maxDepth() * 1.15;
		int nameLengthCount = root.countNameLengthSum();
		drawTree(g,
				root,
				maxDepth,
				nameLengthCount,
				0, this.getWidth(),
				(int)(this.getHeight()*0.05));
	}

	private void drawTree(Graphics g, Vertex v, double maxDepth, int countTotalLength, int x1, int x2, int y){
		if(v.left == null){
			String s = (v.name.length() > 10 ? v.name.substring(0, 10)+"..." : v.name);
			g.drawString(s, (x1+x2)/2-(s.indexOf(' ') == -1 ? s.length() : s.indexOf(' ')) * 4, y+15);
			g.drawLine((x1+x2)/2, y, (x1+x2)/2, y - (int)(v.edgeLength * this.getHeight() / maxDepth));
		}
		else{
			int xmid = x1 + (x2-x1)*v.left.nameLengthSum / (v.left.nameLengthSum + v.right.nameLengthSum);
			if(v.parent != null){
				g.drawLine(xmid, y, xmid, y - (int)(v.edgeLength * this.getHeight() / maxDepth));
			}
			int xleft = (v.left.left != null) ?
						x1 + (xmid - x1) * v.left.left.nameLengthSum / (v.left.left.nameLengthSum + v.left.right.nameLengthSum) :
						x1 + (xmid - x1) / 2;
			int xright = (v.right.left != null) ?
					xmid + (x2 - xmid) * v.right.left.nameLengthSum / (v.right.left.nameLengthSum + v.right.right.nameLengthSum) :
						xmid + (x2 - xmid) / 2;
			g.drawLine(xleft, y, xright, y);
			drawTree(g,v.left,maxDepth,v.left.nameLengthSum,x1,xmid,y+(int)(v.left.edgeLength * this.getHeight() / maxDepth));
			drawTree(g,v.right,maxDepth,v.right.nameLengthSum,xmid,x2,y+(int)(v.right.edgeLength * this.getHeight() / maxDepth));
		}
	}

	/**
	 * Gives the minimum size of the component
	 */
	@Override
	public Dimension getMinimumSize(){
		return getPreferredSize();
	}

	/**
	 * It gives the preferred size of the component
	 */
	@Override
	public Dimension getPreferredSize() {
		return new Dimension(0,100);
	}

}
