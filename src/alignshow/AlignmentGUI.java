package alignshow;

import java.awt.BasicStroke;
import java.awt.Color;
import java.awt.Dimension;
import java.awt.Font;
import java.awt.Graphics;
import java.awt.Graphics2D;
import java.awt.Stroke;
import java.awt.image.BufferedImage;
import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;

import javax.imageio.ImageIO;
import javax.swing.JPanel;


import org.sourceforge.jlibeps.epsgraphics.EpsGraphics2D;

import wvalign.model.SubstitutionModel;


/**
 * This is the graphical interface for showing alignments.
 * 
 * @author miklos,novak
 *
 */
public class AlignmentGUI extends JPanel{

/**
	 * 
	 */
	private static final long serialVersionUID = 1L;
	//	JPanel pan;
//	CurrentAlignment owner;
	SubstitutionModel subst;
	/**
	 *  The alignment that is written onto the screen is stored in this array
	 */
	public String[] alignment = null;

    /**
       The range (specified in terms of the index along refSeq) to print
    */
    public int[] range = new int[2];
    public String refSeq;
    
	/**
	 * The title of the current analysis
	 */
	public String title;
	/**
	 * These are the posterior decoding values
	 */
	public List<List<Double>> tracks = new ArrayList<List<Double>>();
	public List<Color> colors = new ArrayList<Color>();
	
	public List<List<Integer>> grouping = new ArrayList<List<Integer>>();
	
	private boolean showTitle = true;
	
	private int aligNum = 15;
	private int aligLen = 200;

	static final int COLUMN_WIDTH = 11;
	static final int FONT_HEIGHT = 15;
	static final int OFFSET_X = 10;
	static final int OFFSET_Y = 10;
	static final int TITLE_Y = 12;

	//static final String magentaCharacters = "RHK";
	//static final String     redCharacters = "AVFPMILW";
	//static final String    blueCharacters = "DE";
	//static final String   greenCharacters = "STYCNGQ";

/**
 * It initializes a new panel for alignments.
 * 
 * @param title This is the title of the current analysis
 * @param subst The substitution model defines the background colors of characters
 */
	public AlignmentGUI(String title, SubstitutionModel subst) {
		this.title = title;
		this.subst = subst;
	}

    public void setAlignment(String[] ali) {
	alignment = ali;
	range[0] = 1;
	range[1] = ali[0].length();
    }
    public void setRange(int a, int b) {
	range[0] = Math.max(a,1);
	range[1] = Math.min(b,alignment[0].length());
    }
	private static boolean allTab(String[] s, int p){
		boolean b = true;
		for(int i = 0; i < s.length && b; i++){
			b = s[i].charAt(p) == '\t';
		}
		return b;
	}

	/**
	 * This function updates the graphics
	 */
	@Override
	public void paintComponent(Graphics gr){
		paintFunc(gr, true);
	}
	
	private void paintFunc(Graphics gr, boolean paintBackground) {
		double max = 1.0;
		//setSize((int)(owner.pan.getSize().width*0.9),(int)(owner.pan.getSize().height*0.9));
		// text.setRows(pan.getHeight()/13);
		// text.setColumns((int)(pan.getWidth()/6.6));
		Graphics2D g = (Graphics2D)gr;

		if(paintBackground) {
			g.setBackground(Color.white);
			g.clearRect(0, 0, getWidth(), getHeight());
		}
//		alignment = owner.allAlignment;
//		title = owner.title;
		if(alignment != null && alignment[0] != null) {
			int colHeight = (tracks.size() == 0 ? -OFFSET_Y :
				Math.min(200, getHeight() - (4 * OFFSET_Y + 
				(showTitle?TITLE_Y:-OFFSET_Y) + alignment.length * FONT_HEIGHT)));
			
			//Find the first all space column
			int tab =alignment[0].length()-2;
			while(!allTab(alignment, tab)){
				tab--;
			}
			if (range[1]==alignment[0].length()) range[1] -= (tab+1);

			if(showTitle) {
				g.setColor(Color.BLACK);
				g.setFont(new Font("SANS_SERIF", Font.BOLD, 16));
				g.drawString(title, OFFSET_X, OFFSET_Y+TITLE_Y);
			}
			g.setFont(new Font("MONOSPACED", Font.PLAIN, 12));
			
			int startY = colHeight + 3*OFFSET_Y + (showTitle?TITLE_Y:-OFFSET_Y);
			
			// print alignment
			for (int i = 0; i < alignment.length; i++) {
				//System.out.println(alignment[i]);
			    for (int j = 0; j < Math.min(alignment[i].length(),range[1]+tab+1); j++) {
				int jPos = j;
				if(j>tab){
				    jPos -= range[0]-1;
					    if (j<range[0]+tab) continue;
						Color color = subst.getColor(alignment[i].charAt(j));
						if(grouping.size() > 0) {
							double mod = grouping.get(j-tab-1).contains(i) ? 1.3 : 0.65;
							color = new Color(lims((int)(color.getRed()*mod+.5)),
									lims((int)(color.getGreen()*mod+.5)),
									lims((int)(color.getBlue()*mod+.5)));
						}
						g.setColor(color);
						g.fillRect(OFFSET_X + COLUMN_WIDTH * (j-range[0]+1),
								startY + FONT_HEIGHT * i + 3, 
								COLUMN_WIDTH+1,
								FONT_HEIGHT+1);
						//System.out.println((OFFSET_X + COLUMN_WIDTH * j)+" "+(colHeight + OFFSET_Y + TITLE_Y + FONT_HEIGHT * i + 3));
					}
					g.setColor(Color.BLACK);
					g.drawString(alignment[i].charAt(j) + "",
						     OFFSET_X + COLUMN_WIDTH * jPos + 2,
							startY + FONT_HEIGHT * (i+1));

				}
			}
			
			// colour groupings
			g.setStroke(new BasicStroke(2, BasicStroke.CAP_ROUND, BasicStroke.JOIN_ROUND));
			if(grouping != null && grouping.size() > 0) {
				int len = alignment[0].length()-tab-1;
				int hg = 1; // horizontal group
//				System.out.println(len+","+grouping.size());
				for(int i = 0; i < len; i += hg) {
					hg = 1;
					while(i+hg < len && grouping.get(i).equals(grouping.get(i+hg)))
						hg++;
					int gs = grouping.get(i).size();
					int vg = 1;	// vertical group
					for(int j = 0; j < gs; j += vg) {
						vg = 1;
						while(j+vg < gs && grouping.get(i).get(j)+vg == grouping.get(i).get(j+vg))
							vg++;
						g.setColor(Color.black);
						g.drawRect(OFFSET_X + COLUMN_WIDTH * (tab+1+i),
								startY + FONT_HEIGHT * grouping.get(i).get(j) + 2,
								COLUMN_WIDTH * hg, FONT_HEIGHT * vg + 2);
//						g.setPaint(new Color(80, 80, 160, 60));
//						g.fillRect(OFFSET_X + COLUMN_WIDTH * (tab+1+i),
//								colHeight + OFFSET_Y + TITLE_Y + FONT_HEIGHT * grouping.get(i).get(j) + 2,
//								COLUMN_WIDTH * hg, FONT_HEIGHT * vg + 2);
					}
				}
			}

			// draw graph axes
			Stroke st = g.getStroke();
			Stroke basic = new BasicStroke(1);
			Stroke dotted = new BasicStroke(1, BasicStroke.CAP_BUTT, BasicStroke.JOIN_MITER, 1, new float[] {2, 2}, 0);
			if(tracks.size() != 0) {
				g.setStroke(dotted);
				g.drawLine(OFFSET_X + (tab+1) * COLUMN_WIDTH - 10, startY - OFFSET_Y - colHeight,
					   OFFSET_X + (tab+1+Math.min(range[1]-range[0]+1,tracks.get(0).size())) * COLUMN_WIDTH, startY - OFFSET_Y - colHeight);
				g.drawLine(OFFSET_X + (tab+1) * COLUMN_WIDTH - 10, startY - OFFSET_Y,
					   OFFSET_X + (tab+1+Math.min(range[1]-range[0]+1,tracks.get(0).size())) * COLUMN_WIDTH, startY - OFFSET_Y);
				g.drawLine(OFFSET_X + (tab+1) * COLUMN_WIDTH, startY - OFFSET_Y - colHeight,
						OFFSET_X + (tab+1) * COLUMN_WIDTH, startY - OFFSET_Y);
				g.setStroke(basic);
				g.drawString("0.0", OFFSET_X + (tab+1) * COLUMN_WIDTH - 35, startY - OFFSET_Y + 3);
				g.drawString(""+max, OFFSET_X + (tab+1) * COLUMN_WIDTH - 35, startY - OFFSET_Y - colHeight + 5);
				g.setStroke(st);
			}
			
			// draw tracks
			dotted = new BasicStroke(2, BasicStroke.CAP_BUTT, BasicStroke.JOIN_MITER, 2, new float[] {4, 3}, 0);
			for(int track = 0; track < tracks.size(); track++) {
//				System.out.println(decoding[1]);
				g.setColor(colors.get(track));
				//int colHeight = Math.min(100, g.getClipBounds().height - (2 * OFFSET_Y + TITLE_Y + alignment.length * FONT_HEIGHT));
				List<Double> tr = tracks.get(track);
				for(int i = range[0]; i < Math.min(tr.size(),range[1]); i++){
					if (tr.get(i-1) == null || Double.isNaN(tr.get(i-1))) { continue; }
					if (tr.get(i) == null || Double.isNaN(tr.get(i))) { ++i; continue; }
					g.drawLine( OFFSET_X + (tab + i - range[0] + 1) * COLUMN_WIDTH + COLUMN_WIDTH / 2,
							(int)(startY - OFFSET_Y - Math.round(tr.get(i-1) * (colHeight) / max)),
//							(int)(OFFSET_Y + TITLE_Y - Math.round(decoding[i-1] * (TITLE_Y) / max)),
							OFFSET_X + (tab + i - range[0] + 2) * COLUMN_WIDTH + COLUMN_WIDTH / 2,
							(int)(startY - OFFSET_Y - Math.round(tr.get(i) * (colHeight) / max))
//							(int)(OFFSET_Y + TITLE_Y - Math.round(decoding[i] * (TITLE_Y) / max))
							);
				}
			}

		}
		else{
			g.setColor(Color.BLACK);
			g.setFont(new Font("SANS_SERIF", Font.BOLD, 16));
			g.drawString("Waiting for data...", OFFSET_X, TITLE_Y);
		}
	}

	private int lims(int i) {
		return Math.max(Math.min(i, 255), 0);
	}

	/**
	 * This function tells the minimum size of the panel
	 */
	@Override
	public Dimension getMinimumSize(){
		return getPreferredSize();
	}

	/**
	 * This function tells the preferred size of the panel
	 */
	@Override
	public Dimension getPreferredSize() {
	    int len = alignment == null ? aligLen : Math.min(range[1]-range[0]+5,alignment[0].length());
		int num = alignment == null ? aligNum : alignment.length;
		return new Dimension(OFFSET_X + COLUMN_WIDTH * len + 30, 4*OFFSET_Y + (showTitle?TITLE_Y:-OFFSET_Y) + FONT_HEIGHT * num + (tracks.size() == 0 ? -OFFSET_Y : 150));
	}
	
	public void setShowTitle(boolean selected) {
		showTitle = selected;
		revalidate();
		repaint();
	}

	public void saveAsImage(File file) throws IOException {
		BufferedImage im = new BufferedImage(getWidth(), getHeight(), BufferedImage.TYPE_INT_RGB);
		Graphics g = im.getGraphics();
		paintFunc(g, true);
		ImageIO.write(im, "png", file);
	}
	
	public void saveAsEps(File file) throws IOException {
		EpsGraphics2D g = new EpsGraphics2D("eps", file, 0, 0, getWidth(), getHeight());
		paintFunc(g, false);
		g.flush();
		g.close();
	}
}
