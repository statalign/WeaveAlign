package mpd;

import java.awt.BorderLayout;
import java.io.IOException;

import javax.swing.Icon;
import javax.swing.ImageIcon;
import javax.swing.JPanel;
import javax.swing.JScrollPane;

import statalign.postprocess.gui.AlignmentGUI;

/**
 * This is the postprocessmanager for showing the current alignment.
 * 
 * @author miklos, novak
 *
 */
public class CurrentAlignment extends statalign.postprocess.Postprocess{

	public JPanel pan;
	public String[] allAlignment;
	public String[] leafAlignment;
	public String title;
	//private boolean sampling = true;
	
	AlignmentGUI gui;
	
	/**
	 * It construct a postprocess that is screenable (=can be shown on the GUI),
	 * outputable (= can be written into th elog file) but not postprocessable
	 * (= cannot generate its own output file).
	 */
	public CurrentAlignment(){
		screenable = true;
		outputable = true;
		postprocessable = false;
	}
	
	/**
	 * It constructs a new JPanel, and returns with it
	 */
	@Override
	public JPanel getJPanel() {
		pan = new JPanel(new BorderLayout());
		return pan;
	}

	/**
	 * It generates a new icon based on the figure in file icons/calignment.gif
	 */
	@Override
	public Icon getIcon() {		
		return new ImageIcon(ClassLoader.getSystemResource("icons/calignment.gif"));
	}

	/**
	 * It returns with its tab name, 'Alignment'
	 */
	@Override
	public String getTabName() {
		return "Alignment";
	}

	/**
	 * It returns with the tip of the tabulated panel 'Current alignment in the Markov chain'
	 * (shown when mouse is moved over the panel)
	 */
	@Override
	public String getTip() {
		return "Current alignment in the Markov chain";
	}

	@Override
	public String getFileExtension() {
		return "aln";
	}

	private int alignmentIndex = 100;
	/**
	 * After each MCMC step, the current alignment is shown on the GUI.
	 */
	@Override
	public void newStep(){
		alignmentIndex++;
		if(alignmentIndex >= 100){
			alignmentIndex = 0;
			allAlignment = mcmc.tree.printedAlignment("StatAlign");
			int ind = 0, i;
			for(i = 0; i < allAlignment.length; i++)
				if(allAlignment[i].charAt(0) != ' ')
					leafAlignment[ind++] = allAlignment[i];
			if(show) {
				gui.alignment = allAlignment;
				gui.repaint();
			}
	//	alignment = "";
		//for(int i = 0; i < a.length; i++){
			//alignment += a[i]+"\n";
	//	}
		}
	}

	/**
	 * It initialize the graphical interface.
	 */
	@Override
	public void beforeFirstSample() {
		if(show) {
			pan.removeAll();
			title = mcmc.tree.title;
			JScrollPane scroll = new JScrollPane();
			scroll.setViewportView(gui = new AlignmentGUI(title,mcmc.tree.substitutionModel));//, mcmc.tree.printedAlignment()));
			pan.add(scroll, BorderLayout.CENTER);
		}
		leafAlignment = new String[(mcmc.tree.vertex.length+1)/2];		
	}

	/**
	 * At a new MCMC sampling point, it writes the current alignment into the log file, if the
	 * sampling mode is on.
	 */
	public void newSample(int no, int total) {
		if(sampling){
			try {
				String[] alignment = mcmc.tree.printedAlignment(alignmentType);
				for(int i = 0; i < alignment.length; i++){
					file.write("Sample "+no+"\tAlignment:\t"+alignment[i]+"\n");
				}
			} catch (IOException e) {
				// TODO Auto-generated catch block
				e.printStackTrace();
			}
		}
		else{
			//System.out.println("Not sampling alignment!!!");
		}
	}
	
	/**
	 * It switches on or off the sampling mode
	 */
	@Override
	public void setSampling(boolean enabled){
		sampling = enabled;
	}
	
	/**
	 * Empty function, since there is nothing to do after the last sample in this postprocess thread.
	 */
	@Override
	public void afterLastSample() {
	
		
	}

}


/*
package statalign.postprocess.plugins;

import javax.swing.Icon;
import javax.swing.ImageIcon;
import javax.swing.JPanel;

import statalign.postprocess.*;
import statalign.postprocess.gui.CurrentAlignmentGUI;

public class CurrentAlignment extends Postprocess{

	JPanel panel;
	public String[] alignment;
	CurrentAlignmentGUI cagui;
	
	@Override
	public void beforeFirstSample() {
		cagui = new CurrentAlignmentGUI(panel, this);
		String[] alignment = mcmc.tree.printedAlignment();
		String s = "";
		for(int i = 0; i < alignment.length; i++){
			s += alignment[i]+"\n";
		}
		cagui.text.append(s);
		panel.add(cagui);
		
	}

	@Override
	public Icon getIcon() {
		return new ImageIcon("icons/calignment.gif");	}

	@Override
	public JPanel getJPanel() {
		panel = new JPanel();
		return panel;
	}

	@Override
	public String getTabName() {
		return "Alignment";
	}

	@Override
	public String getTip() {
		return "Current alignment";
	}

	@Override
	public void newSample() {
		
		
	}

	@Override
	public void newStep() {
		alignment = mcmc.tree.printedAlignment();
		
	}

}

*/