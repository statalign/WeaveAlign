package statalign.postprocess.plugins;

import java.awt.BorderLayout;
import java.io.IOException;

import javax.swing.Icon;
import javax.swing.ImageIcon;
import javax.swing.JPanel;
import javax.swing.JScrollPane;

import statalign.postprocess.Postprocess;
import statalign.postprocess.gui.TreeVisualisationGUI;

public class TreeVisualisation extends Postprocess{

	JPanel panel;
	String[] newickTrees;
	int newickTreesIndex;
	private TreeVisualisationGUI gui;
	//private boolean sampling;

	public TreeVisualisation(){
		screenable = true;
		outputable = true;
		postprocessable = true;
	}

	@Override
	public void beforeFirstSample() {
		if(show) {
			panel.removeAll();
			JScrollPane scroll = new JScrollPane();
			gui = new TreeVisualisationGUI(panel, this);
			scroll.setViewportView(gui);
			panel.add(scroll, BorderLayout.CENTER);
		}
		newickTrees = new String[mcmc.mcmcpars.cycles/mcmc.mcmcpars.sampRate];
		newickTreesIndex = 0;
	}

	@Override
	public Icon getIcon() {
		return new ImageIcon(ClassLoader.getSystemResource("icons/tree.gif"));
	}

	@Override
	public JPanel getJPanel() {
		panel = new JPanel(new BorderLayout());
		return panel;
	}

	@Override
	public String getTabName() {
		return "Tree";
	}

	@Override
	public String getTip() {
		return "Current tree in the Markov chain";
	}

	@Override
	public String getFileExtension() {
		return "tree";
	}
	
	@Override
	public void newStep() {
		if(show)
			gui.repaint();
	}

	@Override
	public void newSample(int no, int total) {
		newickTrees[newickTreesIndex] = mcmc.tree.printedTree();
		newickTreesIndex++;
		if(sampling){
			try {
				file.write("Sample "+no+"\tTree:\t"+mcmc.tree.printedTree()+"\n");

			} catch (IOException e) {
				// TODO Auto-generated catch block
				e.printStackTrace();
			}
		}
		else{
			//System.out.println("Not sampling alignment!!!");
		}

	}

	@Override
	public void afterLastSample() {
		try{
			  //FileWriter fw = new FileWriter("output/"+mcmc.tree.title+".tree.nexus", false);
			//implement nexus file

			  outputFile.write("#nexus\n\nBEGIN Taxa;\nDIMENSIONS ntax="+((mcmc.tree.vertex.length+1)/2)+";\nTAXLABELS\n");
			  for(int i = 0; i < mcmc.tree.names.length; i++){
				  outputFile.write("["+(i+1)+"] '"+mcmc.tree.names[i].replaceAll(" ", "")+"'\n");
			  }
			  outputFile.write(";\nEND; [Taxa]\nBEGIN Trees;\n");
			  for(int i = 0; i < newickTrees.length; i++){
				  outputFile.write("["+i+"] tree 'tree_"+i+"'= "+newickTrees[i]+"\n");
			  }
			  outputFile.write("END; [Trees]\n"+
					  "BEGIN st_Assumptions;\n"+
					  "\ttreestransform=TreeSelector;\n"+
					  "\tsplitstransform=EqualAngle;\n"+
					  "\tSplitsPostProcess filter=dimension value=4;\n"+
					  "\tautolayoutnodelabels;\nEND; [st_Assumptions]");
			  outputFile.close();
			}
			catch(IOException e){

			}


	}

	/* (non-Javadoc)
	 * @see statalign.postprocess.Postprocess#setSampling(boolean)
	 */
	@Override
	public void setSampling(boolean enabled) {
		sampling = enabled;

	}

}
