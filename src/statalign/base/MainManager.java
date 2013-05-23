package statalign.base;

import java.io.FileWriter;
import java.io.IOException;

import javax.swing.SwingUtilities;

import statalign.base.thread.MainThread;
import statalign.io.RawSequences;
import statalign.model.subst.SubstitutionModel;
import statalign.postprocess.Postprocess;
import statalign.postprocess.PostprocessManager;
import statalign.postprocess.gui.InputGUI;
import statalign.ui.MainFrame;

/**
 * 
 * This is the central class in the information flow amongst classes.
 * 
 * This class also manages the main thread that runs the MCMC.
 * 
 * @author miklos, novak
 *
 */
public class MainManager {
	
	/**
	 * The loaded sequences are stored in the RawSequences class 
	 */
	public RawSequences seqs = new RawSequences();
	
	/**
	 * This class stores the parameters (number of burn-in cycles, number of steps, number of
	 * samplings) of the MCMC
	 */
	public MCMCPars pars = new MCMCPars(10000,100000,1000,1);
	
	/**
	 * This is the postprocess-manager of our program, which handles the postprocesses applied
	 * on the MCMC run.
	 */
	public PostprocessManager postProcMan;

	/**
	 * Array of substitution model classes that can be selected for an analysis.
	 */
	public Class<? extends SubstitutionModel>[] substModels;
	
	/**
	 * The main window of the program.
	 */
	public MainFrame frame;
	
	/**
	 * The current substitution model that is used to analyse the sequences.
	 */
	public SubstitutionModel model;
	
	/**
	 * The graphical interface of the Input panel
	 */
	public InputGUI inputgui;

	/**
	 * Main (background) calculation thread of the application
	 */
	public MainThread thread;

	/**
	 * The full path of the input file from which we read the sequences.
	 */
	public String fullPath;

	/**
	 * Alignment formats in which StatAlign can generate output
	 * 
	 * Implemented formats currently are <tt>StatAlign</tt> (our own format),
	 * <tt>Clustal</tt>, <tt>Fasta</tt>, <tt>Phylip</tt>, <tt>Nexus</tt>
	 */
	public static String[] alignmentTypes = new String[] {"StatAlign", "Clustal", "Fasta", "Phylip", "Nexus"};
	
	/**
	 * This integer stores the index of the current alignment type, as it is in <tt>alignmentTypes</tt>
	 */
	public int currentAlignmentType = 0;
	
	/**
	 * A trivial constructor that only sets <tt>MainFrame</tt> and substModels.
	 * 
	 * The MainFrame creates its MainManager, and the MainManager knows who is
	 * its owner MainFrame, so it can access GUIs on the MainFrame
	 * 
	 * @param frame The owner of the MainManager, the main window of the graphical interface.
	 */
	public MainManager(MainFrame frame) {
		this.frame = frame;
		postProcMan = new PostprocessManager(this);
		
	}
	
	/**
	 * This function starts a new MCMC run.
	 * 
	 * It asks files for writing outputs, and it launches a <tt>MainThread</tt>
	 * that handles the MCMC run.
	 */
	public void start() {
		
		try {
			postProcMan.logFile = new FileWriter(fullPath+".log");

			for(Postprocess p : postProcMan.plugins){
				if(p.postprocessWrite){
					String name = fullPath+"."+p.getFileExtension();
					System.out.println("Output file for "+p.getTabName()+": "+name);
					p.outputFile = new FileWriter(name);
				}
			}

			thread = new MainThread(this);
			thread.start();

		} catch (IOException e) {
			e.printStackTrace();
		}

	}

	/**
	 * Called when the MCMC thread terminates, signals end of the process back to MainFrame.
	 */
	public void finished() {
		if(frame != null) {
			SwingUtilities.invokeLater(new Runnable() {
				public void run() {
					frame.finished();
				}
			});
		}		
	}
}
