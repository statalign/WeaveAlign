package statalign.base.thread;

import java.io.File;

import statalign.base.MainManager;
import statalign.base.Mcmc;
import statalign.base.Tree;
import statalign.io.RawSequences;
import statalign.ui.ErrorMessage;

/**
 * The main (suspendable) thread for background MCMC calculation.
 *
 * @author novak
 */
public class MainThread extends StoppableThread {
	
	/**
	 * Reference to the (singleton) MainManager object encapsulating all settings and data
	 * that an MCMC run depends on.
	 */
	public MainManager owner;
	
	/**
	 * Constructs a new MainThread that can be used to fire a background MCMC calculation.
	 * 
	 * @param owner Reference to the MainManager object.
	 */
	public MainThread(MainManager owner) {
		this.owner = owner;
	}
	/**
	 * Start background MCMC calculation.
	 */
	@Override
	public synchronized void run() {
		try {
			RawSequences seqs = owner.seqs;
			
			if(owner.frame != null) {
				owner.frame.statusText.setText(" Generating initial tree and alignment...");
			}

			System.out.println("\nPreparing initial tree and alignment...\n");

			Tree tree = new Tree(seqs.sequences.toArray(new String[seqs.sequences.size()]), seqs.seqNames.toArray(new String[seqs.seqNames.size()]), 	
					owner.model,
					owner.model.attachedScoringScheme,
					new File(owner.fullPath).getName());
			Mcmc mcmc = new Mcmc(tree, owner.pars, owner.postProcMan);
			mcmc.doMCMC();
		} catch(StoppedException e) {
			// stopped during tree construction
		} catch(Exception e) {
			if(owner.frame != null)
				ErrorMessage.showPane(owner.frame,e,true);
			else
				e.printStackTrace();
		}
		System.out.println("Ready.");
		owner.finished();
	}
}
