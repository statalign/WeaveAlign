package statalign.base;

import java.io.IOException;
import java.util.Random;

import statalign.base.thread.Stoppable;
import statalign.base.thread.StoppedException;
import statalign.model.score.plugins.Blosum62;
import statalign.model.subst.plugins.Dayhoff;
import statalign.postprocess.PostprocessManager;
import statalign.ui.ErrorMessage;
import statalign.ui.MainFrame;

/**
 * 
 * This class handles an MCMC run.
 * 
 * The class extends <tt>Stoppable</tt>, it may be terminated/suspended
 * in graphical mode.
 * 
 * @author miklos, novak
 *
 */
public class Mcmc extends Stoppable {

	/**
	 * Current tree in the MCMC run
	 */
	public Tree tree;
	
	/**
	 * MCMC parameters including the number of burn-in steps, the total
	 * number of steps in the MCMC and the sampling rate.
	 */
	public MCMCPars mcmcpars;
	
	/**
	 * PostprocessManager that handles the postprocessing modules
	 */
	public PostprocessManager postprocMan;
	
	//	int samplingMethod = 1; //0: random sampling, 1: total sampling
	double[] weights; //for selecting internal tree node
	final static double LEAFCOUNT_POWER = 1.0;
	final static double SELTRLEVPROB[] = {0.9,0.6,0.4,0.2,0};
	final static int FIVECHOOSE[] = {35, 5, 15, 35, 10}; //edge, topology, indel parameter, alignment, substitutionparameter
	final static int FOURCHOOSE[] = {35, 5, 25, 35}; //edge, topology, indel parameter, alignment
	
//	BufferedWriter mpd;
	/**
	 * True while the MCMC is in the burn-in phase.
	 */
	public boolean burnin;
	
	//**
	// * Default constructor with default values
	 //*/
	//Mcmc(){
	//	mcmcpars = new MCMCPars(10000,100000,1000);
	//}

	 
	 public Mcmc(Tree tree, MCMCPars pars, PostprocessManager ppm){
		postprocMan = ppm;
//		ppm.mcmc = this;
//		try {
//			mpd = new BufferedWriter(new FileWriter("alignments.txt"));
//			mpd.write((pars.cycles/pars.sampRate)+" "+(tree.vertex.length/2+1)+"\n\n");
//		} catch (IOException e) {
//			e.printStackTrace();
//		}
		this.tree = tree;
		weights = new double[tree.vertex.length];
		mcmcpars = pars;
	}

	//public Mcmc(Tree tree, int burnin, int cycles, int sampleRate, PostprocessManager ppm){
	//	postprocMan = ppm;
	//	ppm.mcmc = this;
	//	this.tree = tree;
	//	weights = new double[tree.vertex.length];
	//	mcmcpars = new MCMCPars(burnin, cycles, sampleRate);
	//}

	 private int alignmentSampled = 0;
	 private int alignmentAccepted = 0;
	 private int edgeSampled = 0;
	 private int edgeAccepted = 0;
	 private int topologySampled = 0;
	 private int topologyAccepted = 0;
	 private int indelSampled = 0;
	 private int indelAccepted = 0;
	 private int substSampled = 0;
	 private int substAccepted = 0;
	 
	 /**
	  * In effect starts an MCMC run. It first performs a prescribed
	  * number of burn-in steps, then makes the prescribed number of steps after burn-in,
	  * drawing samples with the prescribes frequency. It also calls the appropriate functions
	  * of the PostpocessManager <tt>postprocMan</tt> to trigger data transfer to
	  * postprocessing modules when necessary
	  */
	public void doMCMC(){
		System.out.println("Starting MCMC...\n");

		long start = System.currentTimeMillis();
		long currentTime; 
		Utils.generator = new Random(mcmcpars.seed);
		postprocMan.beforeFirstSample();
		MainFrame frame = postprocMan.mainManager.frame;
		try {
			int burnIn = mcmcpars.burnIn;
			burnin = true;
			for(int i = 0; i < burnIn; i++){
				sample(0);
				postprocMan.newStep();
				currentTime = System.currentTimeMillis();
				if(frame != null) {
					String text = "Burn In: " + (i+1);
					if(i > 10)
						text += "  "+remainingTime((currentTime - start)*(burnIn - i - 1 + mcmcpars.cycles)/(i + 1));
					frame.statusText.setText(text);
				} else if(i % 1000 == 999) {
					System.out.println("Burn in: "+(i+1));
				}
			}
			burnin = false;
			int period = mcmcpars.cycles/mcmcpars.sampRate;
			int sampRate = mcmcpars.sampRate;
			for(int i = 0; i < period; i++){
				for(int j = 0; j < sampRate; j++){
					sample(0);
					postprocMan.newStep();
					currentTime = System.currentTimeMillis();
					if(frame != null) {
						frame.statusText.setText("Sample: "+(i+1)+"  "+remainingTime((currentTime - start)*((period - i - 1) * sampRate + sampRate - j - 1)/(burnIn + i*sampRate + j + 1)));
					}
				}
				if(frame == null) {
					System.out.println("Sample: "+(i+1));
				}
				report(i, period);
			}
		} catch(StoppedException ex) {
			// stopped: report and save state
			// should we still call afterLastSample?
		}
		postprocMan.afterLastSample();
		if(frame != null)
			frame.statusText.setText(MainFrame.IDLE_STATUS_MESSAGE);
	}

	private static String remainingTime(long x) {
		x /= 1000;
		return String.format("Estimated time left: %d:%02d:%02d", x/3600, (x/60)%60, x%60);
	}
	
	
	private void sample(int samplingMethod) throws StoppedException {
		if(samplingMethod == 0){
			stoppable();
			switch(tree.substitutionModel.params != null && tree.substitutionModel.params.length > 0 ? 
					Utils.weightedChoose(FIVECHOOSE) : 
					Utils.weightedChoose(FOURCHOOSE)) {
			case 0:
				sampleEdge();
				break;
			case 1:
				sampleTopology();
				break;
			case 2:
				sampleIndelParameter();
				break;
			case 3:
				sampleAlignment();
				break;
			case 4:
				sampleSubstParameter();
				break;
			}
		} else {
			stoppable();
			sampleEdge();
			sampleTopology();
			sampleIndelParameter();
			sampleSubstParameter();
			sampleAlignment();
		}
	}

	private void sampleEdge(){
		edgeSampled++;
		//System.out.print("Edge: ");
		int i = Utils.generator.nextInt(tree.vertex.length - 1);
		double oldEdge = tree.vertex[i].edgeLength;
		double oldLogLikelihood = tree.getLogLike();
		while((tree.vertex[i].edgeLength = 
			oldEdge + Utils.generator.nextDouble() * Utils.EDGE_SPAN - (Utils.EDGE_SPAN/2.0)) < 0.01);
		tree.vertex[i].edgeChangeUpdate();
		//	Vertex actual = tree.vertex[i];
		//while(actual != null){
		//	actual.calcFelsen();
		//	actual.calcOrphan();
		//	actual.calcIndelLogLike();
		//	actual = actual.parent;
		//}
		tree.vertex[i].calcAllUp();
		double newLogLikelihood = tree.getLogLike();
		if(Utils.generator.nextDouble() <
				(Math.exp(newLogLikelihood - oldLogLikelihood - tree.vertex[i].edgeLength + oldEdge) * 
						(Math.min(oldEdge - 0.01,Utils.EDGE_SPAN/2.0) + Utils.EDGE_SPAN/2.0))/
						(Math.min(tree.vertex[i].edgeLength - 0.01, Utils.EDGE_SPAN/2.0) + Utils.EDGE_SPAN/2.0)){
			//acceptance, do nothing
			//System.out.println("accepted (old: "+oldLogLikelihood+" new: "+newLogLikelihood+")");
			edgeAccepted++;
		} else {
			//reject, restore
			//System.out.print("Rejected! i: "+i+"\tOld likelihood: "+oldLogLikelihood+"\tNew likelihood: "+newLogLikelihood);
			tree.vertex[i].edgeLength = oldEdge;
			tree.vertex[i].edgeChangeUpdate();
			//	actual = tree.vertex[i];
			//while(actual != null){
			//	actual.calcFelsen();
			//	actual.calcOrphan();
			//	actual.calcIndelLogLike();
			//	actual = actual.parent;
			//}
			tree.vertex[i].calcAllUp();
			//System.out.println("rejected (old: "+oldLogLikelihood+" new: "+newLogLikelihood+")");

		}
	}

//this is the old	
	/*
	private void sampleTopology(){
		int vnum = tree.vertex.length;

		if(vnum <= 3)
			return;

		System.out.print("Topology: ");
		double oldLogLi = tree.getLogLike();

		int vertId, rnd = Utils.generator.nextInt(vnum-3);
		vertId = tree.getTopVertexId(rnd);
		if(vertId != -1) {
			int lastId[] = new int[3], num = 0, newId = vertId;

			for(int i = vnum-3; i < vnum; i++) {
				int id = tree.getTopVertexId(i);
				if(id == -1)
					lastId[num++] = i;
				else if(id < vertId)
					newId--;
			}
			rnd = lastId[newId];
		}
		Vertex nephew = tree.vertex[rnd];
		Vertex uncle = nephew.parent.brother();

		//		for(vertId = 0; vertId < vnum; vertId++) {
		//		if(tree.getTopVertexId(vertId) == -1) {				// vertex eligible
		//		if(rnd-- == 0)
		//		break;
		//		}
		//		}
		//		Vertex nephew = tree.vertex[vertId];

		double bpp = nephew.swapWithUncle();

		double newLogLi = tree.getLogLike();

		//	tree.root.calcFelsRecursivelyWithCheck();
		//tree.root.calcIndelRecursivelyWithCheck();

		if(Math.log(Utils.generator.nextDouble()) < bpp+newLogLi-oldLogLi) {
			// accepted
			System.out.println("accepted (old: "+oldLogLi+" new: "+newLogLi+")");
		} else {
			// refused
			uncle.swapBackUncle();
			System.out.println("rejected (old: "+oldLogLi+" new: "+newLogLi+")");
		}

		//tree.root.calcFelsRecursivelyWithCheck();
		//tree.root.calcIndelRecursivelyWithCheck();
	}
*/
	private void sampleTopology(){
		int vnum = tree.vertex.length;

		if(vnum <= 3)
			return;

		topologySampled++;
		//	System.out.println("\n\n\t***\t***\t***\n\n\n");
		//System.out.print("Topology: ");
		//	tree.printAllPointers();
		double oldLogLi = tree.getLogLike();

		int vertId, rnd = Utils.generator.nextInt(vnum-3);
		vertId = tree.getTopVertexId(rnd);
		if(vertId != -1) {
			int lastId[] = new int[3], num = 0, newId = vertId;

			for(int i = vnum-3; i < vnum; i++) {
				int id = tree.getTopVertexId(i);
				if(id == -1)
					lastId[num++] = i;
				else if(id < vertId)
					newId--;
			}
			rnd = lastId[newId];
		}
		Vertex nephew = tree.vertex[rnd];
		Vertex uncle = nephew.parent.brother();

		//		for(vertId = 0; vertId < vnum; vertId++) {
		//		if(tree.getTopVertexId(vertId) == -1) {				// vertex eligible
		//		if(rnd-- == 0)
		//		break;
		//		}
		//		}
		//		Vertex nephew = tree.vertex[vertId];

		//	String[] s = tree.root.printedMultipleAlignment();
		//System.out.println("Alignment before topology changing: ");
		//for(int i = 0; i < s.length; i++){
		//  System.out.println(s[i]);
		//}
		double bpp = nephew.fastSwapWithUncle();
		//double bpp = nephew.swapWithUncle();
		//s = tree.root.printedMultipleAlignment();
		//System.out.println("Alignment after topology changing: ");
		//for(int i = 0; i < s.length; i++){
		//  System.out.println(s[i]);
		//}

		double newLogLi = tree.getLogLike();

		//	tree.root.calcFelsRecursivelyWithCheck();
		//tree.root.calcIndelRecursivelyWithCheck();

		if(Math.log(Utils.generator.nextDouble()) < bpp+newLogLi-oldLogLi) {
			// accepted
			//System.out.println("accepted (old: "+oldLogLi+" new: "+newLogLi+")");
			topologyAccepted++;
		} else {
			// rejected
			//System.out.println("Checking pointer integrity before changing back topology: ");
			for(int i = 0; i < tree.vertex.length; i++){
				if(tree.vertex[i].left != null && tree.vertex[i].right != null){
					tree.vertex[i].checkPointers();
					AlignColumn p;
					//checking pointer integrity
					for(AlignColumn c=tree.vertex[i].left.first; c != null; c = c.next){
						p = tree.vertex[i].first;
						while(c.parent != p && p != null)
							p = p.next;
						if(p == null)
							throw new Error("children does not have a parent!!!"+tree.vertex[i]+" "+tree.vertex[i].print());
					}
					for(AlignColumn c=tree.vertex[i].right.first; c != null; c = c.next){
						p = tree.vertex[i].first;
						while(c.parent != p && p != null)
							p = p.next;
						if(p == null)
							throw new Error("children does not have a parent!!!"+tree.vertex[i]+" "+tree.vertex[i].print());
					}


				}
			}


			uncle.fastSwapBackUncle();
			//System.out.println("Checking pointer integrity after changing back topology: ");
			for(int i = 0; i < tree.vertex.length; i++){
				if(tree.vertex[i].left != null && tree.vertex[i].right != null){
					tree.vertex[i].checkPointers();
					AlignColumn p;
					//checking pointer integrity
					for(AlignColumn c=tree.vertex[i].left.first; c != null; c = c.next){
						p = tree.vertex[i].first;
						while(c.parent != p && p != null)
							p = p.next;
						if(p == null)
							throw new Error("children does not have a parent!!!"+tree.vertex[i]+" "+tree.vertex[i].print());
					}
					for(AlignColumn c=tree.vertex[i].right.first; c != null; c = c.next){
						p = tree.vertex[i].first;
						while(c.parent != p && p != null)
							p = p.next;
						if(p == null)
							throw new Error("children does not have a parent!!!"+tree.vertex[i]+" "+tree.vertex[i].print());
					}
				}
			}
			//uncle.swapBackUncle();
			// s = tree.root.printedMultipleAlignment();
			//System.out.println("Alignment after changing back the topology: ");
			//for(int i = 0; i < s.length; i++){
			//	System.out.println(s[i]);
			//}
			// System.out.println("rejected (old: "+oldLogLi+" new: "+newLogLi+")");
		}

		//tree.printAllPointers();
		//	System.out.println("\n\n\t***\t***\t***\n\n\n");
		tree.root.calcFelsRecursivelyWithCheck();
		tree.root.calcIndelRecursivelyWithCheck();
	}

	private void sampleIndelParameter(){
		indelSampled++;
		switch(Utils.generator.nextInt(3)) {
		case 0:
			//System.out.print("Indel param R: ");
			double oldR = tree.hmm2.params[0];
			double oldLogLikelihood = tree.root.orphanLogLike + tree.root.indelLogLike;
			while((tree.hmm2.params[0] = oldR + Utils.generator.nextDouble() * Utils.R_SPAN - Utils.R_SPAN/2.0) <= 0.0 || 
					tree.hmm2.params[0] >= 1.0);
			for(int i = 0; i < tree.vertex.length; i++){
				tree.vertex[i].updateHmmMatrices();
			}
			tree.root.calcIndelLikeRecursively();
			double newLogLikelihood = tree.root.orphanLogLike + tree.root.indelLogLike;
			if(Utils.generator.nextDouble() < 
					Math.exp(newLogLikelihood - oldLogLikelihood) *
					(Math.min(1.0-oldR,Utils.R_SPAN/2.0) + Math.min(oldR, Utils.R_SPAN/2.0))/
					(Math.min(1.0-tree.hmm2.params[0],Utils.R_SPAN/2.0)+Math.min(tree.hmm2.params[0],Utils.R_SPAN/2.0))){
				//accept, do nothing
				//System.out.println("accepted (old: "+oldLogLikelihood+" new: "+newLogLikelihood+")");
				indelAccepted++;
			}
			else{
				//restore
				tree.hmm2.params[0] = oldR;
				for(int i = 0; i < tree.vertex.length; i++){
					tree.vertex[i].updateHmmMatrices();
				}
				tree.root.calcIndelLikeRecursively();
				//System.out.println("rejected (old: "+oldLogLikelihood+" new: "+newLogLikelihood+")");
			}

			break;
		case 1:
			/////////////////////////////////////////////////
			//System.out.print("Indel param Lambda: ");
			double oldLambda = tree.hmm2.params[1];
			oldLogLikelihood = tree.root.orphanLogLike + tree.root.indelLogLike;
			while((tree.hmm2.params[1] = oldLambda + Utils.generator.nextDouble() * Utils.LAMBDA_SPAN - 
					Utils.LAMBDA_SPAN/2.0) <= 0.0 || tree.hmm2.params[1] >= tree.hmm2.params[2]);
			for(int i = 0; i < tree.vertex.length; i++){
				tree.vertex[i].updateHmmMatrices();
			}
			tree.root.calcIndelLikeRecursively();
			newLogLikelihood = tree.root.orphanLogLike + tree.root.indelLogLike;
			if(Utils.generator.nextDouble() < 
					Math.exp(newLogLikelihood - oldLogLikelihood - tree.hmm2.params[1] + oldLambda) *
					(Math.min(Utils.LAMBDA_SPAN/2.0,tree.hmm2.params[2] - oldLambda) + Math.min(oldLambda, Utils.LAMBDA_SPAN/2.0))/
					(Math.min(Utils.LAMBDA_SPAN/2.0,tree.hmm2.params[2] - tree.hmm2.params[1]) + Math.min(tree.hmm2.params[1],Utils.LAMBDA_SPAN/2.0))){
				//accept, do nothing
				//System.out.println("accepted (old: "+oldLogLikelihood+" new: "+newLogLikelihood+" oldLambda: "+oldLambda+" newLambda: "+tree.hmm2.params[1]+")");
				indelAccepted++;
			}
			else{
				//restore
				tree.hmm2.params[1] = oldLambda;
				for(int i = 0; i < tree.vertex.length; i++){
					tree.vertex[i].updateHmmMatrices();
				}
				tree.root.calcIndelLikeRecursively();
				//System.out.println("rejected (old: "+oldLogLikelihood+" new: "+newLogLikelihood+" oldLambda: "+oldLambda+" newLambda: "+tree.hmm2.params[1]+")");
			}
			break;
		case 2:
			/////////////////////////////////////////////////////////
			//System.out.print("Indel param Mu: ");
			double oldMu = tree.hmm2.params[2];
			oldLogLikelihood = tree.getLogLike();
			while((tree.hmm2.params[2] = oldMu + Utils.generator.nextDouble() * Utils.MU_SPAN - 
					Utils.MU_SPAN/2.0) <= tree.hmm2.params[1]);
			for(int i = 0; i < tree.vertex.length; i++){
				tree.vertex[i].updateHmmMatrices();
			}
			tree.root.calcIndelLikeRecursively();
			newLogLikelihood = tree.getLogLike();
			if(Utils.generator.nextDouble() < 
					Math.exp(newLogLikelihood - oldLogLikelihood - tree.hmm2.params[2] + oldMu) *
					(Utils.MU_SPAN/2.0 + Math.min(oldMu-tree.hmm2.params[1], Utils.MU_SPAN/2.0))/
					(Utils.MU_SPAN/2.0 + Math.min(tree.hmm2.params[2]-tree.hmm2.params[1],Utils.MU_SPAN/2.0))){
				//accept, do nothing
				//System.out.println("accepted (old: "+oldLogLikelihood+" new: "+newLogLikelihood+")");
				indelAccepted++;
			}
			else{
				//restore
				tree.hmm2.params[2] = oldMu;
				for(int i = 0; i < tree.vertex.length; i++){
					tree.vertex[i].updateHmmMatrices();
				}
				tree.root.calcIndelLikeRecursively();
				//System.out.println("rejected (old: "+oldLogLikelihood+" new: "+newLogLikelihood+")");
			}
			break;
		}
	}

	private void sampleSubstParameter(){
		substSampled++;
		if(tree.substitutionModel.params.length == 0) return;
		else{
			double mh = tree.substitutionModel.sampleParameter();
			double oldlikelihood = tree.root.orphanLogLike;
			for(int i = 0; i < tree.vertex.length; i++){
				tree.vertex[i].updateTransitionMatrix();
			}
			tree.root.calcFelsRecursively();
			double newlikelihood = tree.root.orphanLogLike;
			if(Utils.generator.nextDouble() < Math.exp(mh + newlikelihood - oldlikelihood)){
				//System.out.println("Substitution parameter: accepted (old: "+oldlikelihood+" new: "+newlikelihood+")");
				substAccepted++;
			}
			else{
				tree.substitutionModel.restoreParameter();
				for(int i = 0; i < tree.vertex.length; i++){
					tree.vertex[i].updateTransitionMatrix();
				}
				tree.root.calcFelsRecursively();
				//System.out.println("Substitution parameter: rejected (old: "+oldlikelihood+" new: "+newlikelihood+")");
			}
		}
	}

	private void sampleAlignment(){
		alignmentSampled++;
		for(int i = 0; i < tree.vertex.length; i++){
			tree.vertex[i].selected = false;
		}
		//System.out.print("Alignment: ");
		double oldLogLi = tree.getLogLike();
		//		System.out.println("fast indel before: "+tree.root.indelLogLike);
		tree.countLeaves(); // calculates recursively how many leaves we have below this node
		for(int i = 0; i < weights.length; i++){
			weights[i] = Math.pow(tree.vertex[i].leafCount,LEAFCOUNT_POWER);
		}
		int k = Utils.weightedChoose(weights, null);
		//		System.out.println("Sampling from the subtree: "+tree.vertex[k].print());
		tree.vertex[k].selectSubtree(SELTRLEVPROB, 0);
		double bpp = tree.vertex[k].selectAndResampleAlignment();
		double newLogLi = tree.getLogLike();
		//		double fastFels = tree.root.orphanLogLike;
		//		double fastIns = tree.root.indelLogLike;
		//		report();
		//		tree.root.first.seq[0] = 0.0;
		//		System.out.println("Old before: "+tree.root.old.indelLogLike);
		//	tree.root.calcFelsRecursivelyWithCheck();
		//tree.root.calcIndelRecursivelyWithCheck();
		//		tree.root.calcIndelLikeRecursively();
		//		System.out.println("Old after: "+tree.root.old.indelLogLike);
		//		System.out.println("Check logli: "+tree.getLogLike()+" fastFels: "+fastFels+" slowFels: "+tree.root.orphanLogLike+
		//		" fastIns: "+fastIns+" slowIns: "+tree.root.indelLogLike);
		//		System.out.println("selected subtree: "+tree.vertex[k].print());
		if(Math.log(Utils.generator.nextDouble()) < bpp+newLogLi-oldLogLi) {
			// accepted
			//System.out.println("accepted (old: "+oldLogLi+" new: "+newLogLi+")");
			alignmentAccepted++;
		} else {
			// refused
			// String[] s = tree.printedAlignment();
			tree.vertex[k].alignRestore();
			//s = tree.printedAlignment();
			//System.out.println("rejected (old: "+oldLogLi+" new: "+newLogLi+")");
			//			System.out.println("after reject fast: "+tree.root.indelLogLike);
			//			tree.root.calcIndelRecursivelyWithCheck();
			//			System.out.println(" slow: "+tree.root.indelLogLike);
		}
		//tree.root.calcFelsRecursivelyWithCheck();
		//tree.root.calcIndelRecursivelyWithCheck();
	}

	void report(int no, int total){
		postprocMan.newSample(no, total);

		try {
			postprocMan.logFile.write("Acceptances\tAlignment\t"+(alignmentSampled == 0 ? 0 : (double)alignmentAccepted/(double)alignmentSampled)+"\t"+
									       "Edge\t"+(     edgeSampled == 0 ? 0 : (double)     edgeAccepted/(double)     edgeSampled)+"\t"+
									   "Topology\t"+( topologySampled == 0 ? 0 : (double) topologyAccepted/(double) topologySampled)+"\t"+
									      "Indel\t"+(    indelSampled == 0 ? 0 : (double)    indelAccepted/(double)    indelSampled)+"\t"+
								   "Substitution\t"+(    substSampled == 0 ? 0 : (double)    substAccepted/(double)    substSampled)+"\n");
		} catch (IOException e) {
			new ErrorMessage(null,e.getLocalizedMessage(),true);
		}
		alignmentSampled = 0;
		alignmentAccepted = 0;
		edgeSampled = 0;
		edgeAccepted = 0;
		topologySampled = 0;
		topologyAccepted = 0;
		indelSampled = 0;
		indelAccepted = 0;
		substSampled = 0;
		substAccepted = 0;

		// move this to a postprocessing plugin
		try {
			postprocMan.logFile.write("Report\tLogLikelihood\t"+(tree.root.orphanLogLike+tree.root.indelLogLike)+
				       "\tR\t"+tree.hmm2.params[0]+
				       "\tLamda\t"+tree.hmm2.params[1]+
				       "\tMu\t"+tree.hmm2.params[2]+"\t"+tree.substitutionModel.print()+"\n");
		} catch (IOException e) {
			e.printStackTrace();
		}
		//System.out.println("Tree\t"+tree.printedTree());
		//String[] s = tree.printedAlignment("StatAlign");
		//System.out.println("Alignment\t");
		//for(int i = 0; i < s.length; i++){
		//	System.out.println("Alignment\t"+s[i]);
		//}

		//double check
		//	if(s[0].charAt(s[0].length() - 1) == 'K'){
		//  for(int i = 0; i < tree.vertex.length; i++){
		//	String[] s1 = tree.vertex[i].printedAlignment();
		//	System.out.println("\n"+tree.vertex[i].print(0)+"\n"+s1[0]+"\n"+s1[1]+"\n");
		//  }
		//}
	}

	/**
	 * This function is only for testing and debugging purposes.
	 * 
	 * @param args Actually, we do not use these parameters, as this function is for testing
	 * and debugging purposes. All necessary input data is directly written into the function.
	 * 
	 */
	public static void main(String[] args) {
		try {
			Tree tree = new Tree(new String[] {
							"kkkkkkwwwwwwwwlidwwwwwkkk", 
							"kkkwwwwwwwlidwwwwwkkk", 
							"kkkwwwwwwwlidwwwwwkkk", 
							"kkkwwwwwwwlidwwwwwkkk",
							"kkkwwwwwlidwwwwwkkkddkldkl",
							"kkkwwwwwlidwwwwwkkkeqiqii",
							"kkkwwwwwlidwwwwwkkkddkidkil",
							"kkkwwwwwlidwwwwwkkkeqiq",
							"kkkwwwwwlidwwwwwkkkddkldkll",
							"kkkwwwwwlidwwwwwkkkddkldkil"},
			    new String[] {"A", "B", "C", "D", "E", "F", "G", "H", "I", "J"},
			new Dayhoff(),
			new Blosum62(),"");
			for(int i = 0; i < tree.vertex.length; i++){
				//    tree.vertex[i].edgeLength = 0.1;
				tree.vertex[i].updateTransitionMatrix();
			}
			tree.root.calcFelsRecursively();
			System.out.println(tree.printedTree());
			Mcmc mcmc = new Mcmc(tree, new MCMCPars(0, 10000, 10, 1L), new PostprocessManager(null));
			mcmc.doMCMC();
		} catch(StoppedException e) {
			// stopped during tree construction
		} catch(IOException e) {
		}
	}

}