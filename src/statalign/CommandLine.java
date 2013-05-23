package statalign;

import java.io.IOException;
import java.lang.reflect.Modifier;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.LinkedList;
import java.util.List;
import java.util.ListIterator;
import java.util.Map;

import ml.options.OptionSet;
import ml.options.Options;
import ml.options.Options.Multiplicity;
import ml.options.Options.Separator;
import statalign.base.MainManager;
import statalign.base.Utils;
import statalign.io.RawSequences;
import statalign.io.input.plugins.FastaReader;
import statalign.model.subst.RecognitionError;
import statalign.model.subst.SubstitutionModel;
import statalign.model.subst.plugins.Dayhoff;
import statalign.model.subst.plugins.Kimura3;
import statalign.postprocess.Postprocess;
import statalign.postprocess.PostprocessManager;

public class CommandLine {

	private static String getUsageString(MainManager man) {
		return "Usage:\n\n" +
				"    java -Xmx512m -jar statalign.jar [options] seqfile1 [seqfile2 ...]\n\n\n"
				+
				"Description:\n\n" +
				"    StatAlign can be used for Bayesian analysis of protein, DNA and RNA sequences.\n" +
				"    Multiple alignments, phylogenetic trees and evolutionary parameters are co-estimated\n" +
				"    in a Markov Chain Monte Carlo framework. The input sequence files must be in\n" +
				"    Fasta format.\n\n\n"
				+
				"Options:\n\n" +
				"    -subst=MODEL\n" +
				"        Lets you select from the present substitution models (see list below)\n" +
				"        Default: " + Kimura3.class.getSimpleName() + " (for DNA/RNA data), "+Dayhoff.class.getSimpleName()+ " (for protein data)\n\n"
				+
				"    -mcmc=burn,cycl,rate\n" +
				"        Sets MCMC parameters: burn-in, cycles after burn-in, sampling rate.\n" +
				"          Abbreviations k and m mean 1e3 and 1e6 factors.\n" +
				"        Default: 10k,100k,1k\n\n"
				+
				"    -seed=value\n" +
				"        Sets the random seed (same value will reproduce same results for\n" +
				"          identical input and settings)\n" +
				"        Default: 1\n\n"
				+
				"    -ot=OUTTYPE\n" +
				"        Sets output alignment type (one of "+Utils.joinStrings(MainManager.alignmentTypes, ", ")+")\n" +
				"        Default: "+MainManager.alignmentTypes[0]+"\n\n"
				+
				"    -log=[" + Utils.joinStrings(postprocAbbr.keySet().toArray(), "][,") + "]\n" +
				"        Lets you customise what is written into the log file (one entry for each sample).\n"+
							buildPpListStr(man, "          ")+
				"        Default: "+buildDefPpList(man)+"\n";
	}

	private static List<String> substModNames = new ArrayList<String>();
	private static Map<String, Integer> postprocAbbr = new HashMap<String, Integer>();

	/**
	 * Fills run-time parameters using a list of command-line arguments. If
	 * error occurs displays error message or usage information.
	 * 
	 * @param args
	 *            list of command-line arguments
	 * @param params
	 *            Parameters object to fill out
	 * @return 0 on success, 1 if usage info and 2 if error msg has been
	 *         displayed
	 */
	public static int fillParams(String[] args, MainManager manager) {
		initArrays(manager);
		
		Options opt = new Options(args, Multiplicity.ZERO_OR_ONE, 1,
				Integer.MAX_VALUE);
		opt.addSet("run")
				.addOption("subst", Separator.EQUALS)
				.addOption("mcmc", Separator.EQUALS)
				.addOption("seed", Separator.EQUALS)
				.addOption("ot", Separator.EQUALS)
				.addOption("log", Separator.EQUALS);

		OptionSet set;
		if ((set = opt.getMatchingSet(false, false)) == null) {
			return usage(manager);
		}

		FastaReader f = new FastaReader();
		try {
			for (String inputFile : set.getData()) {
				try {
					RawSequences seq = f.read(inputFile);
					int errorNum = f.getErrors();
					if (errorNum > 0)
						warning(errorNum + " errors occured reading " + inputFile
								+ ", please check file format");
					// TODO this should be printed only when requested
					System.out.println("INFO: read " + seq.size()
							+ " sequences from " + inputFile);
					manager.seqs.add(seq);
				} catch (IOException e) {
					return error("error reading input file " + inputFile);
				}
			}

			manager.fullPath = set.getData().get(0);

			if(set.isSet("subst")) {
				String modelName = set.getOption("subst").getResultValue(0);
				for(String model : substModNames) {
					Class<?> cl = Class.forName(model);
					if(cl.getSimpleName().equalsIgnoreCase(modelName)) {
						manager.model = (SubstitutionModel) cl.newInstance();
						try {
							manager.model.acceptable(manager.seqs);
						} catch(RecognitionError e){
							return error("Substitution model "+modelName+" does not accept the given input sequences!");
						}
						break;
					}
				}
				if(manager.model == null) {
					return error("Unknown substitution model: "+modelName+"\n");
				}
			} else {
				manager.model = new Kimura3();
				try {
					manager.model.acceptable(manager.seqs);
					System.out.println("Automatically selected "+manager.model.getClass().getSimpleName()+" as substitution model.");
				} catch(RecognitionError e){
					try {
						manager.model = new Dayhoff();
						manager.model.acceptable(manager.seqs);
						System.out.println("Automatically selected "+manager.model.getClass().getSimpleName()+" as substitution model.");
					} catch(RecognitionError ee){
						return error("Default substitution model "+manager.model.getClass().getSimpleName()+" does not accept the given input sequences!");
					}
				}
			}

			if(set.isSet("mcmc")) {
				String mcmcPars = set.getOption("mcmc").getResultValue(0);
				String[] pars = mcmcPars.split(",");
				if(pars.length != 3) {
					return error("MCMC parameters not recognized: "+mcmcPars);
				}
				manager.pars.burnIn = parseValue(pars[0]);
				manager.pars.cycles = parseValue(pars[1]);
				manager.pars.sampRate = parseValue(pars[2]);
				if(manager.pars.burnIn < 0 || manager.pars.cycles < 0 || manager.pars.sampRate < 0) {
					return error("MCMC parameters not recognized: "+mcmcPars);
				}
			}

			if(set.isSet("seed")) {
				String seedPar = set.getOption("seed").getResultValue(0);
				try {
					manager.pars.seed = Integer.parseInt(seedPar);
				} catch (NumberFormatException e) {
					return error("error parsing seed parameter: "+seedPar);
				}
			}
			
			if(set.isSet("ot")) {
				String outType = set.getOption("ot").getResultValue(0);
				int i;
				for(i = 0; i < MainManager.alignmentTypes.length; i++){
					if(outType.equalsIgnoreCase(MainManager.alignmentTypes[i])){
						manager.currentAlignmentType = i;
						break;
					}
				}
				if(i == MainManager.alignmentTypes.length) {
					return error("Unknown output type: "+outType+"\n");
				}
			}

			if(set.isSet("log")) {
				String log = set.getOption("log").getResultValue(0);
				String[] keys = log.split(",");
				Postprocess[] pps = manager.postProcMan.plugins;
				
				for(Postprocess pp : pps)
					pp.sampling = false;
				
				for(String key : keys) {
					try {
						pps[postprocAbbr.get(key.toUpperCase())].sampling = true;
					} catch (NullPointerException e) {
						return error("Log file entry code list not recognised: "+log);
					}
				}
			}
			
		} catch (Exception e) {
			e.printStackTrace();
			return error("Unknown error: "+e.getLocalizedMessage());
		}

		return 0;
	}

	private static int parseValue(String string) {
		if(string.isEmpty())
			return -1;
		int factor = 1;
		switch(Character.toUpperCase(string.charAt(string.length()-1))) {
		case 'K': factor = 1000; break;
		case 'M': factor = 1000000; break;
		}
		if(factor > 1)
			string = string.substring(0, string.length()-1);
		int result = -1;
		try {
			result = factor*Integer.parseInt(string);
		} catch (NumberFormatException e) {
		}
		return result;
	}

	private static int usage(MainManager man) {
		System.out.println(getUsageString(man));
		System.err.println("\nList of available substitution models:");
		for(String model : substModNames) {
			try {
				System.err.println("\t"+Class.forName(model).getSimpleName());
			} catch (Exception e) {}
		}
		return 1;
	}

	private static int error(String msg) {
		System.out.println("statalign: " + msg);
		return 2;
	}

	private static void warning(String msg) {
		System.out.println("warning: " + msg);
	}

	private static void initArrays(MainManager man) {
		findSubstMods();
		fillPostprocAbbr(man.postProcMan);
	}
	
	private static void fillPostprocAbbr(PostprocessManager man) {
		LinkedList<Integer> list = new LinkedList<Integer>();
		for(int i = 0; i < man.plugins.length; i++)
			list.add(i);
		
		for(int len = 1; list.size() > 0; len++) {
			for(ListIterator<Integer> it = list.listIterator(); it.hasNext(); ) {
				int pp = it.next();
				String str = man.plugins[pp].getTabName().substring(0, len).toUpperCase();
				Integer val;
				if((val = postprocAbbr.get(str)) == null) {		// empty slot
					postprocAbbr.put(str, pp);
					it.remove();			// pp is done
				} else if(val >= 0) {		// first collision
					postprocAbbr.put(str, -1);
//					it.add(val);			// pp is not done
				}
			}
		}
		
		for(String key : postprocAbbr.keySet()) {		// remove mappings for collisions
			if(postprocAbbr.get(key) == -1)
				postprocAbbr.remove(key);
		}
	}

	private static void findSubstMods() {
		for(String model : Utils.classesInPackage(SubstitutionModel.class.getPackage().getName()+".plugins")) {
			try {
				Class<?> cl = Class.forName(model);
				if (!Modifier.isAbstract(cl.getModifiers()) && SubstitutionModel.class.isAssignableFrom(cl))
					substModNames.add(model);
			} catch (Exception e) {		// handle class access exceptions etc.
			}
		}
	}
	
	private static String buildPpListStr(MainManager man, String linePrefix) {
		StringBuilder build = new StringBuilder();
		for(String key : postprocAbbr.keySet()) {
			build.append(linePrefix);
			build.append(key);
			build.append(": ");
			build.append(man.postProcMan.plugins[postprocAbbr.get(key)].getTip());
			build.append("\n");
		}
		return build.toString();
	}

	private static String buildDefPpList(MainManager man) {
		StringBuilder build = new StringBuilder();
		for(String key : postprocAbbr.keySet())
			if(man.postProcMan.plugins[postprocAbbr.get(key)].sampling) {
				build.append(key);
				build.append(",");
			}
		build.deleteCharAt(build.length()-1);
		return build.toString();
	}

}
