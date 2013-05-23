package mpd;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;

import ml.options.OptionData;
import ml.options.OptionSet;
import ml.options.Options;
import ml.options.Options.Multiplicity;
import ml.options.Options.Separator;
import mpd.io.ModReader;
import mpd.tree.Tree;
import statalign.model.subst.plugins.CustomSubstModel;


public class MpdMain {
	public static final String MPD_VERSION = "v2.1";

	private static final double DEF_G = 0;
	private static final String DEF_MPD_EXTENSION = ".mpd";
	private static final String DEF_FASTA_EXTENSION = ".fsa";
	private static final String DEF_ANNOT_EXTENSION = ".pred";
	private static final String DEF_SCORE_EXTENSION = ".scr";
	
	private static final String USAGE =
		"MPD "+MPD_VERSION+" (C) Adam Novak, 2010-12.\n\n"
		+
		"Usage:\n\n" +
		"    java -jar mpd.jar [options] input_1.log [-out output.ext] OR\n" +
		"    java -jar mpd.jar [options] input_1.fsa,... input_n.fsa [-out output.ext]\n\n"
		+
		"Description:\n\n" +
		"    Creates consensus alignment using the maximum posterior decoding (MPD)\n" +
		"    method from a set of alignments in FASTA format or from alignment samples\n" +
		"    in StatAlign's log format. Assigns reliability scores to each column.\n"+
		"    Default output file name is input_1.ext"+DEF_FASTA_EXTENSION+" for the alignment and\n" +
		"    input_1.ext"+DEF_SCORE_EXTENSION+" for the column scores. When -out is used, score output\n" +
		"    file is named output"+DEF_SCORE_EXTENSION+"\n"+
		"\n\n"
		+
		"Options:\n\n" +
		"    -g=VAL\n" +
		"        Sets the g parameter's value in the MPD algorithm.\n" +
		"        Default: "+DEF_G+"\n\n"
		+
//		disabled - these apply to MpdInterfaceTree only
//		"    -cm\n" +
//		"        Creates files input.log.fwd and input.log.bwd with forward/backward\n" +
//		"        conditional marginals for the MPD alignment's columns.\n\n"
//		+
//		"    -t treefile\n" +
//		"        For each column of the MPD alignment finds the split within the tree\n" +
//		"        where the fwd/bwd conditional marginal of the subsequences is\n" +
//		"        maximal and outputs to input.log.split.fwd and .bwd\n\n"
//		+
		"    -n=SAMPNO\n"+
		"        Limits the number of samples read and used to construct the MPD\n"+
		"        alignment to SAMPNO. Default: all available\n\n"
		+
		"    -r=SAMPRATE\n"+
		"        Selects a subset of the samples read, with rate SAMPRATE, to decrease\n"+
		"        autocorrelation effect. Default: 1\n\n"
		+
		"    -mpdout\n"+
		"        Instead of separate FASTA and scores files, creates a single MPD file\n" +
		"        that contains both. Default name is input_1.ext"+DEF_MPD_EXTENSION+"\n\n"
		+
		"    -optgi\n"+
		"        Optimise using the gap insensitive scoring of the columns\n" +
		"        Default: score by the gap sensitive coding which the DAG is based on\n\n"
		+
		"    -outgi\n" +
		"        Output the gap insensitive column scores along with the default scores\n" +
		"        in the second column of the score output (mpd file or score file)\n\n"
		+
		"  Annotation options:\n\n"
		+
		"    -mod modfile\n"+
		"        Reads mod file that defines a model for annotation. You can specify\n" +
		"        multiple mod files (one per annotation state) or alternatively just\n" +
		"        a single one and edge length multipliers (rhos) for the rest of the\n" +
		"        states\n\n"
		+
		"    -rho=RHO\n"+
		"        Specifies an edge length multiplier for annotation model. If one mod\n"+
		"        file is used then at least one rho must be specified.\n\n"
		+
		"    -hmm trprobfile\n"+
		"        Reads annotation HMM transition probability matrix from file. Steady\n" +
		"        state probabilities are assumed to follow the matrix in the file.\n\n"
		+
		"    -gr=INS,DEL\n"+
		"        Substitution models of all states are extended with a gap character\n"+
		"        using the specified average insertion and deletion rate. Default: Gaps\n"+
		"        are treated as missing data and averaged over.\n\n"
		+
		"    -mpd=MODE\n"+
		"        Defines the way the MPD alignment is created (the annotation for which\n"+
		"        is then printed in the predfile). Mode 1: based on column probabilities\n"+
		"        Mode 2 (default): column and annotation joint probabilities.\n\n"
		+
		"    -pred predfile\n"+
		"        Specifies the output file for annotation probabilities. Rows in the\n"+
		"        file correspond to MPD alignment columns, columns to annotation states.\n"+
		"        Default: input_1.ext"+DEF_ANNOT_EXTENSION+"\n\n"
		+
		"    -pseq=SEQID\n"+
		"        Annotation probabilites are projected down to the given sequence and\n"+
		"        listed in the file input_1.ext"+DEF_ANNOT_EXTENSION+"1\n\n";

	private static final boolean scoreSamples = false;

	public static double lastDataProb;
		
	public static void main(String[] args) {
//		run();
//		System.exit(0);
		
		Options opt = new Options(args, Multiplicity.ZERO_OR_ONE, 1, Integer.MAX_VALUE);
		opt.addSet("run")
				.addOption("out", Separator.BLANK)
				.addOption("g", Separator.EQUALS)
//				.addOption("cm")
//				.addOption("t", Separator.BLANK)
				.addOption("n", Separator.EQUALS)
				.addOption("r", Separator.EQUALS)
				.addOption("mpdout")
				.addOption("optgi")
				.addOption("outgi")
				.addOption("mod", Separator.BLANK, Multiplicity.ZERO_OR_MORE)
				.addOption("rho", Separator.EQUALS, Multiplicity.ZERO_OR_MORE)
				.addOption("hmm", Separator.BLANK)
				.addOption("gr", Separator.EQUALS)
				.addOption("mpd", Separator.EQUALS)
				.addOption("pred", Separator.BLANK)
				.addOption("pseq", Separator.EQUALS);
		
		OptionSet set;
		if ((set = opt.getMatchingSet(false, false)) == null) {
			System.out.println(USAGE);
			System.exit(1);
		}

		ArrayList<String> data = set.getData();
		String input0 = data.get(0);

		for(String input : data)
			if(!new File(input).exists())
				error("input file '"+input+"' does not exist.");
		
		try {
			boolean mpdOut = set.isSet("mpdout");
			String output;
			if(set.isSet("out"))
				output = set.getOption("out").getResultValue(0);
			else
				output = input0 + (mpdOut?DEF_MPD_EXTENSION:DEF_FASTA_EXTENSION);
			
			String scoreOutput = output;
			if(!mpdOut) {
				int pos = scoreOutput.lastIndexOf('.');
				if(pos == -1)
					pos = scoreOutput.length();
				scoreOutput = scoreOutput.substring(0, pos).concat(DEF_SCORE_EXTENSION);
			}
			
			double g = DEF_G;
			if(set.isSet("g")) {
				String gOpt = set.getOption("g").getResultValue(0);
				try {
					g = Double.parseDouble(gOpt);
				} catch (NumberFormatException e) {
					error("wrong number format for parameter g: "+gOpt);
				}
			}
			
			MpdInterface mpdIf = new MpdInterface(g, set.isSet("optgi"), set.isSet("outgi"));
			
//			if(set.isSet("cm")) {
//				mpdIf.setCondMargFiles(input+".fwd", input+".bwd");
//			}
//			
//			if(set.isSet("t")) {
//				String treeFile = set.getOption("t").getResultValue(0);
//				List<List<Integer>> splits = TreeSplits.getSplits(treeFile);
//				mpdIf.addSplits(splits, input+".split.fwd", input+".split.bwd");
//			}
			
			if(set.isSet("n")) {
				int value = Integer.parseInt(set.getOption("n").getResultValue(0));
				mpdIf.setMaxNoSamples(value);
			}
			
			if(set.isSet("r")) {
				int value = Integer.parseInt(set.getOption("r").getResultValue(0));
				mpdIf.setSampleRate(value);
			}
			
			if(set.isSet("mod")) {	// set up annotator
				OptionData option = set.getOption("mod");
				
				ModReader r = new ModReader();
				ArrayList<CustomSubstModel> models = new ArrayList<CustomSubstModel>();
				for(int i = 0; i < option.getResultCount(); i++) {
					r.readFile(new File(option.getResultValue(i)));
					models.add(r.getModel());
				}
				Tree tree = r.getTree();
				
				ArrayList<Double> rhos = null;
				if(set.isSet("rho")) {
					if(models.size() > 1)
						error("-rho can only be used with a single mod file");
					rhos = new ArrayList<Double>();
					rhos.add(1.0);
					option = set.getOption("rho");
					for(int i = 0; i < option.getResultCount(); i++) {
						try {
							rhos.add(Double.parseDouble(option.getResultValue(i)));
						} catch (NumberFormatException e) {
							error("wrong number format for parameter rho: "+option.getResultValue(i));
						}
					}
				} else if(models.size() == 1) {
					error("at least two models (-mod) or one -rho parameter is required!");
				}
				
				if(!set.isSet("hmm"))
					error("-hmm is not specified and is required for annotation");
				BufferedReader br = new BufferedReader(new FileReader(set.getOption("hmm").getResultValue(0)));
				double[][] transMat = ModReader.readMatrix(br);
				String line;
				do {
					line = br.readLine().trim();
				} while(line.isEmpty());
				double[] initState = ModReader.getDblVector(line);
				
				if(set.isSet("gr")) {
					String val = set.getOption("gr").getResultValue(0);
					String[] vals = val.split(",");
					try {
						if(vals.length != 2)
							throw new NumberFormatException();

						double ins = Double.parseDouble(vals[0]);
						double del = Double.parseDouble(vals[1]);
						for(int i = 0; i < models.size(); i++)
							models.get(i).addGapChar(ins, del);
						
					} catch (NumberFormatException e) {
						error("wrong format for gap rates: "+val);
					}
				}
				
				MpdAnnotator annotator;
				if(rhos == null)
					annotator = new MpdAnnotator(tree, transMat, initState, models);
				else
					annotator = new MpdAnnotator(tree, transMat, initState, models.get(0), rhos);
				
				if(set.isSet("mpd")) {
					String val = set.getOption("mpd").getResultValue(0);
					try {
						int mode = Integer.parseInt(val);
						if(mode < 1 || mode > 2)
							throw new NumberFormatException();
						annotator.setMpdMode(mode-1);
					} catch (NumberFormatException e) {
						error("wrong number for mpd parameter "+ val);
					}
				}
				
				String annOut;
				if(set.isSet("pred")) {
					annOut = set.getOption("pred").getResultValue(0);
				} else {
					annOut = input0 + DEF_ANNOT_EXTENSION;
				}
				annotator.setOutFile(new File(annOut));
				
				if(set.isSet("pseq")) {
					String pseq = set.getOption("pseq").getResultValue(0);
					if(!tree.getNames().contains(pseq)) {
						error("could not find sequence id "+pseq);
					}
					annotator.setProjectSeq(pseq, new File(input0+DEF_ANNOT_EXTENSION+"1"));
				}
				
				mpdIf.setAnnotator(annotator);
				
			} else if(set.isSet("rho") || set.isSet("hmm") || set.isSet("gr") || set.isSet("mpd")
					|| set.isSet("pred") || set.isSet("pseq")) {
				error("-rho, -hmm, -gr, -mpd, -pred and -pseq can only be used together with -mod");
			}
		
//			System.out.println(input+" "+output+" "+g);
			if(!scoreSamples) {
				mpdIf.doMpd(data, output, scoreOutput, 0);
			} else {
				mpdIf.doMpd(data, output, scoreOutput, 1);
				mpdIf.doMpd(data, output, scoreOutput, 2);
			}
			if(mpdIf.getAnnotator() != null)
				lastDataProb = mpdIf.getAnnotator().getDataProb();
		} catch (IOException e) {
			e.printStackTrace();
			System.exit(1);
		}
	}

	static void run() {
		String input = "files/families/1sbp_subfamilies.fasta.log";
		String output = input+".mpd";
		
		try {
			MpdInterface mpdIf = new MpdInterface(0, false, false);
			if(!scoreSamples) {
				mpdIf.doMpd(input, output, output, 0);
			} else {
				mpdIf.doMpd(input, output, output, 1);
				mpdIf.doMpd(input, output, output, 2);
			}
		} catch (Exception e) {
			e.printStackTrace();
			System.exit(1);
		}
	}

	private static void error(String msg) {
		System.out.println("mpd: " + msg);
		System.exit(1);
	}

}
