package wvalign.eval;


import java.io.IOException;
import java.io.InputStream;
import java.io.OutputStream;
import java.nio.charset.Charset;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

import org.apache.commons.math3.stat.descriptive.rank.Percentile;

import wvalign.io.FastaAlignment;
import wvalign.io.FastaAlignmentReader;

public class AlignmentsEvaluator {

	private List<FastaAlignment> baseAlignments;
	private FastaAlignment wvaAlignment;
	private FastaAlignment refAlignment;
	private OutputStream output;
	private final String sep = ",";
	private Percentile percentileFsa;
	private double[] fsaScores;

	private FastaAlignmentReader alignReader = new FastaAlignmentReader();

	public void setReferenceInput(InputStream refStream) throws IOException {
		refAlignment = alignReader.readAlignment(refStream);
		System.out.println("Ref alignment size: " + refAlignment.getNumOfSequences() + 
				"  width: " + refAlignment.getWidth());
	}

	public void setWvaInput(InputStream wvaStream) throws IOException {
		wvaAlignment = alignReader.readAlignment(wvaStream);
	}

	public void setBaseInputs(List<InputStream> baseInputs) throws IOException {
		baseAlignments = new ArrayList<FastaAlignment>();
		for (InputStream is : baseInputs) {
			FastaAlignment al = alignReader.readAlignment(is);
			baseAlignments.add(al);
		}
	}

	public void setOutputStream(OutputStream outStream) {
		output = outStream;
	}

	public void evaluate() throws Exception {
		writeHeader();
		writeWvaScores();
		for (int i = 0; i < baseAlignments.size(); i++) {
			FastaAlignment testAlign = baseAlignments.get(i);
			writeBaseAlignScore(testAlign);
			
		}
		printScoreStatsToScreen();
		printSummaryToScreen();
	}

	private void computeFsaScoreStats() throws Exception {
		int alignmentNum = baseAlignments.size();
		percentileFsa = new Percentile();
		fsaScores = new double[alignmentNum];
		for (int i = 0; i < alignmentNum; i++) {
			FastaAlignment testAlign = baseAlignments.get(i);
			double fsaScore = refAlignment.getFsaScore(testAlign);
			fsaScores[i] = fsaScore;
		}
		percentileFsa.setData(fsaScores);
	}

	private int computeWvaFsaRank() {
		Arrays.sort(fsaScores);
		double targetScore = refAlignment.getFsaScore(wvaAlignment);
		int i = 0;
		while (i < fsaScores.length && targetScore > fsaScores[i]) {
			i++;
		}
		// how many of the samples were better? 0 is the best
		int wvaFsaRank = i;
		return wvaFsaRank;
	}

	private void printScoreStatsToScreen() throws Exception {
		computeFsaScoreStats();
		System.out.println("WVA FSA score = " + refAlignment.getFsaScore(wvaAlignment));
		System.out.println("WVA SCOL score = " + refAlignment.getScolScore(wvaAlignment));
		System.out.println("FSA scores 5 percentile = " + percentileFsa.evaluate(5.0));
		System.out.println("FSA scores 50 percentile = " + percentileFsa.evaluate(50.0));
		System.out.println("FSA scores 95 percentile = " + percentileFsa.evaluate(95.0));
		System.out.println("SCOL scores 5 percentile = " + percentileFsa.evaluate(5.0));
		System.out.println("SCOL scores 50 percentile = " + percentileFsa.evaluate(50.0));
		System.out.println("SCOL scores 95 percentile = " + percentileFsa.evaluate(90.0));
		
	}

	private void printSummaryToScreen() {
		int rank = computeWvaFsaRank();
		int inputAlignNum = fsaScores.length;
		double wvaPercentile = 100.0 * (1.0 - (0.1* rank / inputAlignNum));
		System.out.println("\n------------------------");
		System.out.println("WeaveAlign's alignment scores better than " +
				+ wvaPercentile
				+ " percent of the input alignments (FSA score).");
		System.out.println("------------------------");
	}

	private void writeBaseAlignScore(FastaAlignment testAlign) throws Exception {
		output.write(getScoresRow(testAlign).getBytes(Charset.forName("UTF-8")));		
	}

	private void writeHeader() throws IOException {
		StringBuilder sb = new StringBuilder();
		sb.append("alignlen,fsascore,scolscore,alignid\n");
		output.write(sb.toString().getBytes(Charset.forName("UTF-8")));
	}

	private void writeWvaScores() throws Exception {
		String wvaScores = getScoresRow(wvaAlignment);
		output.write(wvaScores.getBytes(Charset.forName("UTF-8")));
	}

	private String getScoresRow(FastaAlignment testAlign) throws Exception {
		StringBuilder sb = new StringBuilder();
		double scolScore = refAlignment.getScolScore(testAlign);
		double fsaScore = refAlignment.getFsaScore(testAlign);
		sb.append(testAlign.getWidth());
		sb.append(sep);
		sb.append(fsaScore);
		sb.append(sep);
		sb.append(scolScore);
		sb.append("\n");
		return sb.toString();
	}
	
}
