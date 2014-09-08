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
	}

	private Percentile computeScoreStats(String type) throws Exception {
		int alignmentNum = baseAlignments.size();		
		Percentile percentile = new Percentile();
		double[] scores = new double[alignmentNum];
		for (int i = 0; i < alignmentNum; i++) {
			FastaAlignment testAlign = baseAlignments.get(i);			
			if (type.equals("fsa")) {
				scores[i] = refAlignment.getFsaScore(testAlign);
			}
			else if (type.equals("scol")) {
				
			}
			else throw new RuntimeException("Score type `"+type+"' not recognised.");
		}
		percentile.setData(scores);
		return percentile;
	}

	private double computeWvaRank(double wvaScore, double[] scores) {
		Arrays.sort(scores);		
		int rank = 0;
		while (rank < scores.length && wvaScore > scores[rank]) rank++;
		return (100.0 * rank) / scores.length;
	}

	private void printScoreStatsToScreen() throws Exception {
		double wvaFsaScore = refAlignment.getFsaScore(wvaAlignment);
//		double wvaScolScore = refAlignment.getScolScore(wvaAlignment);
		
		System.out.printf("WeaveAlign sum-of-pairs (AMA) score = %5.3f\n", wvaFsaScore);
		//System.out.printf("WeaveAlign column score = %5.3f\n", wvaScolScore);
		
		Percentile percentileFsa = computeScoreStats("fsa");		
		System.out.println("AMA scores for inputted alignment samples:");
		System.out.printf("%-10s%12s%12s%12s\n","Percentile","5%","50%","95%");
		System.out.printf("%-10s%12.3f%12.3f%12.3f\n","Score",percentileFsa.evaluate(5.0),percentileFsa.evaluate(50.0),percentileFsa.evaluate(95.0));
//		Percentile percentileScol = computeScoreStats("scol");
//		System.out.println("Column scores for inputted alignment samples:");
//		System.out.printf("%12s%12s%12s\n","5%","50%","95%");
//		System.out.printf("%12.3f%12.3f%12.3f\n",percentileScol.evaluate(5.0),percentileScol.evaluate(50.0),percentileScol.evaluate(95.0));
		
		printSummaryToScreen(wvaFsaScore,percentileFsa.getData(),"sum-of-pairs");
	}

	private void printSummaryToScreen(double wvaScore, double[] scores,String scoreName) {
		double rankScore = computeWvaRank(wvaScore,scores);		
		//double wvaPercentile = 100.0 * (1.0 - (0.1* rank / inputAlignNum));
		System.out.print("\nRank score = "+rankScore);
		System.out.println(" (WeaveAlign's alignment scores better than " + rankScore 
				+ "% of the input alignments, using the "+ scoreName + " score)");
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
