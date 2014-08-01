package wvalign.eval;


import java.io.IOException;
import java.io.InputStream;
import java.io.OutputStream;
import java.nio.charset.Charset;
import java.util.ArrayList;
import java.util.List;

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
		for (int i = 0; i< baseAlignments.size(); i++) {
			FastaAlignment testAlign = baseAlignments.get(i);
			writeBaseAlignScore(testAlign);
		}
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
