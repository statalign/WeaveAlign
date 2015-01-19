package wvalign.io;

import java.io.BufferedReader;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.util.ArrayList;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

import org.apache.commons.lang3.StringUtils;

public class FastaAlignmentReader {
	private BufferedReader reader;
	private Pattern newSample = Pattern.compile("Sample (\\d+)\\s+Alignment:\\s+");
	
	public FastaAlignment readAlignment(InputStream inputStream1) 
			throws IOException {
		reader = new BufferedReader(new InputStreamReader(inputStream1));
		FastaAlignment ret = new FastaAlignment();		
		String[] fields = new String[3];
		String line;
		while ( (line=reader.readLine()) != null) {
			if (line.charAt(0) == '>') {
				//System.out.println(line);
				if (fields.length < 3) ret.addSeq(fields[0], fields[1]);				
				fields = new String[2];
				fields[0] = line.substring(1);
				fields[1] = "";
				continue;
			}
			else {
				fields[1] += StringUtils.chomp(line); 
			}								
		}
		if (fields.length < 3) ret.addSeq(fields[0], fields[1]);
		//ret.sortByName();
		return ret;
	}
	
	public ArrayList<FastaAlignment> readLogFile(InputStream inputStream1)
			throws IOException {
		return readLogFile(inputStream1, 0, Integer.MAX_VALUE); 
	}

	public ArrayList<FastaAlignment> readLogFile(InputStream inputStream1, int firstSample, int lastSample) 
			throws IOException {
		reader = new BufferedReader(new InputStreamReader(inputStream1));
		ArrayList<FastaAlignment> alignments = new ArrayList<FastaAlignment>();				
		int sample = firstSample;
		FastaAlignment currentAlignment = new FastaAlignment();
		String[] fields = new String[3]; // when length 3, indicates not initialised
		String item;
		while ( (item=reader.readLine()) != null) {						
			Matcher lineMatch = newSample.matcher(item);			
			if (lineMatch.find()) {
				int nextSample = Integer.parseInt(lineMatch.group(1));
				if (nextSample < firstSample) continue;				
				if (nextSample > lastSample) break;				
				if (nextSample != sample) {
					//System.out.println("Sample "+nextSample);	
					if (fields.length < 3) currentAlignment.addSeq(fields[0], fields[1]);				
					alignments.add(currentAlignment);
					fields = new String[3];
					currentAlignment = new FastaAlignment();
					sample = nextSample;
				}				
				String line = lineMatch.replaceAll("");
				if (line.charAt(0) == '>') {
					//System.out.println(line);
					if (fields.length < 3) currentAlignment.addSeq(fields[0], fields[1]);				
					fields = new String[2]; // length 2 indicates initialised
					fields[0] = line.substring(1);
					fields[1] = "";
					continue;
				}
				else {
					fields[1] += StringUtils.chomp(line); 
				}						
			}			
			else throw new IOException("Incorrect logfile format.");
		}
		if (fields.length < 3) currentAlignment.addSeq(fields[0], fields[1]);
		alignments.add(currentAlignment);
		return alignments;
	}

	public void closeReader() throws IOException {
		reader.close();
	}

}
