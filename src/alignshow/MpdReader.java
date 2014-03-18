package alignshow;

import java.io.BufferedReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;

import wvalign.io.FileFormatReader;
import wvalign.io.RawSequences;



/**
 * 
 * Class to read files in Fasta format.
 * 
 * @author novak
 *
 */
public class MpdReader extends FileFormatReader {

	private static final String SCORES_SEPARATOR = "#scores";
	private int errors;
	private List<Double> scores;

	/**
	 * Reads the contents (aligned/non-aligned sequences) of the given data source in
	 * Fasta format.
	 * 
	 * @param reader Data source
	 * @return RawSequences representation of the contents
	 * @throws IOException if an I/O error occurs
	 */
	@Override
	public RawSequences read(BufferedReader reader) throws IOException {
		RawSequences result = new RawSequences();

		String line;
		boolean inSeq = false;
		StringBuilder actSeq = new StringBuilder();
		int errors = 0;
		while(true) {
			line = reader.readLine();
			if(line != null && line.length() == 0)
				continue;
			if(line == null || line.startsWith(SCORES_SEPARATOR) || line.charAt(0) == '>') {
				if(inSeq) {
					if(actSeq.length() == 0) {
						errors++;
						result.seqNames.remove(result.seqNames.size()-1);
					} else {
						result.sequences.add(actSeq.toString());
					}
				}
				if(line == null)
					break;
				if(line.startsWith(SCORES_SEPARATOR)) {
					scores = new ArrayList<Double>();
					for(;;) {
						line = reader.readLine();
						if(line == null)
							break;
						if(line.length() == 0)
							continue;
						try {
							scores.add(Double.parseDouble(line));
						} catch (NumberFormatException e) {
							if(line.equals("*"))
								scores.add(null);	// missing value
							else
								errors++;
						}
					}
					if(result.sequences.size() == 0 || result.sequences.get(0).length() != scores.size())
						errors++;
					break;
				}
				actSeq.setLength(0);
				inSeq = true;
				int start = 1, index;
				if((index = line.indexOf(' ', 1)) == 1) {
					index = line.indexOf(' ', 2);
					start = 2;
				}
				if(index == -1)
					line = line.substring(start);
				else{
					//line = line.substring(start, index);
					line = line.substring(start);			
				}
				line = line.replaceAll("[ \t]+", "_");
				//line.replaceAll(" ", "_");
				line = line.replaceAll("\\(", "{");
				line = line.replaceAll("\\)", "}");
//				System.out.println("new name: "+line);
				result.seqNames.add(line);
			} else if(line.indexOf('\t') != -1 || inSeq) {
				if(!inSeq) {
					String[] parts = line.split("\t");
					if(parts.length != 2) {
						errors++;
					}
					result.seqNames.add(parts[0].trim());
					line = parts[1];
					actSeq.setLength(0);
				}
				int len = line.length();
				char ch;
				for(int i = 0; i < len; i++) {
					if(!Character.isWhitespace(ch = line.charAt(i))) {
						if(Character.isLetter(ch) || ch == '-')
							actSeq.append(ch);
						else
							errors++;
					}
				}
				if(!inSeq) {
					result.sequences.add(actSeq.toString());
				}
			} else {
				errors++;
			}
		}

		return result;
	}
	
	/**
	 * Return the error count of the last read operation.
	 * 
	 * @return number of errors
	 */
	public int getErrors() {
		return errors;
	}
	
	public List<Double> getScores() {
		return scores;
	}

}