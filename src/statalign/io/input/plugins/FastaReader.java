package statalign.io.input.plugins;

import java.io.BufferedReader;
import java.io.IOException;

import statalign.io.RawSequences;
import statalign.io.input.FileFormatReader;

/**
 * 
 * Class to read files in Fasta format.
 * 
 * @author novak
 *
 */
public class FastaReader extends FileFormatReader {

	private int errors;

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
		errors = 0;
		while(true) {
			line = reader.readLine();
			if(line != null && line.length() == 0)
				continue;
			if(line == null || line.charAt(0) == '>') {
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
				actSeq.setLength(0);
				inSeq = true;
				line = line.substring(1).trim();	// remove leading > and leading/trailing whitespace
//				line = line.replaceAll("[ \t]+", "_");
//				line = line.replaceAll("\\(", "{");
//				line = line.replaceAll("\\)", "}");
				result.seqNames.add(line);
			} else if(inSeq) {
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
			} else {
				errors++;
			}
		}

//		for(int i = 0; i < result.sequences.size(); i++){
//			System.out.println(">"+result.seqNames.get(i)+"\n"+result.sequences.get(i));
//		}
//		if(errors > 0)
//			System.out.println("Errors: "+errors);
//		else
//			System.out.println("FastaReader: successfully read "+result.sequences.size()+" sequences.");
		
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

	/**
	 * Only for testing/debugging purposes: reads the given Fasta file
	 * 
	 * @param args First element must be the name of the Fasta file to read
	 * @throws IOException if an I/O error occurs
	 */
	public static void main(String[] args) throws IOException {
		FastaReader f = new FastaReader();
		f.read(args[0]);
	}
}