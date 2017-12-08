package alignshow;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;

/**
 * Reads alignment annotation tracks (one number for alignment column).
 * Normalizes automatically when largest number exceeds 1.
 * Supports missing values, denoted by '*' (places <code>null</code>s in the list)
 * 
 * @author novak
 *
 */
public class TrackReader {

	private int errors;

	public List<Double> read(File file) throws IOException {
		return read(new BufferedReader(new FileReader(file)));
	}
	
	public List<Double> read(String fileName) throws IOException {
		return read(new BufferedReader(new FileReader(fileName)));
	}

	/**
	 * @return The annotation track read from {@link BufferedReader}, potentially with nulls in
	 * place of missing values, denoted by '*' in the stream
	 */
	public List<Double> read(BufferedReader reader) throws IOException {
		errors = 0;
		ArrayList<Double> values = new ArrayList<Double>();
		String line;
		double max = Double.MIN_VALUE;
		
		line = reader.readLine();
		while(line != null) {
			if(line.length() != 0)
				try {
					double val = Double.parseDouble(line);
					values.add(val);
					if(val > max)
						max = val;
				} catch (NumberFormatException e) {
					if(line.equals("*"))
						values.add(null);
					else
						errors++;
				}
			line = reader.readLine();
		}
		// if(max > 1)
		// 	for(int i = 0; i < values.size(); i++)
		// 		if(values.get(i) != null)
		// 			values.set(i, values.get(i)/max);
		return values;
	}
	
	/**
	 * Return the error count of the last read operation.
	 * 
	 * @return number of errors
	 */
	public int getErrors() {
		return errors;
	}
}
