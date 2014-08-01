package wvalign.io;

import java.io.BufferedReader;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;

public class FastaAlignmentReader {
	private BufferedReader reader;
	private String nextHeaderRow;

	public FastaAlignment readAlignment(InputStream inputStream1) 
			throws IOException {
		reader = new BufferedReader(new InputStreamReader(inputStream1));
		FastaAlignment ret = new FastaAlignment();
		String item = getNextFastaItem();
		while (item != null) {
			String[] fields = item.split("\\|");
			ret.addSeq(fields[0], fields[1]);
			item = getNextFastaItem();			
		}
		return ret;
	}

	private String getNextFastaItem() throws IOException {
		String line;
		String headerRow;
		boolean isFirst = true;
		String fastaItem = null;
		while ((line = reader.readLine()) != null) {
			if (line.charAt(0) == '>' && isFirst) {
				headerRow = line;
				fastaItem = new String(headerRow + "|");
				isFirst = false;
			} else if (line.charAt(0) != '>' && isFirst) {
				headerRow = nextHeaderRow;
				fastaItem = new String(headerRow + "|");
				fastaItem += line;
				isFirst = false;
			} else if (line.charAt(0) == '>' && !isFirst) {
				nextHeaderRow = line;
				return fastaItem;
			} else {
				fastaItem += line;
			}
		}
		return fastaItem;
	}

	public void closeReader() throws IOException {
		reader.close();
	}

}
