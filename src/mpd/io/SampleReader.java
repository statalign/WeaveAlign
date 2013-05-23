package mpd.io;

import java.io.BufferedReader;
import java.io.IOException;
import java.io.Reader;

public class SampleReader extends BufferedReader {
	
	private int currentSample;
	private String startString;
	
	private String line;
	private boolean eof;
	
	public SampleReader(Reader r) {
		super(r);
		makeString();
	}
	
	@Override
	public String readLine() throws IOException {
		if(eof)
			return null;
		if(line == null) {
			do
				line = super.readLine();
			while(line != null && line.indexOf("\tAlignment:\t") == -1);
			if(line == null) {
				eof = true;
				return null;
			}
		}
		if(line.startsWith(startString)) {
			String ret = line.substring(startString.length());
			line = null;
			return ret;
		}
		return null;
	}
	
	public int getCurrentSample() {
		return currentSample;
	}
	
	public boolean isEof() {
		return eof;
	}
	
	public void nextSample() {
		currentSample++;
		makeString();
	}
	
	private void makeString() {
		startString = String.format("Sample %d\tAlignment:\t", currentSample);
	}
}
