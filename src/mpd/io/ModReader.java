package mpd.io;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;

import mpd.io.NewickReader.FormatException;
import mpd.model.CustomSubstModel;
import mpd.tree.Tree;

public class ModReader {
	
	CustomSubstModel model;
	Tree tree;
	
	public ModReader() {
	}
	
	public void readFile(File file) throws IOException {
		BufferedReader r = new BufferedReader(new FileReader(file));
		double[][] rateMat = null;
		double[] background = null;
		String alphabet = null;

		String line;
		while((line = r.readLine()) != null) {
			if(line.matches("^ALPHABET:.+")) {
				alphabet = line.split(":")[1].replaceAll("\\s+", "");
			} else if(line.matches("^BACKGROUND:.+")) {
				background = getDblVector(line.split(":")[1].trim());
			} else if(line.matches("^RATE_MAT:.*")) {
				rateMat = readMatrix(r);
			} else if(line.matches("^TREE:.+")) {
				line = line.replaceAll("TREE:", "");
				NewickReader tr = new NewickReader(line, 0);
				try {
					tree = tr.parseTree();
				} catch (FormatException e) {
					throw new Error("ModReader: tree parse exception: "+e);
				}
			}
		}
		// consistency check
		if(background == null || alphabet == null || rateMat == null)
			throw new Error("ModReader: missing one of the required fields BACKGROUND, ALPHABET or RATE_MAT");
		int l = background.length;
		if(alphabet.length() != l || rateMat.length != l)
			throw new Error("ModReader: ALPHABET of RATE_MAT inconsistent with BACKGROUND!");
		model = new CustomSubstModel(alphabet, rateMat, background);
	}
	
	public Tree getTree() {
		return tree;
	}
	
	public CustomSubstModel getModel() {
		return model;
	}
	
	/**
	 * Reads a square matrix.
	 * @throws Error if matrix rows are not of consistent length
	 * @throws IOException on I/O error
	 */
	public static double[][] readMatrix(BufferedReader reader) throws IOException {
		double[][] rateMat = null;
		double[] row = null;
		int i = 0;
		do {
			String line = reader.readLine();
			if(line == null)
				throw new Error("ModReader: unexpected end of file while reading matrix");
			line = line.trim();
			if(line.isEmpty())
				continue;
			row = getDblVector(line);
			if(rateMat == null)
				rateMat = new double[row.length][];
			else if(row.length != rateMat[0].length)
				throw new Error("ModReader: number of columns inconsistent in matrix");
			rateMat[i++] = row;
		} while(row == null || i < row.length);
		return rateMat;
	}
	
	public static double[] getDblVector(String str) {
		String[] strs = str.split("\\s+");
		double[] dbls = new double[strs.length];
		for(int i = 0; i < strs.length; i++)
			dbls[i] = Double.parseDouble(strs[i]);
		return dbls;
	}
}
