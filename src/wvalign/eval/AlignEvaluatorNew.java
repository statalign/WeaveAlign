package wvalign.eval;

import java.io.File;
import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.FileOutputStream;
import java.io.FilenameFilter;
import java.io.IOException;
import java.io.InputStream;
import java.io.OutputStream;
import java.io.PrintStream;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Comparator;
import java.util.List;
import java.util.regex.Matcher;
import java.util.regex.Pattern;
import java.util.Properties;

import wvalign.io.FastaAlignment;
import wvalign.io.FastaAlignmentReader;

public class AlignEvaluatorNew {

	private static final String[] requiredProps = { "referenceAlignment",
		"weaveAlignment", "dataDirOfBaseAligns", "outputCsv"};

	public static void main(String args[]) throws IOException {	
		
		if (args.length < 2) {
			System.out.println("Usage: AlignEvaluatorNew refAli.fasta testAli.fasta [output.file]");
			return;
		}
		
		String refAliFile = args[0];
		String testAliFile = args[1];
		String output = null;
		if (args.length>2) output = args[2];
		
		try {
			
		FileInputStream refStream = new FileInputStream(refAliFile);
		FileInputStream wvaStream = new FileInputStream(testAliFile);;
		OutputStream outStream;
		if (output==null) outStream = System.out;
		else outStream = new FileOutputStream(output);
		
		FastaAlignmentReader alignReader = new FastaAlignmentReader();

		FastaAlignment refAli = alignReader.readAlignment(refStream);			
		FastaAlignment testAli = alignReader.readAlignment(wvaStream);

		final PrintStream printStream = new PrintStream(outStream);
		printStream.print(refAli.getFsaScore(testAli)+"\n");
		printStream.close();

		} catch (FileNotFoundException e) {
			e.printStackTrace();
		} catch (IOException e) {
			e.printStackTrace();
		} catch (Exception e) {
			e.printStackTrace();
		}
	}
}
