package wvalign.eval;

import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.OutputStream;


public class AlignSetEval {
	
	public static void main(String args[]) throws IOException {	
		
		if (args.length < 2) {
			System.out.println("Usage: AlignSetEval refAli.fasta wvaSummary.fasta aliSamples.log outputFile [firstSample] [lastSample]");
			return;
		}
		
		String refAliFile = args[0];
		String wvaAliFile = args[1];
		String logFile = args[2];
		String output = args[3];
		int first = 0;
		int last = Integer.MAX_VALUE;
		
		if (args.length>4) first = Integer.parseInt(args[4]);
		if (args.length>5) last = Integer.parseInt(args[5]);
				
		try {
			
		FileInputStream refStream = new FileInputStream(refAliFile);		
		FileInputStream wvaStream = new FileInputStream(wvaAliFile);
		OutputStream outStream = new FileOutputStream(output);
		
		AlignmentsEvaluator eval = new AlignmentsEvaluator();
		
		eval.setReferenceInput(refStream);
		eval.readLogFile(logFile,first,last);
		eval.setWvaInput(wvaStream);
		eval.setOutputStream(outStream);
		
		boolean printScoresToFile = false;
		eval.evaluate(printScoresToFile);
		
		outStream.flush();
		outStream.close();	
		
		
		} catch (FileNotFoundException e) {
			e.printStackTrace();
		} catch (IOException e) {
			e.printStackTrace();
		} catch (Exception e) {
			e.printStackTrace();
		}
	}
}
