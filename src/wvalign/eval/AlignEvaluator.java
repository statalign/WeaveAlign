package wvalign.eval;

import java.io.File;
import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.FileOutputStream;
import java.io.FilenameFilter;
import java.io.IOException;
import java.io.InputStream;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Comparator;
import java.util.List;
import java.util.regex.Matcher;
import java.util.regex.Pattern;
import java.util.Properties;

public class AlignEvaluator {

	private static final String[] requiredProps = { "referenceAlignment",
		"weaveAlignment", "dataDirOfBaseAligns", "outputCsv"};

	public static void main(String args[]) throws IOException {
		if (args.length != 1) {
			System.out.println("Usage: AlignEvaluator settings.properties");
			return;
		}
		PropertyHandler propHandler = new PropertyHandler(requiredProps);

		FileInputStream refStream;
		FileInputStream wvaStream;
		FileOutputStream outStream;
		
		try {
			Properties properties = propHandler.readPropertiesFile(args[0]);
			if (propHandler.checkProps(properties)) {

				AlignmentsEvaluator eval = new AlignmentsEvaluator();
				
				refStream = new FileInputStream(properties.getProperty("referenceAlignment"));
				wvaStream = new FileInputStream(properties.getProperty("weaveAlignment"));
				outStream = new FileOutputStream(properties.getProperty("outputCsv"));
				
				String path = properties.getProperty("dataDirOfBaseAligns");
				final String postfix = ".fsa";
				File folder = new File(path);
				
				File[] listOfFiles = folder.listFiles(new FilenameFilter() {
					public boolean accept(File dir, String name) {
						return name.toLowerCase().endsWith(postfix);
					}
				});

				Arrays.sort(
							listOfFiles,
							new Comparator<File>() {
								public int compare(File a, File b) {
									return a.getName().compareTo(b.getName());
								}
							});
				
				List<InputStream> baseInputs = new ArrayList<InputStream>();
				for (int i = 0; i < listOfFiles.length; i++) {
					if (listOfFiles[i].isFile()) {
						//System.out.println("File " + listOfFiles[i].getName());
						InputStream is = new FileInputStream(listOfFiles[i]);
						baseInputs.add(is);
					}
				}
				eval.setReferenceInput(refStream);
				eval.setBaseInputs(baseInputs);
				eval.setWvaInput(wvaStream);
				eval.setOutputStream(outStream);
				
				eval.evaluate();
				
				outStream.flush();
				outStream.close();	
			}
		} catch (FileNotFoundException e) {
			e.printStackTrace();
		} catch (IOException e) {
			e.printStackTrace();
		} catch (Exception e) {
			e.printStackTrace();
		}
	}
}
