package wvalign.eval;

import static org.junit.Assert.*;

import java.io.ByteArrayInputStream;
import java.io.ByteArrayOutputStream;
import java.io.InputStream;
import java.io.OutputStream;
import java.util.ArrayList;
import java.util.List;

import org.junit.Before;
import org.junit.Test;


public class AlignmentsEvaluatorTest {

	final String refAlignment = ">emelet\n" + 
			"EMELET--\n" + 
			">elme\n" + 
			"--ELME--\n" + 
			">mellett\n" + 
			"-MELLETT";
	final String wvaAlignment = ">emelet\n" + 
			"EMELET--\n" + 
			">elme\n" + 
			"--ELME--\n" + 
			">mellett\n" + 
			"-MELLETT";
	final String baseAlignment1 = ">elme\n" + 
			"ELME-----\n" + 
			">emelet\n" + 
			"E-MEL--ET\n" + 
			">mellett\n" + 
			"--MELLETT";
	final String baseAlignment2 = ">elme\n" + 
			"ELME-----\n" + 
			">mellett\n" + 
			"--MELLETT\n" +
			">emelet\n" + 
			"E-MEL--ET";
	final String baseAlignment3 = ">emelet\n" + 
			"EMEL--ET\n" + 
			">elme\n" + 
			"--ELME--\n" + 
			">mellett\n" + 
			"-MELLETT";
	private InputStream isRef; 
	private InputStream isBase;
	private InputStream isWva;
	private List<InputStream> baseInputs = new ArrayList<InputStream>();
	
	@Before
	public void setUp() {
	    isRef = new ByteArrayInputStream(refAlignment.getBytes());
	    isWva = new ByteArrayInputStream(wvaAlignment.getBytes());
	    InputStream isBase1 = new ByteArrayInputStream(baseAlignment1.getBytes());
	    InputStream isBase2 = new ByteArrayInputStream(baseAlignment2.getBytes());
	    InputStream isBase3 = new ByteArrayInputStream(baseAlignment3.getBytes());
	    baseInputs.add(isBase1);
	    baseInputs.add(isBase2);
	    baseInputs.add(isBase3);
	}

	@Test
	public void test() throws Exception {
		AlignmentsEvaluator eval = new AlignmentsEvaluator();
		eval.setReferenceInput(isRef);
		eval.setWvaInput(isWva);
		eval.setBaseInputs(baseInputs);
		OutputStream os = new ByteArrayOutputStream();
		eval.setOutputStream(os);
		eval.evaluate();
		String[] rows = os.toString().split("\n");
		assertEquals("alignlen,fsascore,scolscore,alignid", rows[0]);
		assertEquals("8,1.0,1.0", rows[1]);
		
		String[] scores1 = rows[2].split(",");
		double fsa1 = Double.parseDouble(scores1[1]);
		double col1 = Double.parseDouble(scores1[2]);
		assertEquals(0.265, fsa1, 0.001);
		assertEquals(0.0, col1, 0.001);
	
		String[] scores2 = rows[3].split(",");
		double fsa2 = Double.parseDouble(scores2[1]);
		double col2 = Double.parseDouble(scores2[2]);
		assertEquals(0.265, fsa2, 0.001);
		assertEquals(0.0, col2, 0.001);

		String[] scores3 = rows[4].split(",");
		double fsa3 = Double.parseDouble(scores3[1]);
		double col3 = Double.parseDouble(scores3[2]);
		assertEquals(0.706, fsa3, 0.001);
		assertEquals(0.5, col3, 0.001);
	}

}
