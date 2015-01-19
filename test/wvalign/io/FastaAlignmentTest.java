package wvalign.io;

import static org.junit.Assert.*;

import java.io.IOException;
import java.io.InputStream;
import java.util.HashSet;
import java.util.Set;

import org.apache.commons.lang3.tuple.ImmutablePair;
import org.apache.commons.lang3.tuple.Pair;
import org.junit.Before;
import org.junit.Ignore;
import org.junit.Test;

public class FastaAlignmentTest {

	InputStream inputStreamRef;
	InputStream inputStreamRef2;
	InputStream inputStream1;
	InputStream inputStream2;
	InputStream inputStream4;
	InputStream inputStream3;
	InputStream inputStreamSimple;
	InputStream inputStreamSimple2;
	FastaAlignmentReader reader;
	
	
	@Before
	public void setUp() {
		inputStreamRef = this.getClass().getResourceAsStream("/mytest-ref.fsa");
		inputStreamRef2 = this.getClass().getResourceAsStream("/mytest-ref.fsa");
		inputStream1 = this.getClass().getResourceAsStream("/mytest-test1.fsa");
		inputStream2 = this.getClass().getResourceAsStream("/mytest-test2.fsa");
		inputStream3 = this.getClass().getResourceAsStream("/mytest-test3.fsa");
		inputStream4 = this.getClass().getResourceAsStream("/mytest-test4.fsa");
		inputStreamSimple = this.getClass().getResourceAsStream("/simple.fsa");
		inputStreamSimple2 = this.getClass().getResourceAsStream("/simple2.fsa");
		reader = new FastaAlignmentReader();
	}
	
	@Test
	public void testSorted() throws IOException {
		FastaAlignment refAlignment = reader.readAlignment(inputStreamRef);
		refAlignment.sortByName();
		assertEquals("elme", refAlignment.get(0).id);
		assertEquals("meleg", refAlignment.get(2).id);
		FastaAlignment testAlignment = reader.readAlignment(inputStream1);
		testAlignment.sortByName();
		assertEquals("elme", testAlignment.get(0).id);
	}
	
	@Test
	public void testSameScolScore1() throws Exception {
		FastaAlignment refAlignment = reader.readAlignment(inputStreamRef);
		FastaAlignment testAlignment = reader.readAlignment(inputStreamRef2);
		assertEquals(1.0, refAlignment.getScolScore(testAlignment), 0.001);		
	}

	@Test
	public void testSameScolScore2() throws Exception {
		FastaAlignment refAlignment = reader.readAlignment(inputStreamRef);
		FastaAlignment testAlignment = reader.readAlignment(inputStream2);
		assertEquals(1.0, refAlignment.getScolScore(testAlignment), 0.001);		
	}

	@Test
	public void testZeroScolScore() throws Exception {
		FastaAlignment refAlignment = reader.readAlignment(inputStreamRef);
		FastaAlignment testAlignment = reader.readAlignment(inputStream1);
		assertEquals(0.0, refAlignment.getScolScore(testAlignment), 0.001);		
	}

	@Test
	public void testRealScolScore() throws Exception {
		FastaAlignment refAlignment = reader.readAlignment(inputStreamRef);
		FastaAlignment testAlignment = reader.readAlignment(inputStream4);
		assertEquals(0.5, refAlignment.getScolScore(testAlignment), 0.001);		
	}

	@Test
	public void testHomologyNum() throws Exception {
		FastaAlignment refAlignment = reader.readAlignment(inputStreamRef);
		assertEquals(27, refAlignment.getHomologNum());		
	}

	@Test
	public void testNonHomologyNum() throws Exception {
		FastaAlignment refAlignment = reader.readAlignment(inputStreamRef);
		assertEquals(12, refAlignment.getNonHomologNum());		
	}

	@Test
	public void testSimpleFsa() throws Exception {
		FastaAlignment refAlignment = reader.readAlignment(inputStreamSimple);
		FastaAlignment testAlignment = reader.readAlignment(inputStreamSimple2);
		assertEquals(0.667, refAlignment.getFsaScore(testAlignment), 0.001);		
	}

	
	@Test
	public void testFsaScore1() throws Exception {
		FastaAlignment refAlignment = reader.readAlignment(inputStreamRef);
		FastaAlignment testAlignment = reader.readAlignment(inputStream2);
		assertEquals(1.0, refAlignment.getFsaScore(testAlignment), 0.001);
	}
	
	@Test
	public void testFsaScore2() throws Exception {
		FastaAlignment refAlignment = reader.readAlignment(inputStreamRef);
		FastaAlignment testAlignment = reader.readAlignment(inputStream4);
		assertEquals(0.682, refAlignment.getFsaScore(testAlignment), 0.001);
	}

	@Test
	public void testFsaScore3() throws Exception {
		FastaAlignment refAlignment = reader.readAlignment(inputStreamRef);
		FastaAlignment testAlignment = reader.readAlignment(inputStream3);
		assertEquals(0.788, refAlignment.getFsaScore(testAlignment), 0.001);
	}

	@Test
	public void testPairSet() throws Exception {
		Set<Pair<Integer, Integer> > set = new HashSet<Pair<Integer, Integer> >();
		Pair<Integer, Integer> pair1 = new ImmutablePair<Integer, Integer>(1, 2);
		Pair<Integer, Integer> pair2 = new ImmutablePair<Integer, Integer>(1, 2);
		set.add(pair1);
		assertTrue(set.contains(pair2));
		set.add(pair2);
		assertEquals(1, set.size());
	}
}
