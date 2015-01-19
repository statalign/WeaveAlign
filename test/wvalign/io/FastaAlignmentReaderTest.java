package wvalign.io;

import static org.junit.Assert.*;

import java.io.IOException;
import java.io.InputStream;

import org.junit.Before;
import org.junit.Test;
import static wvalign.io.GappedSequence.GAP;

public class FastaAlignmentReaderTest {
	InputStream inputStream1;
	InputStream inputStream2;
	FastaAlignmentReader reader;
	
	
	@Before
	public void setUp() {
		inputStream1 = this.getClass().getResourceAsStream("/mytest-ref.fsa");
		inputStream2 = this.getClass().getResourceAsStream("/mytest-test4.fsa");
		reader = new FastaAlignmentReader();
	}

	@Test
	public void testNotNull() {
		FastaAlignmentReader reader = new FastaAlignmentReader();
		assertNotNull(reader);
		assertNotNull(inputStream1);
	}

	@Test
	public void testRead() throws IOException {
		FastaAlignmentReader reader = new FastaAlignmentReader();
		FastaAlignment align = reader.readAlignment(inputStream1);
		assertNotNull(align);
		assertEquals(4, align.getNumOfSequences());
		// automatic sorting by name!
		assertEquals("emelet", align.get(1).id);
		assertEquals("mellett", align.get(3).id);
	}

	@Test
	public void testDetails() throws IOException {
		FastaAlignmentReader reader = new FastaAlignmentReader();
		FastaAlignment align = reader.readAlignment(inputStream2);
		assertNotNull(align);
		GappedSequence first = align.get(0); // elme as --ELM--E-
		GappedSequence last = align.get(3); // mellett as -MELL-ETT
		assertEquals(9, first.idxs.length);
		assertEquals(9, last.idxs.length);
		assertEquals(GAP, last.idxs[0]);
		assertEquals(0, last.idxs[1]);
		assertEquals(1, last.idxs[2]);
		assertEquals(GAP, first.idxs[0]);
		assertEquals(GAP, first.idxs[1]);
		assertEquals(1, first.idxs[3]);
		assertEquals(3, first.idxs[7]);
	}

	
}
