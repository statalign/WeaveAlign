package wvalign.utils;

import static org.junit.Assert.*;

import org.junit.Test;

public class CircularArrayTest {

	@Test
	public void testNotNull() {
		CircularArray<Integer> array = new CircularArray<Integer>();
		assertNotNull(array);
	}

	@Test
	public void testPushValues() {
		CircularArray<Integer> array = new CircularArray<Integer>();
		Integer first = 2;
		Integer second = 4;
		Integer third = 6;
		
		assertEquals(0, array.length());
		
		array.push(first);
		array.push(second);
		array.push(third);
		
		assertEquals(3, array.length());
		
		assertEquals(third, array.pop());
		assertEquals(second, array.pop());
		assertEquals(first, array.pop());
		
		assertEquals(0, array.length());
	}

	@Test
	public void testPutValues() {
		CircularArray<Integer> array = new CircularArray<Integer>();
		Integer first = 2;
		Integer second = 4;
		Integer third = 6;
		
		assertEquals(0, array.length());
		
		array.put(2, first);
		array.put(4, second);
		array.put(8, third);
		
		// num of elements to store: max_id - min_id + 1 = 8 - 2 + 1
		assertEquals(7, array.length());
		
		assertNull(array.get(23));
		assertEquals(second, array.get(4));
		assertEquals(third, array.get(8));		
	}

}
