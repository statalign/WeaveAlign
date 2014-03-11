package mpd;

import java.util.Arrays;


class Column {
	ColumnKey key;
	ColClass succ;

	ColClass pred;		// not maintained by default, used e.g. for annotation
	double[] scores;
	
	int count = 1;

	Column(ColumnKey _key) {
		key = _key;
	}

	@Override
	public String toString() {
		return Arrays.toString(key.desc);
	}
	

	public static class ColumnKey {
		public int[] desc;

		ColumnKey(int[] arr) {
			desc = arr;
		}

		@Override
		public boolean equals(Object o)	{
			return (o instanceof ColumnKey) && Arrays.equals(desc, ((ColumnKey)o).desc);
		}

		@Override
		public int hashCode()	{
			return Arrays.hashCode(desc);
		}

//			/**
//			 * Calculates key of predecessor column class for this column.
//			 * (Column x can directly precede column y iff x.succ()=y.pred() ) 
//			 */
//			ColumnKey pred() {
//				int[] desc = this.desc, ret = new int[desc.length];
//				for(int i = 0; i < desc.length; i++)
//					ret[i] = desc[i] >> 1;
//				return new ColumnKey(ret);
//			}

		/**
		 * Calculates key of successor column class for this column.
		 */
		ColumnKey succ() {
			int[] desc = this.desc, ret = new int[desc.length];
			for(int i = 0; i < desc.length; i++)
				ret[i] = (desc[i]+1) >> 1;
			return new ColumnKey(ret);
		}
		
		/**
		 * Returns gap insensitive column key variant of the key.
		 */
		ColumnKey giKey() {
			int len = desc.length, d;
			int[] desc = this.desc, sdesc = new int[len]; 
			for(int i = 0; i < len; i++) {
				d = desc[i];
				sdesc[i] = ((d&1)==1)?d:0;
			}
			return new ColumnKey(sdesc);
		}
		


		static int colNext(int n) {
			return n + (n & 1);
		}
		
		public void print() {
			for (int i : desc) {
				System.out.print(i+" ");
			}
			System.out.println();
		}

	}
}