package wvalign.utils;

import java.util.Enumeration;
import java.util.Iterator;
import java.util.Random;

/**
 * 
 * This class contains multi-purpose static functions.
 * 
 * @author miklos, novak
 *
 */
public class Utils{

	private Utils(){
		//... so that no objects of this class can be created
	}
	
	/**
	 * The random number generator used throughout the program.
	 * A new generator is constructed at each MCMC run using the seed in the
	 * corresponding MCMCPars object.
	 */
	public static Random generator = new Random(1);
	
	/**
	 * log(0) is set to Double.NEGATIVE_INFINITY. This is used in logarithmic adding.
	 * The logarithm of an empty sum is set to this value.
	 * 
	 */
	public static final double log0 = Double.NEGATIVE_INFINITY;

	/**
	 * Logarithmically add two numbers
	 * @param a log(x)
	 * @param b log(y)
	 * @return log(x+y)
	 */
	public static double logAdd(double a, double b) {
		if(a == b)
			return Math.log(2)+a;
		if(a < b)
			return b+Math.log(Math.exp(a-b)+1);
		return a+Math.log(Math.exp(b-a)+1);
	}

    static public void calcFelsen(double[] res, double[] fel1, double[][] prob1, double[] fel2, double[][] prob2) {
		double s;
		int i, j, len = res.length;

		for(i = 0; i < len; i++) {
			if(fel1 != null) {
				s = 0.0;
				for(j = 0; j < len; j++)
					s += prob1[i][j]*fel1[j];
				res[i] = s;
			} 
			else{
				res[i] = 1.0;
			}
			if(fel2 != null) {
				s = 0.0;
				for(j = 0; j < len; j++){
					s += prob2[i][j]*fel2[j];
				}
				res[i] *= s;
			}
		}
	}

	/**
	 *  Calculates emission probability from Felsenstein likelihoods
	 */
   static public double calcEmProb(double fel[], double aaEquDist[]) {
		double p = 0;

		for(int i = 0; i < fel.length; i++)
			p += fel[i]*aaEquDist[i];
		return p;
	}
	
	/**
	 * Makes Enumeration iterable.
	 * 
	 * @param <T> Enumeration element type
	 * @param en  the Enumeration
	 * @return an Iterable that can iterate through the elements of the Enumeration
	 */
	public static <T> Iterable<T> iterate(final Enumeration<T> en) {
		final Iterator<T> iterator = new Iterator<T>() {
			@Override
			public boolean hasNext() {
				return en.hasMoreElements();
			}
			@Override
			public T next() {
				return en.nextElement();  
			}
			@Override
			public void remove() {
				throw new UnsupportedOperationException();
			}
		};
		return new Iterable<T>() {
			@Override
			public Iterator<T> iterator() {
				return iterator;
			}
		};
	}
	
	/**
	 * Joins strings using a separator string. Accepts any <code>Object</code>s
	 * converting them to strings using their <code>toString</code> method.
	 * 
	 * @param strs
	 *            strings to join
	 * @param separator
	 *            the separator string
	 * @return a string made up of the strings separated by the separator
	 */
	public static String joinStrings(Object[] strs, String separator) {
		StringBuilder result = new StringBuilder();
		int i, l = strs.length;
		if (l > 0) {
			result.append(strs[0]);
		}
		for (i = 1; i < l; i++) {
			result.append(separator);
			result.append(strs[i]);
		}
		return result.toString();
	}

	/**
	 * Joins strings using a prefix and a separator string. Accepts any
	 * <code>Object</code>s converting them to strings using their
	 * <code>toString</code> method.
	 * 
	 * @param strs
	 *            strings to join
	 * @param prefix
	 *            prefix for each string
	 * @param separator
	 *            the separator string
	 * @return a string made up of the strings with the given prefix and
	 *         separated by the separator
	 */
	public static String joinStrings(Object[] strs, String prefix, String separator) {
		StringBuilder result = new StringBuilder();
		int i, l = strs.length;
		if (l > 0) {
			result.append(prefix);
			result.append(strs[0]);
		}
		for (i = 1; i < l; i++) {
			result.append(separator);
			result.append(prefix);
			result.append(strs[i]);
		}
		return result.toString();
	}

	public static double[][] toLog(double[][] m) {
		double[][] ret = new double[m.length][];
		for(int i = 0; i < m.length; i++)
			ret[i] = toLog(m[i]);
		return ret;
	}
	
	public static double[] toLog(double[] v) {
		double[] ret = new double[v.length];
		for(int i = 0; i < v.length; i++)
			ret[i] = Math.log(v[i]);
		return ret;
	}

}

