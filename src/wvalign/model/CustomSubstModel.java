package wvalign.model;

import java.awt.Color;
import java.util.Arrays;

import wvalign.io.RawSequences;



public class CustomSubstModel extends SubstitutionModel {
	
	public static String type = "custom";
	
	private double[][] rates;
	
	public CustomSubstModel(final String alphabet, double[][] rates, double[] eq) {
		attachedScoringScheme = new SubstitutionScore() {
			{
				which = new int[256][];
				for(int i = 0; i < alphabet.length(); i++) {
					int chl = Character.toLowerCase(alphabet.charAt(i));
					int chu = Character.toUpperCase(chl);
					which[chl] = new int[alphabet.length()];
					which[chu] = new int[alphabet.length()];
					which[chl][i] = which[chu][i] = 1;
				}
			}
		};
		this.alphabet = alphabet.toCharArray();
		this.rates = new double[rates.length][];
		for(int i = 0; i < rates.length; i++)
			this.rates[i] = Arrays.copyOf(rates[i], rates[i].length);
		e = eq;
		d = new double[eq.length];
		v = new double[eq.length][eq.length];
		w = new double[eq.length][eq.length];
		NucleotideModelUtils.numericalDiagonalisation(this, rates);
	}
	
	private static final char GAPCHAR = '-';
	
	public void addGapChar(double insRate, double delRate) {
		int len = e.length+1;
		
		alphabet = (new String(alphabet)+GAPCHAR).toCharArray();
		int[][] wh = attachedScoringScheme.which;
		for(int i = 0; i < wh.length; i++)
			if(wh[i] != null)
				wh[i] = Arrays.copyOf(wh[i], len);
		wh[GAPCHAR] = new int[len];
		wh[GAPCHAR][len-1] = 1;
		
		insRate *= e.length;
		rates = Arrays.copyOf(rates, len);
		for(int i = 0; i < len-1; i++) {
			rates[i] = Arrays.copyOf(rates[i], len);
			rates[i][len-1] = delRate;
			rates[i][i] = 0;
			rates[i][i] = -sum(rates[i]);
		}
		rates[len-1] = new double[len];
		for(int i = 0; i < len-1; i++)
			rates[len-1][i] = insRate*e[i];
		rates[len-1][len-1] = -insRate;

		double eg = delRate/(insRate+delRate);		// equilibrium prob of gap char
		e = Arrays.copyOf(e, len);
		for(int i = 0; i < len-1; i++)
			e[i] -= e[i]*eg;
		e[len-1] = eg;
//		System.out.println(Arrays.toString(e));

		d = new double[len];
		v = new double[len][len];
		w = new double[len][len];
		NucleotideModelUtils.numericalDiagonalisation(this, rates);
	}
		
	private double sum(double[] arr) {
		double s = 0;
		for(double v : arr)
			s += v;
		return s;
	}
	
	@Override
	public double acceptable(RawSequences r) {
		// TODO Auto-generated method stub
		return 0;
	}

	@Override
	public double sampleParameter() {
		return 0;
	}

	@Override
	public void restoreParameter() {
	}

	@Override
	public String print() {
		return null;
	}

	@Override
	public Color getColor(char c) {
		return null;
	}

	@Override
	public char mostLikely(double[] seq) {
		return 0;
	}
}
