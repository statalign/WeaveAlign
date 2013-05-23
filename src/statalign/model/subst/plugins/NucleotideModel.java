package statalign.model.subst.plugins;

import java.awt.Color;
import java.io.BufferedReader;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.io.InputStreamReader;

import statalign.io.RawSequences;
import statalign.model.score.plugins.*;
import statalign.model.subst.RecognitionError;
import statalign.model.subst.SubstitutionModel;
//import statalign.base.*;

/**
 * This is the abstract class for nucleotide models.
 * 
 * @author lyngsoe
 */
public abstract class NucleotideModel extends SubstitutionModel{

	/* Old model parameters */
	protected double oldparams[];
	
	public static String type = "nucleotide";

	/**
	 * This constructor reads the alphabet from data/DNAalphabet.dat
	 * @throws IOException
	 */
	public NucleotideModel() throws IOException{
		try{
			attachedScoringScheme = new DNAScore();
		}
		catch(FileNotFoundException e){
		}
		catch(IOException e){
		}   	
		String alphabetFile = "data/DNAalphabet.dat";
		ClassLoader cl = getClass().getClassLoader();
		BufferedReader bf = new BufferedReader(new InputStreamReader(cl.getResource(alphabetFile).openStream()));
		String a = bf.readLine();
		alphabet = new char[a.length()];
		for(int i = 0; i < alphabet.length; i++){
			alphabet[i] = a.charAt(i);
		}
		int size = alphabet.length;
		v = new double[size][size];
		w = new double[size][size];
		d = new double[size];
		e = new double[size];
	}

	/**
	 * This function decides if this model can analyze input sequences.
	 * accepted characters are 'acgutryswmkbdhvn'. Our program is
	 * case insensitive.
	 */
	@Override
	public double acceptable(RawSequences r) {
		String acceptableCharacters = "acgutryswmkbdhvn";

		int[] count = new int[alphabet.length];
		for(int i = 0; i < r.sequences.size(); i++){
			String sequence = r.sequences.get(i);
			for(int j = 0; j < sequence.length(); j++){
				if(sequence.charAt(j) != '-' && sequence.charAt(j) != ' '){
					if(acceptableCharacters.indexOf(Character.toLowerCase(sequence.charAt(j))) == -1){
						throw new RecognitionError(getMenuName()+" cannot accept the sequences because it contains character '"+sequence.charAt(j)+"'!\n");
					}
					else if(acceptableCharacters.indexOf(Character.toLowerCase(sequence.charAt(j))) < 4){
						count[acceptableCharacters.indexOf(Character.toLowerCase(sequence.charAt(j)))] = 1;
					}

				}
			}
		}
		int sum = 0;
		for(int i = 0; i < count.length; i++){
			sum += count[i];
		}
		return (double)sum/(double)count.length;
	}

	/**
	 * This function assigns colors to characters. 'A' is red, 'C' is blue,
	 * 'G' is orange, 'T' or 'U' is green, the remaining characters (ambiguous characters
	 * and the gap symbol) are grey 
	 */
	@Override
	public Color getColor(char c) {
		if(c == 'A' || c == 'a'){
			return Color.RED;
		}
		if(c == 'C' || c == 'c'){
			return Color.BLUE;
		}
		if(c == 'G' || c == 'g'){
			return Color.ORANGE;
		}
		if(c == 'T' || c == 't' || c == 'U' || c == 'u'){
			return Color.GREEN;
		}
		return Color.GRAY;
	}

	/**
	 * returns the description of the current parameters of the model
	 */
	@Override
	public String print() {
		return "";
	}

	void SetDiagonal(){}

	/**
	 * restore the parameters to the old values when a parameter-changing
	 * proposal is not accepted.
	 */
	@Override
	public void restoreParameter() {
		for(int i = 0; i < params.length; i++){
			params[i] = oldparams[i];
		}
		SetDiagonal();
	}

	/**
	 * Returns with the most likely character, given a Felsentein
	 * likelihood array.  It can handle ambiguous characters. If the
	 * probabilities in the Felsentein array are all zero, it
	 * returns with a *.
	 */
	@Override
	public char mostLikely(double[] seq) {
		/* Create binary encoding of characters with probability
		 * higher than 0.5.
		 */
		int code = 0;
		int radix = 1;
		double max = 0.0;
		char character = '*';
		for(int i = 0; i < seq.length; i++){
			if(seq[i] > 0.5){
				/* Character has probability higher than 0.5 */
				if (max > 0.5)
					/* ...but is not the first high probability
					 * character seen.
					 */
					code += radix;
				else{
					/* ...and is the first high probability character seen */
					max = seq[i];
					code = radix;
				}
			}
			else{
				/* Character has probability at most 0.5 */
				if (seq[i] > max){
					/* ...but is most probable character seen so far */
					max = seq[i];
					code = radix;
				}
				else if (seq[i] == max)
					/* ...but is among the most probable characters
					 * seen so far.
					 */
					code += radix;
			}
			radix *= 2;
		}
		if (max == 0.0)
			code = 0;
		/* Look up corresponding character */
		String conversionTable = "*ACMGRSVTWYHKDBN";
		character = conversionTable.charAt(code);
		return character;
	}
}

