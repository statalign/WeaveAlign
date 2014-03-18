package wvalign.model;

import java.awt.Color;
import java.io.BufferedReader;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.io.InputStreamReader;

import wvalign.io.RawSequences;
import wvalign.utils.Utils;


/**
 * Implements the Kimura 3 parameter model.
 * 
 * @author miklos
 */
public class Kimura3 extends SubstitutionModel{

	public static final String type = "nucleotide";
	public static final String menuName = "Kimura 3 parameters";
	
	double[] oldparams;

	static final double span = 0.1;

	/**
	 * This constructor reads transition rates from the file data/kimura3_rate.dat,
	 * the alphabet from data/DNAalphabet.dat, and the equilibrium distribution from 
	 * data/kimura3_equilibrium.dat.
	 * 
	 * @throws IOException
	 */
	public Kimura3() throws IOException{
		try{
			attachedScoringScheme = new DNAScore();
		}
		catch(FileNotFoundException e){
		}
		catch(IOException e){
		}   	
		double alpha = 0.5;
		double beta  = 0.5;
		double gamma = 0.5;
		String alphabetFile = "data/DNAalphabet.dat";
		String rate = "data/kimura3_rate.dat";
		String equilibrium = "data/kimura3_equilibrium.dat";
		ClassLoader cl = getClass().getClassLoader();
		//System.out.println(cl+" "+cl.getResource(alphabetFile)+" "+alphabetFile);
		BufferedReader bf = new BufferedReader(new InputStreamReader(cl.getResource(alphabetFile).openStream()));
//		BufferedReader bf = new BufferedReader(new FileReader(alphabetFile));
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
		params = new double[3];
		oldparams = new double[3];
		params[0] = alpha; params[1] = beta; params[2] = gamma;
		bf = new BufferedReader(new InputStreamReader(cl.getResource(rate).openStream()));
		//        bf = new BufferedReader(new FileReader(rate));
		String[] temp;
		for(int i = 0; i < size; i++){
			temp = (bf.readLine()).split(" ");
			for(int j = 0; j < size; j++){
				v[i][j] = Double.parseDouble(temp[j]);
			}
		}
		String t = bf.readLine(); //reading an empty line
		for(int i = 0; i < size; i++){
			temp = (bf.readLine()).split(" ");
			for(int j = 0; j < size; j++){
				//    System.out.println(i+" "+j+" "+temp[j]);
				w[i][j] = Double.parseDouble(temp[j]);
			}
		}
		t = bf.readLine(); //reading an empty line
		bf = new BufferedReader(new InputStreamReader(cl.getResource(equilibrium).openStream()));
		//    bf = new BufferedReader(new FileReader(equilibrium));
		for(int i = 0; i < size; i++){
			t = bf.readLine();
			e[i] = Double.parseDouble(t);
		}
		setDiagonal();
	}

	void setDiagonal(){
		d[0] = 0;
		d[1]= -2*params[0]-2*params[1];
		d[2] = -2*params[1]-2*params[2];
		d[3] = -2*params[0]-2*params[2]; 

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
					//	int k;
					//	for(k = 0; k < alphabet.length && 
					//		Character.toLowerCase(sequence.charAt(j)) != Character.toLowerCase(alphabet[k]); k++);
					//	if(k == alphabet.length){
					//		System.out.println("Kimura cannot accept the sequences beacause it contains character '"+sequence.charAt(j)+"' !");
					//		return -1.0;
					//	}
					//	else{
					//		count[k] = 1;
					//	}
					if(acceptableCharacters.indexOf(Character.toLowerCase(sequence.charAt(j))) == -1){
						//System.out.println("Kimura3 cannot accept the sequences beacause it contains character '"+sequence.charAt(j)+"'!");
						throw new RecognitionError("Kimura3 cannot accept the sequences because it contains character '"+sequence.charAt(j)+"'!\n");
						//return -1.0;
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
	 *  Prints the parameters of the model
	 */
	@Override
	public String print(){
		return "\talpha\t"+params[0]+"\tbeta\t"+params[1]+"\tgamma\t"+params[2];
	}

	/**
	 * restore the parameters to the old values when a parameter-changing
	 * proposal is not accepted.
	 */
	@Override
	public void restoreParameter() {
		for(int i = 0; i < params.length; i++){
			params[i] = oldparams[i];
		}
		setDiagonal();
	}

	/**
	 * This function implements a proposal for new parameter values.
	 * Returns with the logarithm of the Metropolis-Hastings ratio.
	 */
	@Override
	public double sampleParameter() {
		for(int i = 0; i < 3; i++){
			oldparams[i] = params[i];
		}
		int k = Utils.generator.nextInt(3);
		double delta = 0.0;
		do{
			delta = Utils.generator.nextDouble()*span - (span/2.0);
		} while(params[k] + delta <= 0.0);
		params[k] += delta;
		setDiagonal();
		return Math.log((params[k] > span/2.0 ? span : span/2.0 + params[k])/
				(oldparams[k] > span/2.0 ? span : span/2.0 + oldparams[k]));
	}

	/**
	 * Returns with the most likely character, given a Felsentein likelihood array.
	 * It can handle ambiguous characters. If the probabilities in the Felsentein array
	 * are too ambiguous, it returns with a *. 
	 */
	@Override
	public char mostLikely(double[] seq) {
		int[] code = new int[4];
		char character = '*';
		for(int i = 0; i < seq.length; i++){
			if(seq[i] > 0.5){
				code[i] = 1;
			}
		}
		switch(code[0]+code[1]+code[2]+code[3]){
		case 1:{
			if(code[0] == 1){
				character = 'A';
			}
			else if(code[1] == 1){
				character = 'C';
			}
			else if(code[2] == 1){
				character = 'G';
			}
			else{
				character = 'T';
			}
			break;
		}
		case 2:{
			if(code[0] == 1 && code[1] == 1){
				character = 'M';
			}
			else if(code[0] == 1 && code[2] == 1){
				character = 'R';				
			}
			else if(code[0] == 1 && code[3] == 1){
				character = 'W';
			}
			else if(code[1] == 1 && code[2] == 1){
				character = 'S';
			}
			else if(code[1] == 1 && code[3] == 1){
				character = 'Y';
			}
			else if(code[2] == 1 && code[3] == 1){
				character = 'K';				
			}
			break;
		}
		case 3:{
			if(code[0] == 0){
				character = 'B';
			}
			else if(code[1] == 0){
				character = 'D';
			}
			else if(code[2] == 0){
				character = 'H';
			}
			else{
				character = 'V';
			}
			break;
		}
		case 4:{
			character = 'N';
			break;
		}
		default:{
			double max = 0.0;
			for(int i = 0; i < seq.length; i++){
				if(seq[i] > max){
					character =  alphabet[i];
					max = seq[i];
					//			s[1] += owner.substitutionModel.alphabet[i];
					//found = true;
				}
			}

		}
		}
		return character;
	}

}
