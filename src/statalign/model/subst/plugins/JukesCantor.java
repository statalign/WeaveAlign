package statalign.model.subst.plugins;

import java.awt.Color;
import java.io.BufferedReader;
import java.io.FileNotFoundException;
//import java.io.FileReader;
import java.io.IOException;
import java.io.InputStreamReader;

//import statalign.base.Utils;
import statalign.io.RawSequences;
import statalign.model.score.plugins.DNAScore;
import statalign.model.subst.RecognitionError;
import statalign.model.subst.SubstitutionModel;

/**
 * Implements the Jukes-Cantor model for nucleic acids.
 *
 * @author miklos, novak
 *
 */
public class JukesCantor extends SubstitutionModel{

	public static final String menuName = "Jukes-Cantor";
	public static final String type = "nucleotide";

	static final double span = 0.1;

	/**
	 * This constructor reads transition rates from the file data/jukescantor_rate.dat,
	 * the alphabet from data/DNAalphabet.dat, and the equilibrium distribution from 
	 * data/jukescantor_equilibrium.dat.
	 * 
	 * @throws IOException
	 */
	public JukesCantor() throws IOException{
		try{
			attachedScoringScheme = new DNAScore();
		}
		catch(FileNotFoundException e){
		}
		catch(IOException e){
		}
		String alphabetFile = "data/DNAalphabet.dat";
		String rate = "data/jukescantor_rate.dat";
		String equilibrium = "data/jukescantor_equilibrium.dat";
		ClassLoader cl = getClass().getClassLoader();
		//System.out.println(cl+" "+cl.getResource(alphabetFile)+" "+alphabetFile);
		BufferedReader bf = new BufferedReader(new InputStreamReader(cl.getResource(alphabetFile).openStream()));
		//BufferedReader bf = new BufferedReader(new FileReader(alphabetFile));
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
		params = new double[0];
		bf = new BufferedReader(new InputStreamReader(cl.getResource(rate).openStream()));
		//      bf = new BufferedReader(new FileReader(rate));
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
		//     bf = new BufferedReader(new FileReader(equilibrium));
		for(int i = 0; i < size; i++){
			t = bf.readLine();
			e[i] = Double.parseDouble(t);
		}
		d[0] = 0.0;
		d[1] = -4.0/3.0;
		d[2] = d[1];
		d[3] = d[1];
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
					//		System.out.println("Jukes-Cantor cannot accept the sequences beacause it contains character '"+sequence.charAt(j)+"' !");
					//		return -1.0;
					//	}
					//	else{
					//		count[k] = 1;
					//	}
					if(acceptableCharacters.indexOf(Character.toLowerCase(sequence.charAt(j))) == -1){
						//System.out.println("Jukes-Cantor cannot accept the sequences beacause it contains character '"+sequence.charAt(j)+"'!");
						throw new RecognitionError("Jukes-Cantor cannot accept the sequences because it contains character '"+sequence.charAt(j)+"'!\n");
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
	 * dummy function for the possibility of further developments.
	 */
	@Override
	public String print() {
		return "";
	}

	/**
	 * Empty function, since the only parameter of the model is handled by
	 * the edge lengths
	 */
	@Override
	public void restoreParameter() {
	}

	/**
	 * It does nothing, and always return with 0, namely, log-probability 1.
	 */
	@Override
	public double sampleParameter() {
		return 0.0;
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

	/**
	 * For testing purposes only.
	 * @param args No argument is used.
	 */
	public static void main(String[] args){
		try {
			JukesCantor jc = new JukesCantor();
			double[][] p = jc.updateTransitionMatrix(null, 1.0);
			for(int i = 0; i < p.length; i++){
				for(int j = 0; j < p[i].length; j++){
					System.out.print(p[i][j]+" ");
				}
				System.out.println();
			}
		} catch (IOException e) {
			e.printStackTrace();
		}

	}


}
