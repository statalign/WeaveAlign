package mpd.model;

import java.awt.Color;

import mpd.io.RawSequences;


/**
 * Common superclass for amino acid substitution models.
 * 
 * @author novak
 */
public abstract class AminoAcidModel extends SubstitutionModel {
	
	public static final String type = "protein";
	
	static final String magentaCharacters = "RHK";
	static final String     redCharacters = "AVFPMILW";
	static final String    blueCharacters = "DE";
	static final String   greenCharacters = "STYCNGQ";
	
	static final Color char2color[] = new Color[255];
	
	static {
		int i;
		for(i = 0; i < 255; i++)
			char2color[i] = Color.LIGHT_GRAY;
		for(i = 0; i < magentaCharacters.length(); i++)
			char2color[magentaCharacters.charAt(i)] = Color.MAGENTA;
		for(i = 0; i < redCharacters.length(); i++)
			char2color[redCharacters.charAt(i)] = Color.RED;
		for(i = 0; i < blueCharacters.length(); i++)
			char2color[blueCharacters.charAt(i)] = Color.BLUE;
		for(i = 0; i < greenCharacters.length(); i++)
			char2color[greenCharacters.charAt(i)] = Color.GREEN;
	}

	/**
	 * Returns the color of a character.
	 * 
	 * Magenta characters are 'RHK', red characters are 'AVFPMILW',
	 * blue characters are 'DE' and green characters are 'STYCNGQ'.
	 */
	@Override
	public Color getColor(char c) {
		return char2color[c];		// 'cause speed matters...
	}

	/**
	 * Dummy function with the possibility of further developments.
	 */
	@Override
	public String print() {
		return "";
	}

	/**
	 * Empty method provided for models with no parameters (the evolutionary time
	 * is represented in the edge lengths).
	 */
	@Override
	public void restoreParameter() {
	}

	/**
	 * It does nothing and returns 0, i.e log-probability 1.
	 */
	@Override
	public double sampleParameter() {
		return 0.0;
	}

	/**
	 * This function decides if the model can accept the input sequences.
	 */
	@Override
	public double acceptable(RawSequences r) {
		int[] count = new int[alphabet.length];
		for(int i = 0; i < r.sequences.size(); i++){
			String sequence = r.sequences.get(i);
			for(int j = 0; j < sequence.length(); j++){
				if(sequence.charAt(j) != '-' && sequence.charAt(j) != ' '){
					int k;
					for(k = 0; k < alphabet.length && 
						Character.toLowerCase(sequence.charAt(j)) != Character.toLowerCase(alphabet[k]); k++);
					if(k == alphabet.length){
						throw new RecognitionError(getMenuName()+" cannot accept the sequences because it contains character '"+sequence.charAt(j)+"'!\n");
						//return -1.0;
					}
					else{
						count[k] = 1;
					}
				}
			}
		}
		int sum = 0;
		for(int i = 0; i < count.length; i++){
			sum += count[i];
		}
		return (double)sum/(double)count.length + 0.1;
	}

	/**
	 * Returns the most likely character given the Felsentein likelihood array. 
	 */
	@Override
	public char mostLikely(double[] seq) {
		double max = 0.0;
		char character = '*';
		for(int i = 0; i < seq.length; i++){
			if(seq[i] > max){
				character =  alphabet[i];
				max = seq[i];
			}
		}
		return character;
	}
}
