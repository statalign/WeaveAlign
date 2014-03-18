package mpd.model;

import java.io.*;

import mpd.io.FileTokenReader;

/**
 * 
 * An implementation of the BLOSUM62 distances
 * This is the distance and not the similarity matrix!!!
 * 
 * @author miklos, novak
 *
 */
public class Blosum62 extends SubstitutionScore{

	final String aaCodes = "ARNDCQEGHILKMFPSTWYV";   /* no B,Z or U currently */

	/**
	 * This constructor reads the distances from data/aadist.dat
	 * and sets the which[][] array.
	 * @throws IOException
	 * @throws FileNotFoundException
	 */
	public Blosum62() throws IOException, FileNotFoundException{
		which = new int[256][20];   /* tells amino acid no. for each valid char */
		dist = new int[20][20];  /* amino acid distances in Gotoh's alg. */

		//	System.out.println("I'm in Blosum62 constructor!");
		//	System.out.println("I'm in Blosum62 constructor! "+which);
		//System.out.println("I'm in Blosum62 constructor! "+which.length);

		FileTokenReader file = new FileTokenReader("data/aadist.dat");
		for(int i = 0; i < 20; i++){
			for(int j = 0; j <= i; j++){
				dist[i][j] = dist[j][i] = file.readDbl().intValue();
			}
		}
		file.close();

//		for(int i = 0; i < 256; i++){
	//		which[i] = -1;  /* no amino acid associated */
		for(int i = 0; i < aaCodes.length(); i++){
			which[aaCodes.charAt(i)][i] = which[Character.toLowerCase(aaCodes.charAt(i))][i] = 1;
		}
		which[128] = which[135] = which['c'];   /* cystein with S-S bond */

	}


}