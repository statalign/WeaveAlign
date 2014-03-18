package wvalign.io;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;



/**
 * 
 * Ancestor for sequence/alignment file reader classes that handle one specific format.
 * 
 * @author novak
 *
 */
public abstract class FileFormatReader {

	/**
	 * Reads a sequence/alignment file and constructs a RawSequences representation.
	 * 
	 * @param file File to read
	 * @return RawSequences object containing the data read
	 * @throws IOException when an I/O error occurs
	 */
	public RawSequences read(File file) throws IOException {
		return read(new BufferedReader(new FileReader(file)));
	}
	
	/**
	 * Equivalent to read(new File(fileName)).
	 * 
	 * @param fileName Name of file to read
	 * @return RawSequences object containing the data read
	 * @throws IOException when an I/O error occurs
	 */
	public RawSequences read(String fileName) throws IOException {
		return read(new BufferedReader(new FileReader(fileName)));
	}

	/**
	 * Reads sequence/alignment from a Reader and constructs a RawSequences representation.
	 * 
	 * @param reader Data source
	 * @return RawSequences object containing the data read
	 * @throws IOException when an I/O error occurs
	 */
	public abstract RawSequences read(BufferedReader reader) throws IOException;
	
}
