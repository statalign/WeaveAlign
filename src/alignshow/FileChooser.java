/**
 *
 */
package alignshow;

import java.io.File;

import javax.swing.JFileChooser;
import javax.swing.filechooser.FileFilter;

import statalign.base.Utils;

/**
 * File chooser with an easy-to-use interface to set/add file filters.
 * 
 * @author miklos, novak
 */
class FileChooser extends JFileChooser { // implements ActionListener {

	private static final long serialVersionUID = 1L;

	/**
	 * Constucts a {@link FileChooser} without file filters and jumps to the current directory.
	 */
	public FileChooser() {
//		setCurrentDirectory(new File(System.getProperty("user.dir")));
		setCurrentDirectory(new File("."));		// simpler
	}
	
	/**
	 * Constructs a {@link FileChooser} with a given file filter and jumps to the current directory.
	 * <p>See {@link #addFilter(String[], String)} for description of params.
	 */
	public FileChooser(String[] extensions, String description) {
		this();
		setFileFilter(new MyFilter(extensions, description));
	}

	/**
	 * Adds a new filter to the file filter list.
	 * 
	 * @param extensions list of file extensions that are accepted when filter is selected
	 * @param description short description of the file type
	 */
	public void addFilter(String[] extensions, String description) {
		setFileFilter(new MyFilter(extensions, description));
	}

	private static class MyFilter extends FileFilter {

		String[] extensions;
		String description;

		MyFilter(String[] extensions, String description) {
			this.extensions = extensions;
			this.description = description;
		}

		@Override
		public boolean accept(File f) {
			for (String ext : extensions) {
				if (f.isDirectory()
						|| f.getName().toUpperCase()
								.endsWith(ext.toUpperCase()))
					return true;
			}
			return false;
		}

		@Override
		public String getDescription() {
			return description + " ("
					+ Utils.joinStrings(extensions, "*.", ", ") + ")";
		}

	}

	/**
	 * Adds the default (first) extension of the selected filter to the current file name
	 * unless it has an extension recognised by the filter.
	 * 
	 * @return the current (selected) file, after possible concatenation of the extension
	 */
	public File addFilterExtension() {
		File file = getSelectedFile();
		String name = file.getAbsolutePath();
		FileFilter filter = getFileFilter();
		if(!(filter instanceof MyFilter))
			return file;
		String[] exts = ((MyFilter)filter).extensions;
		for(String ext : exts)
			if(name.endsWith("."+ext))
				return file;
		file = new File(name+"."+exts[0]);
		setSelectedFile(file);
		return file;
	}

}
