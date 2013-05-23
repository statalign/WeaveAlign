package alignshow;

import java.awt.BorderLayout;
import java.awt.Color;
import java.awt.Cursor;
import java.awt.Dimension;
import java.awt.FlowLayout;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collections;
import java.util.Comparator;
import java.util.HashMap;
import java.util.List;

import javax.swing.JButton;
import javax.swing.JCheckBox;
import javax.swing.JFileChooser;
import javax.swing.JFrame;
import javax.swing.JOptionPane;
import javax.swing.JPanel;
import javax.swing.JScrollPane;
import javax.swing.SwingUtilities;
import javax.swing.UIManager;

import ml.options.OptionData;
import ml.options.OptionSet;
import ml.options.Options;
import ml.options.Options.Multiplicity;
import ml.options.Options.Separator;
import statalign.io.RawSequences;
import statalign.model.subst.RecognitionError;
import statalign.model.subst.SubstitutionModel;
import statalign.model.subst.plugins.Dayhoff;
import statalign.model.subst.plugins.Kimura3;

public class Show extends JFrame {
	
	private static final long serialVersionUID = 1L;

	public static final String MPDSHOW_VERSION = "v1.1.3";

	private static final String USAGE =
		"MPDSHOW "+MPDSHOW_VERSION+" (C) Adam Novak, 2010-11.\n\n"
		+
		"Usage:\n\n" +
		"    java -jar alignshow.jar [options] file.fsa/mpd\n\n\n"
		+
		"Description:\n\n" +
		"    Shows alignment and annotation tracks (plots above alignment).\n\n\n"
		+
		"Options:\n\n" +
		"    -t FILE\n" +
		"        Adds an annotation track (FILE should contain a number for each column)\n\n"
		+
		"    -c=COLOR\n" +
		"        Sets the plotting color of a track (when used more than once, it\n" +
		"          affects the color of the tracks in the order they are added\n" +
		"          with -t, the MPD posterior track being first when present)\n" +
		"        E.g. -c=BLUE -c=yellow -c=0,255,255\n\n"
		+
		"    -o ORDERFILE\n" +
		"        Specifies the order of the sequences.\n\n"
		+
		"    -g GROUPINGFILE\n" +
		"        Sequence grouping markup per column.\n\n";

	private static Color[] DEF_COLORS = {
		Color.blue, Color.red, Color.green, Color.black, Color.orange, Color.gray
	};
	
	private AlignmentGUI alignGui;
	
	public Show(AlignmentGUI align) {
		super("Annotated alignment");
		alignGui = align;
		
		// set system look & feel
		try {
			String syslook = UIManager.getSystemLookAndFeelClassName();
			if (!UIManager.getLookAndFeel().getClass().getName()
					.equals(syslook)) {
				UIManager.setLookAndFeel(syslook);
				SwingUtilities.updateComponentTreeUI(this);
			}
		} catch (Exception ex) {
		}
		
		JScrollPane scroll = new JScrollPane(align);
		setLayout(new BorderLayout());
		
		getContentPane().add(scroll, BorderLayout.CENTER);
		
		JPanel panel = new JPanel();
		panel.setLayout(new FlowLayout());

		final JCheckBox box = new JCheckBox("Show filename");
		box.setSelected(true);
		box.addActionListener(new ActionListener() {
			public void actionPerformed(ActionEvent e) {
				alignGui.setShowTitle(box.isSelected());
			}
		});
		panel.add(box);
		panel.add(new JPanel());
		
		JButton button = new JButton("Save as PNG image...");
		button.addActionListener(new ActionListener() {
			public void actionPerformed(ActionEvent e) {
				FileChooser choose = new FileChooser(new String[] { "png" }, "PNG files");
				if(choose.showSaveDialog(Show.this) == JFileChooser.APPROVE_OPTION) {
					File file = choose.addFilterExtension();
					if(!checkOverwrite(file))
						return;
					try {
						setCursor(Cursor.getPredefinedCursor(Cursor.WAIT_CURSOR));
						alignGui.saveAsImage(file);
						setCursor(null);
					} catch (IOException ex) {
						JOptionPane.showMessageDialog(Show.this, "An error occured while saving the image.", "Error", JOptionPane.ERROR_MESSAGE);
					}
				}
			}
		});
		panel.add(button);
		panel.add(new JPanel());
		
		button = new JButton("Save as EPS file...");
		button.addActionListener(new ActionListener() {
			public void actionPerformed(ActionEvent e) {
				FileChooser choose = new FileChooser(new String[] { "eps", "ps" }, "EPS files");
				if(choose.showSaveDialog(Show.this) == JFileChooser.APPROVE_OPTION) {
					File file = choose.addFilterExtension();
					if(!checkOverwrite(file))
						return;
					try {
						setCursor(Cursor.getPredefinedCursor(Cursor.WAIT_CURSOR));
						alignGui.saveAsEps(file);
						setCursor(null);
					} catch (IOException ex) {
						JOptionPane.showMessageDialog(Show.this, "An error occured while saving the image.", "Error", JOptionPane.ERROR_MESSAGE);
					}
				}
			}
		});
		panel.add(button);

		getContentPane().add(panel, BorderLayout.SOUTH);

		setDefaultCloseOperation(EXIT_ON_CLOSE);
		pack();
		setLocationRelativeTo(null);
		if(scroll.getHorizontalScrollBar().isVisible()) {// && getContentPane().getPreferredSize().height < Toolkit.getDefaultToolkit().getScreenSize().height) {
			getContentPane().setPreferredSize(addDimY(addDimY(alignGui.getPreferredSize(), panel.getPreferredSize()),
					scroll.getHorizontalScrollBar().getSize()));
			pack();
		}
		setVisible(true);
	}
	
	private Dimension addDimY(Dimension d1, Dimension d2) {
		return new Dimension(Math.max(d1.width,d2.width), d1.height+d2.height+5);
	}
	
	private boolean checkOverwrite(File file) {
		if(file.exists() && JOptionPane.showConfirmDialog(Show.this, "Overwrite existing file?", "Confirm overwrite", JOptionPane.OK_CANCEL_OPTION) != JOptionPane.OK_OPTION)
			return false;
		return true;
	}
	
	public static void main(String[] args) {
		Options opt = new Options(args, Multiplicity.ZERO_OR_ONE, 1, 1);
		opt.addSet("run")
				.addOption("t", Separator.BLANK, Multiplicity.ZERO_OR_MORE)
				.addOption("c", Separator.EQUALS, Multiplicity.ZERO_OR_MORE)
				.addOption("o", Separator.BLANK, Multiplicity.ZERO_OR_ONE)
				.addOption("g", Separator.BLANK, Multiplicity.ZERO_OR_ONE);
		
		OptionSet set;
		if ((set = opt.getMatchingSet(false, false)) == null) {
			System.out.println(USAGE);
			System.exit(1);
		}
		
		try {
			MpdReader reader = new MpdReader();
			File input = new File(set.getData().get(0));
			RawSequences seqs = reader.read(input);
			if(reader.getErrors() > 0) {
				error("error reading Fasta/MPD file "+input);
			}
			if(seqs.sequences.size() == 0 || seqs.sequences.get(0).length() == 0) {
				error("no sequences read or sequences empty");
			}
			
			int alignSize = seqs.sequences.get(0).length();

			SubstitutionModel model = new Kimura3();
			try {
				model.acceptable(seqs);
			} catch (RecognitionError r) {
				model = new Dayhoff();
			}

			AlignmentGUI alignGui = new AlignmentGUI(input.getPath(), model);
	
			HashMap<String, Integer> orderMap = null;
			if(set.isSet("o")) {
				String ofile = set.getOption("o").getResultValue(0);
				BufferedReader oreader = new BufferedReader(new FileReader(ofile));
				
				orderMap = new HashMap<String, Integer>();
				String line;
				int lno = 0;
				while((line = oreader.readLine()) != null) {
					orderMap.put(line.trim(), lno++);
				}
			}
			
			alignGui.alignment = convertAlign(seqs, orderMap);
			
//			System.out.println("as");
			if(set.isSet("g")) {
				int n = alignGui.alignment.length;
				int[] convTab = new int[n];
				if(orderMap != null) {
					List<String> names = new ArrayList<String>();
					names.addAll(orderMap.keySet());
					Collections.sort(names);
					for(int i = 0; i < n; i++)
						convTab[i] = orderMap.get(names.get(i));
				} else {
					for(int i = 0; i < n; i++)
						convTab[i] = i;
				}
				
				
				String ofile = set.getOption("g").getResultValue(0);
				BufferedReader greader = new BufferedReader(new FileReader(ofile));
				
				List<List<Integer>> groups = new ArrayList<List<Integer>>();
				String line;
				while((line = greader.readLine()) != null) {
					List<Integer> list = new ArrayList<Integer>();
					String[] strs = line.split(" ");
					for(String str : strs)
						list.add(convTab[Integer.parseInt(str)]);
					Collections.sort(list);
					groups.add(list);
				}
				
				alignGui.grouping = groups;
			}
			
			List<Double> mpdScoreTrack = reader.getScores();
			if(mpdScoreTrack != null) {
				if(mpdScoreTrack.size() != alignSize)
					System.err.println("warning: MPD annotation track is incompatible with alignment!");
				alignGui.tracks.add(mpdScoreTrack);
				alignGui.colors.add(DEF_COLORS[0]);
			}
			
			OptionData trackData = set.getOption("t");
			TrackReader treader = new TrackReader();
			for(int i = 0; i < trackData.getResultCount(); i++) {
				List<Double> track = treader.read(trackData.getResultValue(i));
				if(track.size() != alignSize)
					System.err.println("warning: track no. "+(i+1)+" is incompatible with alignment!");
				alignGui.tracks.add(track);
				alignGui.colors.add(DEF_COLORS[(alignGui.colors.size())%DEF_COLORS.length]);
			}
			
			OptionData colorData = set.getOption("c");
			for(int i = 0; i < colorData.getResultCount(); i++) {
				alignGui.colors.set(i, getColor(colorData.getResultValue(i)));
			}

			new Show(alignGui);
		} catch (Exception e) {
			e.printStackTrace();
		}
	}
	
	private static Color getColor(String color) {
		if(color.isEmpty())
			return null;
		try {
			if(Character.isDigit(color.charAt(0))) {
				String[] vals = color.split(",");
				return new Color(Integer.parseInt(vals[0]), Integer.parseInt(vals[1]), Integer.parseInt(vals[2]));
			}
		    return (Color)Color.class.getField(color).get(null);
		} catch (Exception e) {
		    return null;
		}
	}
	
	private static void error(String string) {
		System.err.println("mpdshow: "+string);
		System.exit(1);
	}

	private static String[] convertAlign(RawSequences seqs, final HashMap<String, Integer> orderMap) {
		int n = seqs.sequences.size();
		if(n == 0)
			return new String[0];
		
		int maxLen = seqs.seqNames.get(0).length();
		for(String name : seqs.seqNames)
			maxLen = Math.max(maxLen, name.length());
		
		List<String> conv = new ArrayList<String>();
		StringBuilder build = new StringBuilder();
		for(int i = 0; i < n; i++) {
			build.setLength(0);
			String name = seqs.seqNames.get(i);
			build.append(name);
			for(int j = name.length(); j < maxLen; j++)
				build.append(" ");
			build.append("\t");
			build.append(seqs.sequences.get(i));
			conv.add(build.toString());
		}
		if(orderMap != null)
			Collections.sort(conv, new Comparator<String>() {
				@Override
				public int compare(String o1, String o2) {
					return orderMap.get(o1.split("\t")[0].trim())-orderMap.get(o2.split("\t")[0].trim());
				}
			});
		return conv.toArray(new String[n]);
	}

}
