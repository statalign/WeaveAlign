package statalign.ui;

import java.awt.BorderLayout;
import java.awt.Container;
import java.awt.Cursor;
import java.awt.Dimension;
import java.awt.Toolkit;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.awt.event.KeyEvent;
import java.io.File;
import java.io.IOException;
import java.lang.reflect.Modifier;
import java.util.ArrayList;
import java.util.HashMap;

import javax.swing.ButtonGroup;
import javax.swing.JFileChooser;
import javax.swing.JFrame;
import javax.swing.JLabel;
import javax.swing.JMenu;
import javax.swing.JMenuBar;
import javax.swing.JMenuItem;
import javax.swing.JOptionPane;
import javax.swing.JPanel;
import javax.swing.JRadioButtonMenuItem;
import javax.swing.JTabbedPane;
import javax.swing.JWindow;
import javax.swing.KeyStroke;
import javax.swing.SwingUtilities;
import javax.swing.UIManager;

import statalign.StatAlign;
import statalign.base.Input;
import statalign.base.MainManager;
import statalign.base.Utils;
import statalign.io.input.FileFormatReader;
import statalign.io.input.plugins.FastaReader;
import statalign.model.subst.RecognitionError;
import statalign.model.subst.SubstitutionModel;
import statalign.model.subst.plugins.Dayhoff;
import statalign.model.subst.plugins.Kimura3;
import statalign.postprocess.Postprocess;

/**
 * The main frame of the program.
 * 
 * Created by <code>statalign.StatAlign</code> when run in graphical mode.
 * 
 * @author miklos, novak
 *
 */
public class MainFrame extends JFrame implements ActionListener {
	public static final String IDLE_STATUS_MESSAGE = " Phylogeny Cafe :: Statistical Alignment";
	private static final long serialVersionUID = 1L;

	Thread workerThr = null;
	
	JMenuItem openItem;
	JMenuItem ioPreferencesItem;
	JMenuItem runItem;
	JMenuItem pauseItem;
	JMenuItem resumeItem;
	JMenuItem stopItem;
	JWindow message;
	McmcSettingsDlg mcmcSettingsDlg;
	
	JMenuItem[] modelButtons;
	
	/**
	 * The main manager that handles the MCMC run.
	 */
	public MainManager manager;
	
	File inFile;
	
	JPanel statusBar;
	public JLabel statusText;

	private Class<? extends SubstitutionModel>[] substModels;
	
	/**
	 * The only constructor of the class. It launches the main window.
	 */
	@SuppressWarnings("unchecked")
	public MainFrame() throws Exception {
		super("StatAlign "+StatAlign.version);
	
		try {
			String syslook = UIManager.getSystemLookAndFeelClassName();
			if(!UIManager.getLookAndFeel().getClass().getName().equals(syslook)) {
				UIManager.setLookAndFeel(syslook);
				SwingUtilities.updateComponentTreeUI(this);
			}
		} catch(Exception ex) { JOptionPane.showMessageDialog(this, ex); }

		manager = new MainManager(this);

		ArrayList<Class<?>> substModList = new ArrayList<Class<?>>();
		for(String model : Utils.classesInPackage(SubstitutionModel.class.getPackage().getName()+".plugins")) {
			try {
				Class<?> cl = Class.forName(model);
				/* Only include non-abstract substitution models that extend SubstitutionModel */
				if (!Modifier.isAbstract(cl.getModifiers()) && 
					SubstitutionModel.class.isAssignableFrom(cl))
					substModList.add(cl);
			} catch(Exception ex) {
				ErrorMessage.showPane(null, ex, true);
			}
		}
		substModels = (Class<? extends SubstitutionModel>[])substModList.toArray(new Class<?>[substModList.size()]);
		
		mcmcSettingsDlg = new McmcSettingsDlg(this);
		JMenuBar menubar = new JMenuBar();
		
		
		JMenu menu = new JMenu("File");
		menu.setMnemonic(KeyEvent.VK_F);
		JMenuItem item;
		
		openItem = new JMenuItem("Add sequence(s)...");
		openItem.addActionListener(this);
		openItem.setAccelerator(KeyStroke.getKeyStroke("control O"));
		openItem.setMnemonic(KeyEvent.VK_A);
		menu.add(openItem);
		
		ioPreferencesItem = new JMenuItem("Preferences...");
		ioPreferencesItem.addActionListener(this);
		ioPreferencesItem.setAccelerator(KeyStroke.getKeyStroke("control 1"));
		ioPreferencesItem.setMnemonic(KeyEvent.VK_P);
		menu.add(ioPreferencesItem);
		
		menu.addSeparator();
		
		item = new JMenuItem("Exit");
		item.addActionListener(this);
		item.setAccelerator(KeyStroke.getKeyStroke("alt F4"));
		item.setMnemonic(KeyEvent.VK_X);
		menu.add(item);
		menubar.add(menu);
		
//		menu = new JMenu("Edit");
//		menu.setMnemonic(KeyEvent.VK_E);
//		item = new JMenuItem("Cut");
//		item.addActionListener(this);
//		item.setAccelerator(KeyStroke.getKeyStroke("control X"));
//		menu.add(item);
//		item = new JMenuItem("Copy");
//		item.addActionListener(this);
//		item.setAccelerator(KeyStroke.getKeyStroke("control C"));
//		menu.add(item);
//		item = new JMenuItem("Paste");
//		item.addActionListener(this);
//		item.setAccelerator(KeyStroke.getKeyStroke("control V"));
//		menu.add(item);
//		menubar.add(menu);
		
		
		menu = new JMenu("MCMC");
		menu.setMnemonic(KeyEvent.VK_M);
		
		item = new JMenuItem("Settings");
		item.addActionListener(this);
		item.setAccelerator(KeyStroke.getKeyStroke("control M"));
		item.setMnemonic(KeyEvent.VK_S);
		menu.add(item);
		
		runItem = new JMenuItem("Run");
		runItem.setAccelerator(KeyStroke.getKeyStroke("control ENTER"));
		runItem.setEnabled(false);
		runItem.addActionListener(this);
		menu.add(runItem);
		
		pauseItem = new JMenuItem("Pause");
		pauseItem.setEnabled(false);
		pauseItem.addActionListener(this);
		menu.add(pauseItem);
		
		resumeItem = new JMenuItem("Resume");
		resumeItem.setEnabled(false);
		resumeItem.addActionListener(this);
		menu.add(resumeItem);
		
		stopItem = new JMenuItem("Stop");
		stopItem.setEnabled(false);
		stopItem.addActionListener(this);
		menu.add(stopItem);
		menubar.add(menu);
		
		
		menu = new JMenu("Model");
		menu.setMnemonic(KeyEvent.VK_L);

		modelButtons = new JMenuItem[substModels.length];
		ButtonGroup modelGroup = new ButtonGroup();
		HashMap<String,ArrayList<Class<? extends SubstitutionModel>>> substModTypes = 
				new HashMap<String, ArrayList<Class<? extends SubstitutionModel>>>();
		for(Class<? extends SubstitutionModel> cl : substModels) {
			String type = SubstitutionModel.getType(cl);
			ArrayList<Class<? extends SubstitutionModel>> arr = substModTypes.get(type);
			if(arr == null)
				substModTypes.put(type, arr = new ArrayList<Class<? extends SubstitutionModel>>());
			arr.add(cl);
		}
		String[] typeArr = new String[substModTypes.keySet().size()];
		int s = 0;
		if(typeArr.length >= 2) {
			typeArr[0] = "protein";				// amino acid subst models first
			typeArr[1] = "nucleotide";		// then nucleotide models, then the rest
			s = 2;
		}
		for(String type : substModTypes.keySet()) {
			if(!type.equals("protein") && !type.equals("nucleotide"))
				typeArr[s++] = type;
		}
		s = 0;
		for(String type : typeArr) {
			if(s > 0)
				menu.addSeparator();
			for(Class<? extends SubstitutionModel> cl : substModTypes.get(type)) {
				String name = SubstitutionModel.getMenuName(cl);
				item = new JRadioButtonMenuItem(name);
				item.addActionListener(this);
				modelButtons[s++] = item;
				modelGroup.add(item);
				menu.add(item);
			}
		}
		modelGroup.clearSelection();
		menubar.add(menu);
		
		menu = new JMenu("Help");
		menu.setMnemonic(KeyEvent.VK_H);
		item = new JMenuItem("About...");
		item.addActionListener(this);
		menu.add(item);
		menu.addSeparator();
		item = new JMenuItem("Help for users");
		item.addActionListener(this);
		menu.add(item);
		menu.addSeparator();		
		item = new JMenuItem("Html doc for developers");
		item.addActionListener(this);
		menu.add(item);
		item = new JMenuItem("Description of plugins");
		item.addActionListener(this);
		menu.add(item);
		
		menubar.add(menu);
		
		Container cp = getContentPane();
		JPanel pan = new JPanel(new BorderLayout());
		Dimension screenSize = Toolkit.getDefaultToolkit().getScreenSize();
		pan.setMinimumSize(new Dimension(screenSize.width/3,screenSize.height/3));
		pan.setMaximumSize(new Dimension(screenSize.width,screenSize.height));
		pan.setPreferredSize(new Dimension(screenSize.width/2,screenSize.height/2));

		JTabbedPane tab = new JTabbedPane();
		Input input = new Input(manager);
		tab.addTab(input.getTabName(),input.getIcon(),input.getJPanel(),input.getTip());
		manager.inputgui = input.inputgui;
		for(Postprocess plugin : manager.postProcMan.plugins){
			if(plugin.selected){
				tab.addTab(plugin.getTabName(), plugin.getIcon(), plugin.getJPanel(), plugin.getTip());
			}
		}
		
		pan.add(tab, "Center");

		cp.add(pan, "Center");
		statusBar = new JPanel(new BorderLayout());
		statusText = new JLabel(IDLE_STATUS_MESSAGE);
		statusBar.add(statusText, BorderLayout.CENTER);
		cp.add(statusBar,BorderLayout.SOUTH);

		setJMenuBar(menubar);
//		setSize(300, 200);
//		setLocationByPlatform(true);
//    	setLocation(screenSize.width/4,screenSize.height/4);
		pack();
		setBounds(screenSize.width/5-15,screenSize.height/5-15,screenSize.width*3/5+30,screenSize.height*3/5+30);
		setVisible(true);
		setDefaultCloseOperation(EXIT_ON_CLOSE);
	}
	
	/**
	 * An ActioListener is implemented, so we have to implement this function. It handles actions on the menu bar.
	 * 
	 */
	public void actionPerformed(ActionEvent ev) {
		if(ev.getActionCommand() == "Add sequence(s)...") {
			JFileChooser choose = new JFileChooser("Add sequence(s)...");
			choose.setCurrentDirectory(new File(System.getProperty("user.dir")));
			if(choose.showOpenDialog(this) == JFileChooser.APPROVE_OPTION) {
				inFile = choose.getSelectedFile();
				FileFormatReader reader = new FastaReader();
				try {
					manager.seqs.add(reader.read(inFile));
					manager.inputgui.updateSequences();
					manager.fullPath = inFile.getAbsolutePath();
					if(manager.model != null) {
						try {
							manager.model.acceptable(manager.seqs);
						} catch (RecognitionError e) {
							tryModels();
						}
					} else {
						tryModels();
					}
					if(manager.model != null)
						runItem.setEnabled(true);
				} catch (IOException e) {
					JOptionPane.showMessageDialog(this,e.getLocalizedMessage(),"Error reading input file",JOptionPane.ERROR_MESSAGE);
				}
			}
		} else if(ev.getActionCommand() == "Exit") {
			System.exit(0);
		} else if(ev.getActionCommand() == "Preferences..."){
			//System.out.println("here!!!");
			new OutputPreferences(this);
		} else if(ev.getActionCommand() == "Settings") {
			mcmcSettingsDlg.display(this);
		} else if(ev.getActionCommand() == "Run") {
			if(manager.seqs.sequences.size() < 2){
				JOptionPane.showMessageDialog(this,"At least two sequences are needed!!!",
						"Not enough sequences",JOptionPane.ERROR_MESSAGE);
//				manager.finished();
				return;
			}
			openItem.setEnabled(false);
			runItem.setEnabled(false);
			pauseItem.setEnabled(true);
			resumeItem.setEnabled(false);
			stopItem.setEnabled(true);
			manager.start();
		} else if(ev.getActionCommand() == "Pause") {
			pauseItem.setEnabled(false);
			resumeItem.setEnabled(true);
			setCursor(Cursor.getPredefinedCursor(Cursor.WAIT_CURSOR));
			final String savTit = getTitle();
			setTitle("Pausing...");
			manager.thread.suspendSoft();
			setTitle(savTit);
			setCursor(Cursor.getDefaultCursor());
		} else if(ev.getActionCommand() == "Resume") {
			manager.thread.resumeSoft();
			pauseItem.setEnabled(true);
			resumeItem.setEnabled(false);
		} else if(ev.getActionCommand() == "Stop") {
			setCursor(Cursor.getPredefinedCursor(Cursor.WAIT_CURSOR));
			final String savTit = getTitle();
			setTitle("Stopping...");
			manager.thread.stopSoft();
			finished();
			setTitle(savTit);
			setCursor(Cursor.getDefaultCursor());
		} else if(ev.getActionCommand() == "About...") {
			new HelpWindow(this,"About", getClass().getClassLoader().getResource("doc/about/index.html"), false);
		}
		else if(ev.getActionCommand() == "Html doc for developers") {
			new HelpWindow(this,"Html doc for Developers", ClassLoader.getSystemResource("doc/index.html"), true);
		}
		else if(ev.getActionCommand() == "Description of plugins"){
			new HelpWindow(this,"Description of plugins", 
						ClassLoader.getSystemResource("doc/plugin_description/index.html"), true);
		
		}
		else if(ev.getActionCommand() == "Help for users"){
			new HelpWindow(this,"Help for users", 
					ClassLoader.getSystemResource("doc/help/index.html"), true);
		} else {		// new substitution model selected
			for(Class<? extends SubstitutionModel> cl : substModels) {
				try {
					if(ev.getActionCommand().equals(SubstitutionModel.getMenuName(cl))){
						try{
							SubstitutionModel model = cl.newInstance();
							model.acceptable(manager.seqs);
							manager.model = model;
							break;
						}
						catch(RecognitionError e){
							selectModel(manager.model);
							JOptionPane.showMessageDialog(this,e.message,"Cannot apply this model...",JOptionPane.ERROR_MESSAGE);
							break;
						}
					}
				} catch (Exception e) {
					new ErrorMessage(this,e.getLocalizedMessage(),true);
				}
			}

		}
	}

	private String tryModels() {
		String message = "";
		try {
			SubstitutionModel[] defaultSubstList = {
					new Kimura3(),
					new Dayhoff()
			};
			for(SubstitutionModel model : defaultSubstList) {
				try {
					model.acceptable(manager.seqs);
					selectModel(model);
					break;
				} catch (RecognitionError e) {
				}
			}
			if(manager.model == null) {
				double max = 0.0;
				int wrong = 0;
				for(Class<? extends SubstitutionModel> cl : substModels) {
					SubstitutionModel m;
					try {
						m = cl.newInstance();
						if(m.acceptable(manager.seqs) > max){
							manager.model = m;
							max = m.acceptable(manager.seqs);
							selectModel(m);
						}
					} catch (RecognitionError e){
						wrong++;
						message += e.message;
					}
				}
			}
		} catch (Exception e) {
			JOptionPane.showMessageDialog(this, e.getLocalizedMessage(),"Error accessing substitution models", JOptionPane.ERROR_MESSAGE);
		}
		return message;
		
	}

	private void selectModel(SubstitutionModel m) {
		manager.model = m;
		for(JMenuItem mi : modelButtons){
			if(mi.getText().equals(m.getMenuName())){
				mi.setSelected(true);
				break;
			}
		}
	}

	/**
	 * Enables several menu items that were kept disabled during the run.
	 */
   public void finished() {
		openItem.setEnabled(true);
		runItem.setEnabled(true);
		pauseItem.setEnabled(false);
		resumeItem.setEnabled(false);
		stopItem.setEnabled(false);
	}
	
   /**
    * Merely for testing purposes.
    * @param args No argument is used.
    */
	public static void main(String[] args) throws Exception{
		new MainFrame();
	}
}