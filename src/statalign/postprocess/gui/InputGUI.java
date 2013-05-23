package statalign.postprocess.gui;

import javax.swing.*;
import javax.swing.border.EtchedBorder;
import javax.swing.event.ListSelectionEvent;
import javax.swing.event.ListSelectionListener;

import java.awt.*;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.util.ArrayList;

import statalign.base.Input;
import statalign.base.MainManager;

/**
 * This is the graphical interface for showing the input data
 * 
 * @author miklos, novak
 *
 */
public class InputGUI extends JPanel implements ActionListener, ListSelectionListener{

	/**
	 * 
	 */
	private static final long serialVersionUID = 1L;
	JPanel pan;
	Input owner;
	TextArea text;
	MainManager manager;
	
	private JList sequences;
	DefaultListModel dlmSequences;
	JButton jbDelete;
	
	/**
	 * This constructor makes an initial GUI for showing the input sequences and their names.
	 * @param manager The MainManager that handles the MCMC run.
	 */
	public InputGUI(MainManager manager){
		super(new BorderLayout());
		this.manager = manager;
		dlmSequences = new DefaultListModel();
		sequences = new JList(dlmSequences);
		sequences.setBorder(new EtchedBorder());
		sequences.setToolTipText("Input sequences - click on them to view or remove");
		sequences.addListSelectionListener(this);
//		JPanel intermediatePanel = new JPanel(new BorderLayout());
//		intermediatePanel.add(sequences);
//		intermediatePanel.setSize(this.getSize());
//		intermediatePanel.setMaximumSize(this.getSize());
		JScrollPane spSzoveg = new JScrollPane(sequences);//intermediatePanel, JScrollPane.VERTICAL_SCROLLBAR_ALWAYS, JScrollPane.HORIZONTAL_SCROLLBAR_ALWAYS);
		setPreferredSize(new Dimension(500,300));
//		spSzoveg.getViewport().add(sequences);
//		spSzoveg.setMaximumSize(this.getSize());
		add(spSzoveg,"Center");
		
		JPanel actionPanel = new JPanel(new BorderLayout());
		jbDelete = new JButton("Remove");
		jbDelete.addActionListener(this);
		actionPanel.add(jbDelete,"West");
		add(actionPanel,"South");

		updateSequences();
	}
	
	/**
	 * It rereads sequences from MainManager
	 */
	public void updateSequences(){
		if(dlmSequences.size() > 0){
			dlmSequences.removeAllElements();
		}
		if(manager.seqs != null){
//			System.out.println("sequences size: "+manager.seqs.sequences.size()+
//					" names size: "+manager.seqs.seqNames.size());		    
			for(int i = 0; i < manager.seqs.sequences.size(); i++){
				String s1 = manager.seqs.sequences.get(i);
				String s2 = "";
				int length = 0;
				for(int j = 0; j < s1.length(); j++){
					if(s1.charAt(j) != ' ' && s1.charAt(j) != '-'){
						s2 += s1.charAt(j);
						length++;
						if(length % 60 == 0){
							s2+="<br>";
						}
					}
				}
				dlmSequences.addElement("<html><font color=\"000099\">&gt;"+manager.seqs.seqNames.get(i)+"</font>\n<br><font face=\"Courier New\">"+s2+"</font></html>");
		    }
		}
		listListener();
	}
	
/*	
	public InputGUI(JPanel pan, Input inp, String s){
		text = new TextArea(s,pan.getHeight()/13,(int)(pan.getWidth()/6.6),TextArea.SCROLLBARS_BOTH);
		this.pan = pan;
		this.owner = inp;
	  	text.setFont(new Font("Monospaced",Font.PLAIN,10));
	  	text.setEditable(false);
	  	add(text);
	}
	
	  public void paintComponent(Graphics gr){
		  super.paintComponent(gr);
		  text.setColumns((int)(pan.getWidth()/6.6));
		  text.setRows(pan.getHeight()/13);
	  }
*/
	
	void listListener(){
		int index = sequences.getSelectedIndex();
		if(index == -1){
			jbDelete.setEnabled(false);
			updateUI();
		}
		else{
			jbDelete.setEnabled(true);
			updateUI();
		}
	}

	/**
	 * Handles removing sequences.
	 */
	public void actionPerformed(ActionEvent arg0) {
	    int index = sequences.getSelectedIndex();
	    dlmSequences.remove(index);
	    if(dlmSequences.getSize() != 0){
	    	if(index == dlmSequences.getSize()){
	    		index--;
	    	}
	    	sequences.setSelectedIndex(index);
	    }
	    listListener();
	    manager.seqs.seqNames = new ArrayList<String>();
	    manager.seqs.sequences = new ArrayList<String>();
	    for(int i = 0; i < dlmSequences.getSize(); i++){
	    	String s = (String) dlmSequences.get(i);
	    	manager.seqs.seqNames.add(new String(s.substring(31,s.indexOf('\n')-7)));
	    	//System.out.println(s.substring(31,s.indexOf('\n')-7));
	    	
	    	manager.seqs.sequences.add(new String(s.substring(s.indexOf('\n')+30, s.length()-14)).replaceAll("<br>", ""));
	    	//System.out.println(s.substring(s.indexOf('\n')+30, s.length()-14).replaceAll("<br>", ""));
	    }
	}

	/**
	 * It invokes the list listener when a value changed.
	 */
	public void valueChanged(ListSelectionEvent arg0) {
		listListener();
		
	}

}
