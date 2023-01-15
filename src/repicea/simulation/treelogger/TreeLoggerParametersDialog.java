/*
 * This file is part of the repicea-simulation library.
 *
 * Copyright (C) 2009-2012 Mathieu Fortin for Rouge-Epicea
 *
 * This library is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public
 * License as published by the Free Software Foundation; either
 * version 3 of the License, or (at your option) any later version.
 *
 * This library is distributed with the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied
 * warranty of MERCHANTABILITY or FITNESS FOR A
 * PARTICULAR PURPOSE. See the GNU Lesser General Public
 * License for more details.
 *
 * Please see the license at http://www.gnu.org/copyleft/lesser.html.
 */
package repicea.simulation.treelogger;

import java.awt.BorderLayout;
import java.awt.Component;
import java.awt.Dimension;
import java.awt.FlowLayout;
import java.awt.Insets;
import java.awt.Window;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.lang.reflect.ParameterizedType;
import java.util.ArrayList;
import java.util.List;
import java.util.Set;
import java.util.TreeMap;

import javax.swing.Box;
import javax.swing.BoxLayout;
import javax.swing.JButton;
import javax.swing.JLabel;
import javax.swing.JList;
import javax.swing.JMenu;
import javax.swing.JMenuBar;
import javax.swing.JMenuItem;
import javax.swing.JOptionPane;
import javax.swing.JPanel;
import javax.swing.JScrollPane;
import javax.swing.ListSelectionModel;
import javax.swing.SwingConstants;
import javax.swing.border.BevelBorder;
import javax.swing.event.ListSelectionEvent;
import javax.swing.event.ListSelectionListener;

import repicea.gui.CommonGuiUtility;
import repicea.gui.OwnedWindow;
import repicea.gui.REpiceaAWTProperty;
import repicea.gui.REpiceaDialog;
import repicea.gui.REpiceaMemorizerHandler;
import repicea.gui.Resettable;
import repicea.gui.UIControlManager;
import repicea.gui.UIControlManager.CommonControlID;
import repicea.gui.UIControlManager.CommonMenuTitle;
import repicea.gui.WindowSettings;
import repicea.gui.icons.IconFactory;
import repicea.io.IOUserInterface;
import repicea.io.REpiceaIOFileHandlerUI;
import repicea.lang.REpiceaSystem;
import repicea.serial.Memorizable;
import repicea.util.REpiceaTranslator;
import repicea.util.REpiceaTranslator.TextableEnum;

/**
 * The TreeLoggerParametersDialog is the abstract dialog that enables the specification of the different log grades or log categories. This class
 * makes it possible to add or remove species and log grades. The control over this feature can be overriden in derived class. 
 * @author Mathieu Fortin - February 2012
 * @param <P> a TreeLogCategory-extended class
 */
@SuppressWarnings("serial")
public abstract class TreeLoggerParametersDialog<P extends LogCategory> 
				extends REpiceaDialog implements ListSelectionListener, 
												ActionListener,
												OwnedWindow,
												IOUserInterface,
												Resettable {

	/**
	 * This enum variable contains the different messages displayed by the TreeLoggerParametersDialog classes.
	 * @author Mathieu Fortin - April 2012
	 */
	public static enum MessageID implements TextableEnum {
		LogGrade("Log Categories", 
				"Cat\u00E9gories de billon"),
		Species("Species", 
				"Esp\u00E8ces"),
		PleaseEnterTheSpeciesCode("Please enter the species code (3 characters): ",
				"Veuillez saisir le code \u00E0 3 caract\u00E8res de l'esp\u00E8ce :"),
		AddANewSpecies("Add a new species",
				"Ajouter une nouvelle esp\u00E8ce"),
		TheFollowingSpeciesIsAlreadyIncludedInTheList("The following species is already included in the list: ",
				"L'esp\u00E8ce suivante est d\u00E9j\u00E0 comprise dans la liste :"),
		ParamsHaveChanged("The parameters of the tree logger have been changed. Do you want to save them before executing the tree logger task?",
				"Les param\u00E8tres du module de billonnage ont \u00E9t\u00E9 chang\u00E9s. D\u00E9sirez-vous les sauvegarder avant d'ex\u00E9cuter la t\u00E2che de billonnage?"),
		DefaultValues("Default values", "Valeurs par d\u00E9fault"),
		ErrorDifferentParametersClass("An error occurred! The parameters in the file are from a different tree logger.", "Une erreur est survenue! Les param\u00E8tres lus dans le fichier correspondent \u00E0 un module de billonnage diff\u00E9rent."),
		TreeLoggerParametersFileExtension("tree logger parameters file (*.tlp)", "Fichier des param\u00E8tres de billonnage (*.tlp)");
		
		MessageID(String englishText, String frenchText) {
			setText(englishText, frenchText);
		}
		
		public void setText(String englishText, String frenchText) {
			REpiceaTranslator.setString(this, englishText, frenchText);
		}

		@Override
		public String toString() {
			return REpiceaTranslator.getString(this);
		}
	}
	
	
	/*
	 * Just a JList that makes sure the selected index is always visible
	 */
	@SuppressWarnings("rawtypes")
	protected static class CustomJList extends JList {
		@Override
		public void setSelectedIndex(int index) {
			super.setSelectedIndex(index);
			ensureIndexIsVisible(index);
		}
	}
	
		
	protected TreeLoggerParameters<P> params;

	protected JPanel panLogCategory;
	protected CustomJList logCategoryList;	
	protected CustomJList speciesList;				
	
	protected JButton ok;
	protected JButton cancel;
	
	protected JMenu mnFile; 
	private JMenuItem save;
	private JMenuItem saveAs;
	private JMenuItem load;
	
	protected JMenu mnEdit;
	private JMenuItem undo;
	private JMenuItem redo;
	private JMenuItem defaultValues;
	
	protected JMenu mnSpecies; 
	private JMenuItem speciesAdd;
	private JMenuItem speciesRemove;
	
	protected JMenu mnLogGrade;
	private JMenuItem logGradeAdd;
	private JMenuItem logGradeRemove;
	
	protected JMenu mnTools;
	private JMenuItem settingsMenuItem;
	
	protected JMenu mnHelp;
	protected JMenuItem help;

	protected JButton logGradeGoUp;
	protected JButton logGradeGoDown;
	
	protected JPanel controlPanel;
	
	protected boolean logGradePriorityChangeEnabled = true;  // default value
	
	protected final WindowSettings windowSettings;
	
	/**
	 * Constructor. To be called by derived classes.
	 * @param window the parent window
	 * @param params a TreeLoggerParameters instance
	 */
	protected TreeLoggerParametersDialog(Window window, TreeLoggerParameters<P> params) {
		super(window);
		windowSettings = new WindowSettings(REpiceaSystem.getJavaIOTmpDir() + getClass().getSimpleName()+ ".ser", this);
		instantiateVariables(params);
		new REpiceaIOFileHandlerUI(this, params, save, saveAs, load);
		new REpiceaMemorizerHandler(this, undo, redo);
		synchronizeUIWithOwner();
		initUI();
	}

	
	protected void instantiateVariables(TreeLoggerParameters<P> params) {
		this.params = params;
		
		panLogCategory = new JPanel(new BorderLayout());
		speciesList = new CustomJList();
		speciesList.setSelectionMode(ListSelectionModel.SINGLE_SELECTION);
		speciesList.setBorder(new BevelBorder(BevelBorder.LOWERED, null, null, null, null));
		speciesList.setVisibleRowCount(8);
		speciesList.setLayoutOrientation(JList.VERTICAL);
		
		save = UIControlManager.createCommonMenuItem(CommonControlID.Save);
		saveAs = UIControlManager.createCommonMenuItem(CommonControlID.SaveAs);
		load = UIControlManager.createCommonMenuItem(CommonControlID.Open);
		
		redo = UIControlManager.createCommonMenuItem(CommonControlID.Redo);
		undo = UIControlManager.createCommonMenuItem(CommonControlID.Undo);
		defaultValues = new JMenuItem(REpiceaTranslator.getString(MessageID.DefaultValues));

		speciesAdd = UIControlManager.createCommonMenuItem(CommonControlID.Addone);
		speciesRemove = UIControlManager.createCommonMenuItem(CommonControlID.RemoveOne);
		
		logGradeAdd = UIControlManager.createCommonMenuItem(CommonControlID.Addone);
		logGradeRemove = UIControlManager.createCommonMenuItem(CommonControlID.RemoveOne);

		logGradeGoUp = new JButton();
		logGradeGoUp.setIcon(IconFactory.getIcon(CommonControlID.GoUp));
		logGradeGoUp.setMargin(new Insets(1,1,1,1));
		logGradeGoUp.setEnabled(getTreeLoggerParameters().getGUIPermission().isEnablingGranted());
		
		logGradeGoDown = new JButton();
		logGradeGoDown.setIcon(IconFactory.getIcon(CommonControlID.GoDown));
		logGradeGoDown.setMargin(new Insets(1,1,1,1));
		logGradeGoDown.setEnabled(getTreeLoggerParameters().getGUIPermission().isEnablingGranted());


		logCategoryList = new CustomJList();
		logCategoryList.setSelectionMode(ListSelectionModel.SINGLE_SELECTION);
		logCategoryList.setBorder(new BevelBorder(BevelBorder.LOWERED, null, null, null, null));
		logCategoryList.setVisibleRowCount(8);
		logCategoryList.setLayoutOrientation(JList.VERTICAL);

		settingsMenuItem = UIControlManager.createCommonMenuItem(CommonControlID.Options);
		
		help = UIControlManager.createCommonMenuItem(CommonControlID.Help);
		ok = UIControlManager.createCommonButton(CommonControlID.Ok);
		cancel = UIControlManager.createCommonButton(CommonControlID.Cancel);
		ok.setEnabled(getTreeLoggerParameters().getGUIPermission().isEnablingGranted());
		cancel.setEnabled(getTreeLoggerParameters().getGUIPermission().isEnablingGranted());
	}
		

	@Override
	public void okAction() {
		params.setParameterDialogCanceled(false);
		super.okAction();
	}
	
//	protected abstract String getTreeLoggerName();

	protected JPanel getControlPanel() {
		if (controlPanel == null) {
			controlPanel = new JPanel(new FlowLayout(FlowLayout.RIGHT));
			controlPanel.add(ok);
			controlPanel.add(cancel);
		}
		return controlPanel;
	}
	
	protected List<Object> provideSpeciesInNaturalOrder() {
		TreeMap<String, Object> names = new TreeMap<String, Object>();
		Set<Object> speciesSet = params.getLogCategories().keySet();
		for (Object species : speciesSet) {
			names.put(species.toString(), species);
		}
		List<Object> sortedSpeciesList = new ArrayList<Object>();
		sortedSpeciesList.addAll(names.values());
		Object defaultSpecies = params.getDefaultSpecies();
		if (defaultSpecies != null) {
			if (sortedSpeciesList.remove(defaultSpecies)) {		// we remove the default species
				sortedSpeciesList.add(defaultSpecies);		// we add the default species at the end of the list
			};
		}
		return sortedSpeciesList;
	}

	@SuppressWarnings("unchecked")
	protected final void redefineSpeciesList() {
		speciesList.removeListSelectionListener(this);
		
//		speciesList.setListData(params.getLogCategories().keySet().toArray());
		speciesList.setListData(provideSpeciesInNaturalOrder().toArray());
		
		
		speciesList.addListSelectionListener(this);
		if (speciesList.getModel().getSize() <= 1) {
			speciesRemove.setEnabled(false);
		} else {
			speciesRemove.setEnabled(true);
		}

	}
	
	@Override
	protected void initUI() {
		getContentPane().setLayout(new BorderLayout());
		

		getContentPane().add(getControlPanel(), BorderLayout.SOUTH);

		JPanel listPanel = new JPanel(new FlowLayout(FlowLayout.LEFT));
		
		listPanel.add(Box.createHorizontalStrut(10));
		
		JPanel speciesListPanel = new JPanel();
		speciesListPanel.setLayout(new BoxLayout(speciesListPanel, BoxLayout.Y_AXIS));
		speciesListPanel.add(Box.createVerticalStrut(10));
		JLabel speciesLabel = new JLabel(REpiceaTranslator.getString(MessageID.Species));
		speciesLabel.setAlignmentX(0.5f);
		speciesLabel.setHorizontalAlignment(SwingConstants.LEFT);
		speciesListPanel.add(speciesLabel);
		speciesListPanel.add(Box.createVerticalStrut(10));
		
		JScrollPane speciesListScrollPane = new JScrollPane(speciesList);
		speciesListScrollPane.setWheelScrollingEnabled(true);
		
		speciesListPanel.add(speciesListScrollPane);
		speciesList.setAlignmentX(0.5f);
		speciesListPanel.setAlignmentY(1);
		
		JPanel logCategoryListPanel = new JPanel();
		logCategoryListPanel.setLayout(new BoxLayout(logCategoryListPanel, BoxLayout.Y_AXIS));
		logCategoryListPanel.add(Box.createVerticalStrut(10));
		JLabel logCategoryLabel = new JLabel(REpiceaTranslator.getString(MessageID.LogGrade));
		logCategoryLabel.setAlignmentX(0.5f);
		logCategoryLabel.setHorizontalAlignment(SwingConstants.LEFT);
		logCategoryListPanel.add(logCategoryLabel);
		logCategoryListPanel.add(Box.createVerticalStrut(10));
		
		
		JScrollPane logCategoriesScrollPane = new JScrollPane(logCategoryList);
		logCategoriesScrollPane.setWheelScrollingEnabled(true);
		logCategoryListPanel.add(logCategoriesScrollPane);
		logCategoryList.setAlignmentX(0.5f);
		logCategoryListPanel.setAlignmentY(1);
		
		
		listPanel.add(speciesListPanel);
		listPanel.add(Box.createHorizontalStrut(10));
		listPanel.add(logCategoryListPanel);

		JPanel listPanelWrapper = new JPanel(new FlowLayout(FlowLayout.LEFT));
		listPanelWrapper.add(listPanel);
		
		JPanel panel = new JPanel();
		listPanel.add(panel);
		panel.setLayout(new BoxLayout(panel, BoxLayout.Y_AXIS));
		
		panel.add(logGradeGoUp);
	
		panel.add(logGradeGoDown);
		
		Component horizontalStrut = Box.createHorizontalStrut(10);
		listPanel.add(horizontalStrut);
		
		getContentPane().add(listPanelWrapper, BorderLayout.WEST);
		getContentPane().add(panLogCategory, BorderLayout.CENTER);
		
		setMenuBar();
		
		
		pack();
		
		Dimension minDim = new Dimension(600, 350);
		setMinimumSize(minDim);
		setSize(minDim);
		
		ok.setDefaultCapable(true);
		getRootPane().setDefaultButton(ok);
	}


	protected void setMenuBar() {
		JMenuBar menuBar = new JMenuBar();
		setJMenuBar(menuBar);
		mnFile = UIControlManager.createCommonMenu(CommonMenuTitle.File);
		menuBar.add(mnFile);
		mnFile.add(load);
		mnFile.add(save);
		mnFile.add(saveAs);
		mnFile.setEnabled(getTreeLoggerParameters().getGUIPermission().isEnablingGranted());

		mnEdit = UIControlManager.createCommonMenu(CommonMenuTitle.Edit);
		menuBar.add(mnEdit);
		mnEdit.add(redo);
		mnEdit.add(undo);
		mnEdit.add(defaultValues);
		mnEdit.setEnabled(getTreeLoggerParameters().getGUIPermission().isEnablingGranted());
		
		mnSpecies = new JMenu(REpiceaTranslator.getString(MessageID.Species));
		menuBar.add(mnSpecies);
		mnSpecies.add(speciesAdd);
		mnSpecies.add(speciesRemove);
		mnSpecies.setEnabled(getTreeLoggerParameters().getGUIPermission().isEnablingGranted());
		
		mnLogGrade = new JMenu(REpiceaTranslator.getString(MessageID.LogGrade));
		menuBar.add(mnLogGrade);
		mnLogGrade.add(logGradeAdd);
		mnLogGrade.add(logGradeRemove);
		mnLogGrade.setEnabled(getTreeLoggerParameters().getGUIPermission().isEnablingGranted());
		
		mnTools = UIControlManager.createCommonMenu(CommonMenuTitle.Tools);
		menuBar.add(mnTools);
		mnTools.add(settingsMenuItem);
		mnTools.setEnabled(getTreeLoggerParameters().getGUIPermission().isEnablingGranted());
		
		mnHelp = UIControlManager.createCommonMenu(CommonMenuTitle.About);
		menuBar.add(mnHelp);
		mnHelp.add(help);
	}
	
	@SuppressWarnings("unchecked")
	@Override
	public void valueChanged(ListSelectionEvent evt) {
		if (evt.getValueIsAdjusting()) {
			return;
		} else if (evt.getSource().equals(speciesList)) {
			defineLogCategoryList();
		} else if (evt.getSource().equals(logCategoryList)) {
			int selectedIndex = logCategoryList.getSelectedIndex();
			if (logGradePriorityChangeEnabled) {
				logGradeGoUp.setEnabled(true && getTreeLoggerParameters().getGUIPermission().isEnablingGranted());
				logGradeGoDown.setEnabled(true && getTreeLoggerParameters().getGUIPermission().isEnablingGranted());
				if (selectedIndex == 0) {
					logGradeGoUp.setEnabled(false);
				}
				if (selectedIndex == logCategoryList.getModel().getSize() - 1) {
					logGradeGoDown.setEnabled(false);
				}
			}
			P logCategory = (P) logCategoryList.getSelectedValue();
			panLogCategory.removeAll();
			panLogCategory.add(logCategory.getUI(), BorderLayout.CENTER);
			boolean enabled = getTreeLoggerParameters().getGUIPermission().isEnablingGranted();
			CommonGuiUtility.enableAllControls(panLogCategory, enabled);
			panLogCategory.revalidate();
			panLogCategory.repaint();
		}
	}

	@SuppressWarnings("unchecked")
	private void defineLogCategoryList() {
		Object species = speciesList.getSelectedValue();
		logCategoryList.removeListSelectionListener(this);
		logCategoryList.setListData(params.getLogCategories().get(species).toArray());
		logCategoryList.addListSelectionListener(this);
		logCategoryList.setSelectedIndex(0);
		if (logCategoryList.getModel().getSize() <= 1) {
			logGradeRemove.setEnabled(false);
		} else {
			logGradeRemove.setEnabled(true);
		}
		if (speciesList.getModel().getSize() >= 2) {
			speciesRemove.setEnabled(!species.equals(params.getDefaultSpecies()));	// if the species list has two elements or more the species remove button is disabled if the selected species is the default species
		}
	}
	
	private void logGradeGoUpAction() {
		int selectedIndex = logCategoryList.getSelectedIndex();
		Object species = speciesList.getSelectedValue();
		List<P> logList = params.getLogCategories().get(species);
		P removedElement = logList.remove(selectedIndex - 1);
		logList.add(selectedIndex, removedElement);
		defineLogCategoryList();
		logCategoryList.setSelectedIndex(selectedIndex - 1);		//keep track of the selected index
	}

	private void logGradeGoDownAction() {
		int selectedIndex = logCategoryList.getSelectedIndex();
		Object species = speciesList.getSelectedValue();
		List<P> logList = params.getLogCategories().get(species);
		P removedElement = logList.remove(selectedIndex);
		logList.add(selectedIndex + 1, removedElement);
		defineLogCategoryList();
		logCategoryList.setSelectedIndex(selectedIndex + 1);		//keep track of the selected index
	}


	protected void undoAction() {
		ok.requestFocus();			// put the focus on ok to avoid having a field with the focus (otherwise the refresh might not be effective on all components)
		firePropertyChange(REpiceaAWTProperty.UndoClicked, null, this);
		refreshInterface();
	}

	@Override
	public void reset() {
		params.initializeDefaultLogCategories();
		synchronizeUIWithOwner();
		refreshInterface();
	}
	
	@Override
	public void cancelAction() {
		params.setParameterDialogCanceled(true);
		super.cancelAction();	
	}
	

	@Override
	public void listenTo() {
		speciesList.setSelectedIndex(0);
		undo.addActionListener(this);
		defaultValues.addActionListener(this);
		cancel.addActionListener(this);
		ok.addActionListener(this);
		settingsMenuItem.addActionListener(this);
		help.addActionListener(this);
		logGradeGoDown.addActionListener(this);
		logGradeGoUp.addActionListener(this);
		logGradeAdd.addActionListener(this);
		logGradeRemove.addActionListener(this);
		speciesAdd.addActionListener(this);
		speciesRemove.addActionListener(this);
	}

	@Override
	public void doNotListenToAnymore() {
//		speciesList.removeListSelectionListener(this);
//		logCategoryList.removeListSelectionListener(this);
		undo.removeActionListener(this);
		defaultValues.addActionListener(this);
		cancel.removeActionListener(this);
		ok.removeActionListener(this);
		settingsMenuItem.removeActionListener(this);
		help.removeActionListener(this);
		logGradeGoDown.removeActionListener(this);
		logGradeGoUp.removeActionListener(this);
		logGradeAdd.removeActionListener(this);
		logGradeRemove.removeActionListener(this);
		speciesAdd.removeActionListener(this);
		speciesRemove.removeActionListener(this);
	}


	@Override
	public void actionPerformed(ActionEvent arg0) {
		if (arg0.getSource().equals(logGradeGoDown)) {
			logGradeGoDownAction();
		} else if (arg0.getSource().equals(logGradeGoUp)) { 
			logGradeGoUpAction();
		} else if (arg0.getSource().equals(speciesAdd)) {
			speciesAddAction();
		} else if (arg0.getSource().equals(speciesRemove)) {
			speciesRemoveAction();
		} else if (arg0.getSource().equals(logGradeAdd)) {
			logGradeAddAction();
		} else if (arg0.getSource().equals(logGradeRemove)) {
			logGradeRemoveAction();
		} else if (arg0.getSource().equals(undo)) {
			undoAction();
		} else if (arg0.getSource().equals(defaultValues)) {
			reset();
		} else if (arg0.getSource().equals(cancel)) {
			cancelAction();
		} else if (arg0.getSource().equals(ok)) {
			okAction();
		} else if (arg0.getSource().equals(help)) {
			helpAction();
		} else if (arg0.getSource().equals(settingsMenuItem)) {
			settingsAction();
		}
		
	}

	protected abstract void settingsAction();

	@Override
	public void postLoadingAction() {
		synchronizeUIWithOwner();
	}
	
	@Override
	public void postSavingAction() {
		REpiceaIOFileHandlerUI.RefreshTitle(params, this);
	}
	
	@Override
	public WindowSettings getWindowSettings() {return windowSettings;}
	
	private void logGradeRemoveAction() {
		int selectedIndex = logCategoryList.getSelectedIndex();
		Object species = speciesList.getSelectedValue();
		params.getLogCategories().get(species).remove(selectedIndex);
		defineLogCategoryList();
		if (selectedIndex < logCategoryList.getModel().getSize()) {
			logCategoryList.setSelectedIndex(selectedIndex);
		} else {
			logCategoryList.setSelectedIndex(logCategoryList.getModel().getSize() - 1);
		}
		firePropertyChange(TreeLoggerAWTProperty.LogGradeRemoved, null, this);
	}

	@SuppressWarnings("unchecked")
	private void logGradeAddAction() {
		Object species = speciesList.getSelectedValue();
		List<P> logGrades = params.getLogCategories().get(species);
		Class<? extends LogCategory> clazz = (Class<? extends LogCategory>) ((ParameterizedType) this.getClass().getGenericSuperclass()).getActualTypeArguments()[0];
		try {
			P newGrade = (P) clazz.newInstance();
			newGrade.setSpecies(species);
			logGrades.add(newGrade);
		} catch (InstantiationException e) {
			e.printStackTrace();
		} catch (IllegalAccessException e) {
			e.printStackTrace();
		}
		defineLogCategoryList();
		logCategoryList.setSelectedIndex(logGrades.size() - 1);
		firePropertyChange(TreeLoggerAWTProperty.LogGradeAdded, null, this);
	}

	private void speciesRemoveAction() {
		int selectedIndex = speciesList.getSelectedIndex();
		Object species = speciesList.getSelectedValue();
		params.getLogCategories().remove(species);
		redefineSpeciesList();
		if (selectedIndex < speciesList.getModel().getSize()) {
			speciesList.setSelectedIndex(selectedIndex);
		} else {
			speciesList.setSelectedIndex(speciesList.getModel().getSize() - 1);
		}		
		defineLogCategoryList();
		firePropertyChange(TreeLoggerAWTProperty.SpeciesRemoved, null, this);
	}

	/**
	 * The basic implementation of this method assumes that the species are designated by String instance.
	 */
	@SuppressWarnings("unchecked")
	protected void speciesAddAction() {
		Object species = speciesList.getSelectedValue();
		List<P> logGrades = params.getLogCategories().get(species);
		Class<? extends LogCategory> clazz = logGrades.get(0).getClass();
		String newSpecies;
		boolean valid = false;
		do {
			newSpecies = JOptionPane.showInputDialog(this, 
					REpiceaTranslator.getString(MessageID.PleaseEnterTheSpeciesCode),
					REpiceaTranslator.getString(MessageID.AddANewSpecies),
					JOptionPane.QUESTION_MESSAGE);
			if (newSpecies == null) {		// means it has been canceled
				return;
			}
			if (newSpecies.length() >= 3) {
				newSpecies = newSpecies.substring(0, 3).toUpperCase();
				if (params.getLogCategories().keySet().contains(newSpecies)) {
					JOptionPane.showMessageDialog(this, 
							REpiceaTranslator.getString(MessageID.TheFollowingSpeciesIsAlreadyIncludedInTheList) + newSpecies,
							REpiceaTranslator.getString(UIControlManager.InformationMessageTitle.Error), 
							JOptionPane.ERROR_MESSAGE);
				} else {
					valid = true;
				}
			}
		} while (!valid);
		
		List<P> newLogGrades = new ArrayList<P>();
		try {
			newLogGrades.add((P) clazz.newInstance());
			params.getLogCategories().put(newSpecies, newLogGrades);
		} catch (InstantiationException e) {
			e.printStackTrace();
		} catch (IllegalAccessException e) {
			e.printStackTrace();
		}
		redefineSpeciesList();
		int selectedIndex = 0;
		for (Object speciesName : params.getLogCategories().keySet()) {
			if (speciesName.equals(newSpecies)) {
				speciesList.setSelectedIndex(selectedIndex);
				return;
			}
			selectedIndex++;
		}
		speciesList.setSelectedIndex(0);		// in case no match has been found
		firePropertyChange(TreeLoggerAWTProperty.SpeciesAdded, null, this);
	}

	@Override
	public Memorizable getWindowOwner() {return getTreeLoggerParameters();}

	protected TreeLoggerParameters<P> getTreeLoggerParameters() {return params;}

	
	@Override
	public void synchronizeUIWithOwner() {
		redefineSpeciesList();
		speciesList.setSelectedIndex(0);
		defineLogCategoryList();
		REpiceaIOFileHandlerUI.RefreshTitle(params, this);
	}
	
}
