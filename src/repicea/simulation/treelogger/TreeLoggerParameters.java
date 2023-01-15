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

import java.awt.Container;
import java.awt.Window;
import java.io.File;
import java.io.IOException;
import java.io.Serializable;
import java.lang.instrument.IllegalClassFormatException;
import java.security.InvalidParameterException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.TreeMap;

import javax.swing.filechooser.FileFilter;

import repicea.gui.REpiceaShowableUIWithParent;
import repicea.gui.permissions.DefaultREpiceaGUIPermission;
import repicea.gui.permissions.REpiceaGUIPermission;
import repicea.gui.permissions.REpiceaGUIPermissionProvider;
import repicea.io.IOUserInterfaceableObject;
import repicea.serial.Memorizable;
import repicea.serial.MemorizerPackage;
import repicea.serial.xml.PostXmlUnmarshalling;
import repicea.serial.xml.XmlDeserializer;
import repicea.serial.xml.XmlMarshallException;
import repicea.serial.xml.XmlMarshallingUtilities;
import repicea.serial.xml.XmlSerializer;
import repicea.simulation.treelogger.TreeLoggerParametersDialog.MessageID;
import repicea.util.ExtendedFileFilter;

/**
 * The TreeLoggerParameters contains the basic features for defining the
 * parameters of a treelogger.
 * @param <LC> a TreeLogCategory-derived class
 * @author M. Fortin - August 2010
 */
@SuppressWarnings("serial")
public abstract class TreeLoggerParameters<LC extends LogCategory>	implements Memorizable, 
																						IOUserInterfaceableObject, 
																						Serializable, 
																						REpiceaShowableUIWithParent, 
																						PostXmlUnmarshalling,
																						REpiceaGUIPermissionProvider {
	
	protected static class TreeLoggerParametersFileFilter extends FileFilter implements ExtendedFileFilter {

		private String extension = ".tlp";
		
		@Override
		public boolean accept(File file) {
			if (file.isDirectory()) { 
				return true; 
			} else {
				return file.getPath().toLowerCase().endsWith(getExtension());
			}
		}

		@Override
		public String getDescription() {
			return MessageID.TreeLoggerParametersFileExtension.toString();
		}

		@Override
		public String getExtension() {return extension;}
	}

	protected final static TreeLoggerParametersFileFilter TreeLoggerFileFilter = new TreeLoggerParametersFileFilter();

	public static final String ANY_SPECIES = "ANY";

	protected transient TreeLogger<?,?> treeLogger;		// MF20140207 changed to transient to avoid serializing when users save the parameters
	private String treeLoggerClass;
	private final Map<Object, List<LC>> selectedLogCategories;
	private transient boolean isParameterDialogCanceled;
	private String filename;
	private transient REpiceaGUIPermission readWrite = new DefaultREpiceaGUIPermission(true);
	
	/**
	 * General constructor all the TreeLoggerParameters-derived classes.
	 */
	protected TreeLoggerParameters(Class<? extends TreeLogger<?,?>> treeLoggerClass) {
		this.treeLoggerClass = treeLoggerClass.getName();
		isParameterDialogCanceled = false; // default value
		selectedLogCategories = new HashMap<Object, List<LC>>();
		setFilename("");
	}

	/**
	 * This method serves to create a default list of selected tree log
	 * categories. Each TreeLoggerParameter object must implements this method
	 * in order to provide its own default configuration.
	 */
	protected abstract void initializeDefaultLogCategories();

	/**
	 * This method can be override in order to specify a default species. That species cannot
	 * be removed from the list.
	 * @return an Object preferably either a String or an Enum
	 */
	protected Object getDefaultSpecies() {return null;}

	/**
	 * This method returns a specific LogCategory object among the selected
	 * tree log categories contained in this class.
	 * @param speciesName a String that represents the species name
	 * @param name the name of the tree log category
	 * @return a TreeLogCategory-derived instance or null if it cannot be found
	 */
	public LC getLogCategory(String speciesName, String name) {
		List<LC> logCategoryList = getSpeciesLogCategories(speciesName);
		if (logCategoryList != null) {
			for (LC logCategory : logCategoryList) {
				if (logCategory.getName().equals(name)) {
					return logCategory;
				}
			}
		}
		return null;
	}

	/**
	 * This method checks if the parameters are correct for the treelogger.
	 * @return a boolean
	 */
	public abstract boolean isCorrect();

	@Override
	public FileFilter getFileFilter() {return TreeLoggerFileFilter;}

	/**
	 * This method returns the list of the log category names. IMPORTANT: the list is
	 * filtered such that each value is unique. If the log category "sawlog" appears twice, 
	 * let's say once for balsam fir and once for black spruce, the list will only have one occurrence
	 * of "sawlog".
	 * @return a List of String
	 */
	public List<String> getLogCategoryNames() {
		List<String> outputList = new ArrayList<String>();
		for (List<LC> list : getLogCategories().values()) {
			for (LC logCategory : list) {
				if (!outputList.contains(logCategory.getName())) {
					outputList.add(logCategory.getName());
				}
			}
		}
		return outputList;
	}
	
	/**
	 * This method returns the List of LogCategory-derived object for a
	 * particular species. IMPORTANT: If there is a default species as specified through the
	 * getDefaultSpecies() method and the species is not found
	 * in the map, then the log categories of the default species are returned.
	 * @param species an Object that represents the species name
	 * @return a List of TreeLogCategory-derived instances or null if the species was not found
	 */
	public List<LC> getSpeciesLogCategories(Object species) {
		if (getLogCategories().containsKey(species)) {
			return getLogCategories().get(species);
		} else if (getDefaultSpecies() != null) {
			return getLogCategories().get(getDefaultSpecies());
		} else {
			return null;
		}
	}

	/**
	 * This method returns a Map instance with the species names as key and the
	 * corresponding List of LogCategory-derived instances as values.
	 * @return a Map instance
	 */
	public Map<Object, List<LC>> getLogCategories() {
		return selectedLogCategories;
	}
	
	/**
	 * This method returns a list of all the LogCategory instances.
	 * @return a List of LogCategory instances
	 */
	public List<LC> getLogCategoryList() {
		List<LC> outputList = new ArrayList<LC>();
		for (List<LC> list : getLogCategories().values()) {
			outputList.addAll(list);
		}
		return outputList;
	}

	/**
	 * This method returns false if the parameter initialization has been
	 * aborted in GUI mode. Returns true otherwise.
	 * 
	 * @return a boolean
	 */
	public boolean isParameterDialogCanceled() {
		return isParameterDialogCanceled;
	}

	/**
	 * This method makes it possible to specify whether or not the parameter
	 * initialization in GUI mode has been aborted.
	 */
	public void setParameterDialogCanceled(boolean b) {
		isParameterDialogCanceled = b;
	}

	@SuppressWarnings("rawtypes")
	@Override
	public MemorizerPackage getMemorizerPackage() {
		MemorizerPackage mp = new MemorizerPackage();
		mp.add(treeLoggerClass);
		mp.add(getFilename());
		Map logCategoryMap = getLogCategories();
		if (logCategoryMap instanceof TreeMap) { // former implementation
			mp.add((TreeMap) getLogCategories());		
		} else {
			mp.add((HashMap) getLogCategories());
		}
		return mp;
	}

	@SuppressWarnings({ "rawtypes", "unchecked" })
	@Override
	public void unpackMemorizerPackage(MemorizerPackage mp) {
		this.treeLoggerClass = mp.get(0).toString();
		setFilename(mp.get(1).toString());
		getLogCategories().clear();
		getLogCategories().putAll((Map) mp.get(2));
	}

	@Deprecated
	protected void fireEvent() {
		if (treeLogger != null) {
			if (treeLogger.wrapper != null) {
				treeLogger.wrapper.fireTreeLoggerEvent(treeLogger);
			}
		}
	}
	
	/**
	 * This method saves the parameters in the current filename.
	 * @throws IOException
	 */
	protected void save() throws IOException {
		save(getFilename());
	}
	
	@Override
	public void save(String filename) throws IOException {
		XmlSerializer serializer = new XmlSerializer(filename);
		try {
			serializer.writeObject(this);
			setFilename(filename);
			fireEvent();
		} catch (XmlMarshallException e) {
			e.printStackTrace();
			throw new IOException("Error while marshalling XML file!"); 
		}
	}


	public void load() throws IOException, IllegalClassFormatException {
		load(getFilename());
	}
	
	@Override
	public void load(String filename) throws IOException {
		TreeLoggerParameters<?> treeLoggerParameters = TreeLoggerParameters.loadFromFile(filename);
		if (!treeLoggerParameters.getClass().equals(this.getClass())) {
			throw new InvalidParameterException("The TreeLoggerParameters class in the file does not match the current TreeLoggerParameters class!");
		}
		unpackMemorizerPackage(treeLoggerParameters.getMemorizerPackage());
		fireEvent();
	}

	
	@Override
	public abstract TreeLoggerParametersDialog<?> getUI(Container parent);

	@Override
	public void showUI(Window parent) {
		getUI(parent).setVisible(true);
	}

	public String getFilename() {return filename;}

	public void setFilename(String filename) {this.filename = filename;}

	
	@SuppressWarnings("unchecked")
	public boolean equals(Object obj) {
		if (obj == null || !(obj instanceof TreeLoggerParameters)) {
			return false;
		} else {
	 		TreeLoggerParameters<LC> refParams = (TreeLoggerParameters<LC>) obj;
	 		if (!refParams.getClass().equals(getClass())) {
	 			return false;
	 		}
			if (!refParams.getFilename().equals(getFilename())) {
				return false;
			}

			if (refParams.getLogCategories().size() != getLogCategories().size()) {
				return false;
			}

			List<LC> refList;
			List<LC> thisList;
			LC refLogCategory;
			LC thisLogCategory;
			for (Object species : refParams.getLogCategories().keySet()) {
				refList = refParams.getSpeciesLogCategories(species);
				thisList = getSpeciesLogCategories(species);
				if (refList.size() != thisList.size()) {
					return false;
				} else {
					for (int i = 0; i < refList.size(); i++) {
						refLogCategory = refList.get(i);
						thisLogCategory = thisList.get(i);
						if (!refLogCategory.equals(thisLogCategory)) {
							return false;
						}
					}
				}
			}

			return true;
			
		}
	}
	
	/**
	 * This method returns a TreeLogger instance adapted to these parameters.
	 * @return a TreeLogger<?,?> instance
	 */
	@SuppressWarnings("unchecked")
	public TreeLogger<TreeLoggerParameters<?>,?> createTreeLoggerInstance() {
		try {
//			Class<?> clazz = ClassLoader.getSystemClassLoader().loadClass(treeLoggerClass);
			Class<?> clazz = Class.forName(treeLoggerClass);
			TreeLogger<TreeLoggerParameters<?>,?> treeLogger = (TreeLogger<TreeLoggerParameters<?>,?>) clazz.newInstance();
			treeLogger.setTreeLoggerParameters(this);
			return treeLogger;
		} catch (InstantiationException e) {
			e.printStackTrace();
			return null;
		} catch (IllegalAccessException e) {
			e.printStackTrace();
			return null;
		} catch (ClassNotFoundException e) {
			e.printStackTrace();
			return null;
		}
	}
	
	/**
	 * This method returns a TreeLoggerDescription based on the TreeLogger class of this
	 * instance.
	 * @return a TreeLoggerDescription instance
	 */
	public TreeLoggerDescription getTreeLoggerDescription() {
		return new TreeLoggerDescription(treeLoggerClass);
	}
	
	
	/**
	 * This method creates a TreeLoggerParameters instance from a previously saved .xml file
	 * @param filename the path of the xml file
	 * @return a TreeLoggerParameter<?> instance
	 * @throws IOException
	 */
	public static TreeLoggerParameters<?> loadFromFile(String filename) throws IOException {
		XmlDeserializer deserializer = new XmlDeserializer(filename);
		try {
			TreeLoggerParameters<?> treeLoggerParameters = (TreeLoggerParameters<?>) deserializer.readObject();
			treeLoggerParameters.setFilename(filename);
			return treeLoggerParameters;
		} catch (XmlMarshallException e) {
			throw new IOException("Error while unmarshalling the XML file!");
		}
	}

	
	@Override
	public void postUnmarshallingAction() {
		for (Object species : selectedLogCategories.keySet()) {
			for (LC logCategory : selectedLogCategories.get(species)) {
				logCategory.setSpecies(species);
			}
		}
		treeLoggerClass = XmlMarshallingUtilities.getClassName(treeLoggerClass);	// replace the tree logger class if it has been changed mean while
	}
	
	@Override
	public REpiceaGUIPermission getGUIPermission() {return readWrite;}

	/**
	 * This method sets the GUI rights for changing the tree logger parameters. By default, the rights are granted.
	 * @param permission a REpiceaGUIPermission instance
	 */
	public void setReadWritePermissionGranted(REpiceaGUIPermission permission) {this.readWrite = permission;}

	@Override
	public String toString() {
		int index = treeLoggerClass.lastIndexOf(".") + 1;
		return treeLoggerClass.substring(index, treeLoggerClass.length());
	}
	
}
