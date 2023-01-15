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
package repicea.simulation.processsystem;

import java.awt.Container;
import java.awt.Window;
import java.io.File;
import java.io.IOException;
import java.io.InputStream;
import java.util.ArrayList;
import java.util.List;

import javax.swing.filechooser.FileFilter;

import repicea.gui.ListManager;
import repicea.gui.REpiceaShowableUIWithParent;
import repicea.gui.Resettable;
import repicea.gui.permissions.DefaultREpiceaGUIPermission;
import repicea.gui.permissions.REpiceaGUIPermission;
import repicea.gui.permissions.REpiceaGUIPermissionProvider;
import repicea.io.IOUserInterfaceableObject;
import repicea.serial.Memorizable;
import repicea.serial.MemorizerPackage;
import repicea.serial.xml.XmlDeserializer;
import repicea.serial.xml.XmlMarshallException;
import repicea.serial.xml.XmlSerializer;
import repicea.util.ObjectUtility;

public class SystemManager implements ListManager<Processor>, 
										IOUserInterfaceableObject, 
										Resettable, 
										Memorizable, 
										REpiceaShowableUIWithParent,
										REpiceaGUIPermissionProvider {

	private final List<Processor> processors;

	private String filename;
	
	private transient REpiceaGUIPermission readWrite;
	
	protected transient SystemManagerDialog guiInterface;
	
	public SystemManager(REpiceaGUIPermission permission) {
		readWrite = permission;
		processors = new ArrayList<Processor>();
		setFilename(System.getProperty("user.home") + File.separator + SystemManagerDialog.MessageID.Unnamed.toString());
	}

	public SystemManager() {
		this(new DefaultREpiceaGUIPermission(true));
	}

	
	@Override
	public List<Processor> getList() {
		return processors;
	}

	@Override
	public void registerObject(Processor process) {
		if (!processors.contains(process)) {
			processors.add(process);
		}
	}

	@Override
	public void removeObject(Processor process) {
		processors.remove(process);
	}


	@Override
	public SystemManagerDialog getUI(Container parent) {
		if (guiInterface == null) {
			guiInterface = new SystemManagerDialog((Window) parent, this);
		}
		return guiInterface;
	}

	@Override
	public void showUI(Window parent) {
		getUI(parent).setVisible(true);
	}

	@Override
	public void load(String filename) throws IOException {
		XmlDeserializer deserializer;
		try {
			deserializer = new XmlDeserializer(filename);
		} catch (Exception e) {
//			InputStream is = ClassLoader.getSystemResourceAsStream(filename);
			InputStream is = getClass().getResourceAsStream("/" + filename);
			if (is == null) {
				throw new IOException("The filename is not a file and cannot be converted into a stream!");
			} else {
				deserializer = new XmlDeserializer(is);
			}
		}
		SystemManager newManager;
		try {
			newManager = (SystemManager) deserializer.readObject();
			newManager.setFilename(filename);
			unpackMemorizerPackage(newManager.getMemorizerPackage());
		} catch (XmlMarshallException e) {
			throw new IOException("A XmlMarshallException occurred while loading the file!");
		}
	}

	
	@Override
	public void save(String filename) throws IOException {
		setFilename(filename);
		XmlSerializer serializer = new XmlSerializer(filename);
		try {
			serializer.writeObject(this);
		} catch (XmlMarshallException e) {
			throw new IOException("A XmlMarshallException occurred while saving the file!");
		}
	}


	/**
	 * This method returns the list of primary processors, ie. those who are the sons of no one.
	 * @return a List of Processor instances
	 */
	public List<Processor> getPrimaryProcessors() {
		List<Processor> primaryProcessors = new ArrayList<Processor>();
		for (Processor processor : processors) {
			boolean isSonOfSomeOne = false;
			for (Processor otherProcessor : processors) {
				if (otherProcessor.getSubProcessorIntakes().containsKey(processor)) {
					isSonOfSomeOne = true;
					break;
				}
			}
			if (!isSonOfSomeOne) {
				primaryProcessors.add(processor);
			}
		}
		return primaryProcessors;
	}
	

	/*
	 * To be refined in derived class (non-Javadoc)
	 * @see repicea.simulation.Parameterizable#getFileFilter()
	 */
	@Override
	public FileFilter getFileFilter() {
		return null;
	}

	private void setFilename(String filename) {this.filename = filename;}
	
	@Override
	public String getFilename() {return filename;}

	public String getName() {
		File file = new File(filename);
		return ObjectUtility.relativizeTheseFile(file.getParentFile(), file).toString();
	}

	@Override
	public void reset() {
		processors.clear();
		setFilename(System.getProperty("user.home") + File.separator + SystemManagerDialog.MessageID.Unnamed.toString());
	}

	@Override
	public REpiceaGUIPermission getGUIPermission() {return readWrite;}

	
	@SuppressWarnings("rawtypes")
	@Override
	public MemorizerPackage getMemorizerPackage() {
		MemorizerPackage mp = new MemorizerPackage();
		mp.add(filename);
		mp.add((ArrayList) processors);
		return mp;
	}


	@SuppressWarnings({ "unchecked", "rawtypes" })
	@Override
	public void unpackMemorizerPackage(MemorizerPackage wasMemorized) {
		String filename = (String) wasMemorized.remove(0);
		setFilename(filename);
		List<Processor> retrievedProcessors = (ArrayList) wasMemorized.remove(0);
		processors.clear();
		processors.addAll(retrievedProcessors);
	}

	@Override
	public boolean isVisible() {
		return guiInterface != null && guiInterface.isVisible();
	}
	
	@SuppressWarnings({ "rawtypes"})
	protected void checkForEndlessLoops() {
		for (Processor processor : getList()) {
			processor.setPartOfEndlessLoop(false);
		}
//		List<ProcessUnit> inputUnits = new ArrayList<ProcessUnit>();
//		addTestUnits(inputUnits);
////		inputUnits.add(new TestProcessUnit());
//		
//		for (Processor processor : getPrimaryProcessors()) {
//			processor.doProcess(inputUnits);
//		}
		List<ProcessUnit> inputUnits;
		
		for (Processor processor : getPrimaryProcessors()) {
			inputUnits = new ArrayList<ProcessUnit>();
			addTestUnits(inputUnits);
			processor.doProcess(inputUnits);
		}

	}
	
	protected void addTestUnits(List<ProcessUnit> inputUnits) {
		inputUnits.add(new TestProcessUnit());
	}
	
	@SuppressWarnings("rawtypes")
	public static void main(String[] args) {
		Processor unit1 = new Processor("1");
		Processor unit2 = new Processor("2");
		Processor unit3 = new Processor("3");
		Processor unit4 = new Processor("4");
		unit1.addSubProcessor(unit2);
		unit1.addSubProcessor(unit3);
		unit2.addSubProcessor(unit3);
		unit2.addSubProcessor(unit4);
		SystemManager manager = new SystemManager();
		manager.registerObject(unit1);
		manager.registerObject(unit2);
		manager.registerObject(unit3);
		manager.registerObject(unit4);

		List<ProcessUnit> units = new ArrayList<ProcessUnit>();
		units.add(new TestProcessUnit());
		unit1.doProcess(units);
		@SuppressWarnings("unused")
		List<Processor> processors = manager.getList();
		manager.showUI(null);
		System.exit(0);
	}

}
