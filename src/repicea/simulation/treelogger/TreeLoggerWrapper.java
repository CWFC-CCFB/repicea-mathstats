/*
 * This file is part of the repicea library.
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

import java.io.File;
import java.io.IOException;
import java.util.List;
import java.util.concurrent.CopyOnWriteArrayList;

import javax.swing.filechooser.FileFilter;

import repicea.io.Workspaceable;
import repicea.simulation.Parameterizable;

/**
 * The TreeLoggerWrapper class contains a TreeLogger instance that can be changed. Every time a change occurs, the class
 * notifies the listeners.
 * @author Mathieu Fortin - January 2013
 */
@Deprecated
public final class TreeLoggerWrapper implements Parameterizable, Workspaceable {

	private List<TreeLoggerChangeListener> listeners;
	
	private TreeLogger<?, ?> treeLogger;
	
	private String workspace;
	
	public TreeLoggerWrapper() {
		listeners = new CopyOnWriteArrayList<TreeLoggerChangeListener>();
	}
	
	/**
	 * This method sets the TreeLogger instance that this wrapper contains.
	 * @param treeLogger a TreeLogger instance
	 */
	public void setTreeLogger(TreeLogger<?, ?> treeLogger) {
		if (this.treeLogger != null) {
			this.treeLogger.wrapper = null;
		}
		this.treeLogger = treeLogger;
		treeLogger.wrapper = this;
		fireTreeLoggerEvent(treeLogger);
	}
	
	/**
	 * This method returns the TreeLogger instance that is contained in this wrapper.
	 * @return a TreeLogger instance
	 */
	public TreeLogger<?,?> getTreeLogger() {
		return treeLogger;
	}
	
	protected void fireTreeLoggerEvent(TreeLogger<?, ?> treeLogger) {
		TreeLoggerEvent evt = new TreeLoggerEvent(treeLogger);
		for (TreeLoggerChangeListener listener : listeners) {
			listener.treeLoggerChanged(evt);
		}
	}

	/**
	 * This method adds a TreeLoggerChangeListener instance to the listeners only and only
	 * if the listeners do not already contain the instance.
	 * @param listener a TreeLoggerChangeListener instance 
	 */
	public void addTreeLoggerChangeListener(TreeLoggerChangeListener listener) {
		if (!listeners.contains(listener)) {
			listeners.add(listener);
			if (treeLogger != null) {
				TreeLoggerEvent evt = new TreeLoggerEvent(treeLogger);
				listener.treeLoggerChanged(evt);						// update the listener
			}
		} 
	}
	
	
	/**
	 * This method removes a TreeLoggerChangeListener instance to the listeners only and only
	 * if the listeners already contain the instance.
	 * @param listener a TreeLoggerChangeListener instance 
	 */
	public void removeTreeLoggerChangeListener(TreeLoggerChangeListener listener) {
		while (listeners.contains(listener)) {
			listeners.remove(listener);
		}
	}

	@Override
	public void loadFromFile(String filename) throws IOException {
		TreeLoggerParameters<?> parameters = TreeLoggerParameters.loadFromFile(filename);
		TreeLogger<?,?> treeLogger = parameters.createTreeLoggerInstance();
		setTreeLogger(treeLogger);
	}

	@Override
	public String getWorkspace() {
		return workspace;
	}

	@Override
	public void setWorkspace(String workspace) {
		this.workspace = workspace;
	}

	@Override
	public FileFilter getFileFilter() {return TreeLoggerParameters.TreeLoggerFileFilter;}

	@Override
	public String getFilename() {
		if (treeLogger != null) {
			return treeLogger.getTreeLoggerParameters().getFilename();
		} else {
			return "";
		}
	}

	@Override
	public String getName() {
		return getFilename().replace(workspace + File.separator, "");
	}

	
}
