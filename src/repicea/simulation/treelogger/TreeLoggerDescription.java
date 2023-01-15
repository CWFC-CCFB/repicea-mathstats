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

import java.security.InvalidParameterException;

/**
 * The TreeLoggerDescription class has an address for a TreeLogger instance. 
 * @author Mathieu Fortin - November 2012
 */
public class TreeLoggerDescription {
	
	private final String treeLoggerClassName;
	
	/**
	 * Constructor.
	 * @param className a String that represents the class of the tree logger.
	 */
	@SuppressWarnings({ "unchecked", "unused" })
	public TreeLoggerDescription(String className) {
		try {
			Class<? extends TreeLogger<?, ?>> clazz = (Class<? extends TreeLogger<?,?>>) Class.forName(className);
//			Class<? extends TreeLogger<?, ?>> clazz = (Class<? extends TreeLogger<?,?>>) ClassLoader.getSystemClassLoader().loadClass(className);
			this.treeLoggerClassName = className;
		} catch (Exception e) {
			throw new InvalidParameterException("Class name " + className + " is not a valid TreeLogger class");
		}
	}
	
	/**
	 * Constructor 2.
	 * @param clazz a TreeLogger class
	 */
	public TreeLoggerDescription(Class<? extends TreeLogger<?, ?>> clazz) {
		this.treeLoggerClassName = clazz.getName();
	}
	
	
	/**
	 * This method returns the TreeLogger class.
	 * @return a Class instance
	 */
	@SuppressWarnings("unchecked")
	public Class<? extends TreeLogger<?,?>> getTreeLoggerClass() {
		try {
//			return (Class<? extends TreeLogger<?,?>>) ClassLoader.getSystemClassLoader().loadClass(treeLoggerClassName);
			return (Class<? extends TreeLogger<?,?>>) Class.forName(treeLoggerClassName);
		} catch (ClassNotFoundException e) {
			throw new InvalidParameterException("Class name " + treeLoggerClassName + " is not a valid TreeLogger class");
		}
	}
	

	@Override
	public String toString() {
		return getTreeLoggerClass().getSimpleName();
	}

	/**
	 * This method instantiate the tree logger. If the scriptMode parameter is set to true, the TreeLogger instance is set
	 * to its default parameters. Otherwise it is unparameterized.
	 * @param scriptMode a boolean that takes the value true if in scriptMode
	 * @return a TreeLogger instance
	 */
	@SuppressWarnings({ "rawtypes", "unchecked" })
	public TreeLogger instantiateTreeLogger(boolean scriptMode) {
		TreeLogger treeLogger;
		try {
			treeLogger = (TreeLogger) getTreeLoggerClass().newInstance();
			if (scriptMode) {
				TreeLoggerParameters params = treeLogger.createDefaultTreeLoggerParameters();
				treeLogger.setTreeLoggerParameters(params);
			}
			return treeLogger;
		} catch (InstantiationException e) {
			e.printStackTrace();
		} catch (IllegalAccessException e) {
			e.printStackTrace();
		}
		return null;
	}
	
	@Override
	public boolean equals(Object obj) {
		if (obj instanceof TreeLoggerDescription) {
			TreeLoggerDescription tld = (TreeLoggerDescription) obj;
			if (tld.treeLoggerClassName.equals(treeLoggerClassName)) {
				return true;
			}
		} 
		return false;
	}
	
	
}

