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
package repicea.simulation;

import java.io.IOException;

import javax.swing.filechooser.FileFilter;

/**
 * This interface ensures the class can be parameterized through a file.
 * @author Mathieu Fortin - January 2013
 */
@Deprecated
public interface Parameterizable {

	/**
	 * This method reads the parameters of the instance in a file.
	 * @param filename the path of the file 
	 * @throws IOException
	 */
	public void loadFromFile(String filename) throws IOException;
	
	
	/**
	 * This method returns the appropriate FileFilter instance for this class.
	 * @return a FileFilter instance
	 */
	public FileFilter getFileFilter();

	/**
	 * This method returns the filename of the parameters.
	 * @return a String
	 */
	public String getFilename();
	
	
	/**
	 * This method returns a name that defines the parameters. Typically, it can be the filename without
	 * the workspace.
	 * @return a String
	 */
	public String getName();
	
	
}
