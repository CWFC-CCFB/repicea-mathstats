/*
 * This file is part of the repicea-simulation library.
 *
 * Copyright (C) 2009-2014 Mathieu Fortin for Rouge-Epicea
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

import repicea.gui.REpiceaAWTProperty;

/**
 * The TreeLoggerAWTProperty class contains the event that are thrown by the 
 * TreeLoggerParameterDialog class.
 * @author Mathieu Fortin - May 2014
 */
public class TreeLoggerAWTProperty extends REpiceaAWTProperty {

	public static final TreeLoggerAWTProperty LogGradeRemoved = new TreeLoggerAWTProperty("LogGradeRemoved");
	public static final TreeLoggerAWTProperty LogGradeAdded = new TreeLoggerAWTProperty("LogGradeAdded");
	public static final TreeLoggerAWTProperty SpeciesRemoved = new TreeLoggerAWTProperty("SpeciesRemoved");
	public static final TreeLoggerAWTProperty SpeciesAdded = new TreeLoggerAWTProperty("SpeciesAdded");
	
	
	
	protected TreeLoggerAWTProperty(String propertyName) {
		super(propertyName);
	}

}
