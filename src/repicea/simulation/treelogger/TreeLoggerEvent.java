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

import java.util.EventObject;

/**
 * The TreeLoggerEvent class is a basic event that is sent every time the TreeLogger instance
 * is changed in the TreeLoggerWrapper class.
 * @author Mathieu Fortin - January 2013
 */
@SuppressWarnings("serial")
public final class TreeLoggerEvent extends EventObject {

	/**
	 * Protected constructor since it is generated only in this package.
	 * @param source a TreeLogger instance that generated this event
	 */
	protected TreeLoggerEvent(TreeLogger<?,?> source) {
		super(source);
	}


	
	
}
