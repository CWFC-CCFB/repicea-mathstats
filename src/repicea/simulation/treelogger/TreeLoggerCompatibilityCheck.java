/*
 * This file is part of the repicea library.
 *
 * Copyright (C) 2009-2018 Mathieu Fortin for Rouge-Epicea
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

/**
 * The TreeLoggerCompatibilityCheck class provides some parameters than can be included in the compatibility check
 * of the tree loggers.
 * @author Mathieu Fortin - Sept 2018
 */
public class TreeLoggerCompatibilityCheck {
	
	private final Object treeInstance;
//	private final List<String> speciesNameList;
	
	public TreeLoggerCompatibilityCheck(Object treeInstance) {
		this.treeInstance = treeInstance;
//		this.speciesNameList = new ArrayList<String>();
//		this.speciesNameList.addAll(speciesNameList);
	}

	public Object getTreeInstance() {return treeInstance;}

//	public List<String> getSpeciesNameList() {
//		List<String> copyList = new ArrayList<String>();
//		copyList.addAll(speciesNameList);
//		return copyList;
//	}


}
