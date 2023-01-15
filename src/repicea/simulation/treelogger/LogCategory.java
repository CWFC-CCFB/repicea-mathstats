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

import java.io.Serializable;
import java.util.List;

import repicea.gui.REpiceaUIObject;

/**
 * This abstract class ensures that any TreeLogCategory-derived object has
 * at least (1) a name, (2) a list of end use products and (3) a method that
 * returns the yield of the piece.
 * @author Mathieu Fortin - January 2010
 */
public abstract class LogCategory implements Serializable, REpiceaUIObject {

	private static final long serialVersionUID = 20100804L;
	
	private String name;
	private Object species;
	private boolean isFromStump = false;

	/**
	 * General constructor for all the TreeLogCategory-derived classes.
	 * @param name the name of the log category
	 * @param isFromStump true if this log category is extracted from the stump and the roots
	 */
	public LogCategory(String name, boolean isFromStump) {
		this.name = name;
		this.isFromStump = isFromStump;
	}

	/**
	 * Default instantiation.
	 */
	public LogCategory() {
		this("unnamed", false);
	}

	protected void setSpecies(Object species) {
		this.species = species;
	}
	
	/**
	 * This method returns the name of this tree log category.
	 * @return a String
	 */
	public String getName() {return name;}
	
	/**
	 * This method sets the name of this tree log category.
	 */
	public void setName(String str) {name = str;}

	public Object getSpecies() {return species;}
	
	@Override
	public String toString() {return getName();}
	
	/**
	 * This method returns the yield of this piece, i.e. the volume of end use products with respect
	 * to the initial volume of this piece.
	 * @param piece a WoodPiece instance
	 * @return a Double
	 * @throws Exception
	 */
	public abstract double getYieldFromThisPiece(WoodPiece piece) throws Exception;

	
	@Override
	public boolean equals(Object obj) {
		if (obj instanceof LogCategory) {
			LogCategory tlc = (LogCategory) obj;
			if (tlc.getName().equals(getName())) {
				if (tlc.getSpecies() == null && getSpecies() == null) {
					return true;
				} else if (tlc.getSpecies().equals(getSpecies())) {
					return true;
				}
			}
		} 
		return false;
	}
	
	/**
	 * This method returns true if the log category is actually extracted from stumps. Typically, this is
	 * done with maritime pine, where stumps are converted into energy wood.
	 * @return a boolean
	 */
	public boolean isFromStump() {return isFromStump;}
	

	/**
	 * This method returns a wood piece if it can be extract from the tree.
	 * @param tree a LoggableTree instance
	 * @param parms optional parameters
	 * @return a WoodPiece instance or null if cannot be extracted from the tree
	 */
	protected abstract List<? extends WoodPiece> extractFromTree(LoggableTree tree, Object... parms);
	
	
	
}
