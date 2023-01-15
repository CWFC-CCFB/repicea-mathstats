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
package repicea.simulation;

/**
 * This interface defines the method to access growth model predictions. 
 * @author Mathieu Fortin - November 2012
 * @param <S> a class that represents the stand
 * @param <T> a class that represents the tree
 */
public interface GrowthModel<S, T> {

	/**
	 * This method returns the growth of a particular tree represented by the parameter
	 * "tree".
	 * @param stand a S-derived instance
	 * @param tree a T-derived instance
 	 * @param parms some additional parameters
	 * @return the growth
	 */
	public double predictGrowth(S stand, T tree, Object... parms);
	
}
