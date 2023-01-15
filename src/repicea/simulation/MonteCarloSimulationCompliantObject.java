/*
 * This file is part of the repicea  library.
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


public interface MonteCarloSimulationCompliantObject {
	
	/**
	 * This method returns an object that makes it possible to identify
	 * the subject that implements this interface. This id remains constant 
	 * throughout the Monte Carlo iterations in case of stochastic implementation.
	 * @return a String that defines the subject id and that remains constant throughout the simulation
	 */
	public String getSubjectId();

	/**
	 * This method returns the hierarchical levels of the object.
	 * @return a HierarchicalLevel instance
	 */
	public HierarchicalLevel getHierarchicalLevel();
	
	
//	/**
//	 * This method sets the MonteCarlo id of the subject. Some instances might have a different implementation of this and might not have to use this method.
//	 * @param i the MonteCarlo id
//	 */
//	public void setMonteCarloRealizationId(int i);

	
	/**
	 * This method returns the id of the Monte Carlo realization. It is necessary for the implementation 
	 * of the random deviates on the parameter estimates. These deviates remain constant for a particular
	 * Monte Carlo iteration, regardless of the plot.
	 * @return an integer
	 */
	public int getMonteCarloRealizationId();	
	
}
