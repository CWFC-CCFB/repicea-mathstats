/*
 * This file is part of the repicea library.
 *
 * Copyright (C) 2009-2019 Mathieu Fortin for Rouge-Epicea
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
package repicea.simulation.covariateproviders.plotlevel;

import java.util.List;

import repicea.simulation.covariateproviders.StochasticImplementation;

/**
 * The StochasticInformationProvider interface provides basic information about the simulation mode.
 * @author Mathieu Fortin - November 2014
 */
public interface StochasticInformationProvider<Realization> extends StochasticImplementation {

	/**
	 * This method returns the number of realizations to be made.
	 * @return an integer
	 */
	public default int getNumberRealizations() {
		return getRealizationIds().size();
	}

	/**
	 * Returns the ids of the different realizations
	 * @return a List of Integer
	 */
	public List<Integer> getRealizationIds();
	
	
	
	/**
	 * This method returns the realization.
	 * @param realizationID the id of the realization
	 * @return an instance whose class is defined by the Realization parameter
	 */
	public Realization getRealization(int realizationID);
	
	
}
