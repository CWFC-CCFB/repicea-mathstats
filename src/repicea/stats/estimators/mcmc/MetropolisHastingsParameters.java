/*
 * This file is part of the repicea-mathstats library.
 *
 * Copyright (C) 2021-24 His Majesty the King in Right of Canada
 * Author: Mathieu Fortin, Canadian Forest Service
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
package repicea.stats.estimators.mcmc;

/**
 * A class containing the simulation parameters for the MCMC algorithm.
 * @author Mathieu Fortin - September 2021
 */
public class MetropolisHastingsParameters {

	/**
	 * The number of samples in the burn-in period.<p>
	 * Is normally set to 10000.
	 */
	public int nbBurnIn = 10000;
	
	/**
	 * Total number of accepted realizations in the Markov Chain.<p>
	 * Is normally set to 500000 + {@link MetropolisHastingsParameters#nbBurnIn}. 
	 */
	public int nbAcceptedRealizations = 500000 + nbBurnIn;
	
	/**
	 * The number of internal realizations.<p>
	 * The algorithm will stop if no realization has been accepted in the chain after this number of realizations. 
	 * Is normally set to 100000.
	 */
	public int nbInternalIter = 100000;
	
	/**
	 * The selection rate to obtain the final sample.<p>
	 * Once the chain contains the appropriate number of realizations (see {@link MetropolisHastingsParameters#nbAcceptedRealizations}),
	 * then those of the burn-in period are discarded and one realization every x realizations is selected to the final sample. This rate 
	 * is set through this parameter. Is normally set to 50.
	 */
	public int oneEach = 50;
	
	/**
	 * Set the size of the initial grid.<p>
	 * The realizations are generated at random and that with the greatest likelihood is selected as
	 * starting value for the Markov Chain.
	 */
	public int nbInitialGrid = 10000;	

	/**
	 * An empty constructor with default values.
	 */
	public MetropolisHastingsParameters() {}

	@Override
	public MetropolisHastingsParameters clone() {
		try {
			return (MetropolisHastingsParameters) super.clone();
		} catch (CloneNotSupportedException e) {
			e.printStackTrace();
			return null;
		}
	}

}
