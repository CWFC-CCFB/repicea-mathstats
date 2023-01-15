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
package repicea.simulation.climate;

import repicea.simulation.MonteCarloSimulationCompliantObject;
import repicea.simulation.covariateproviders.plotlevel.GeographicalCoordinatesProvider;

/**
 * The REpiceaClimateChangeGenerator interface ensures that the climate generator also predicts 
 * the climate change and as such it can produce annual change for different climate variables.
 * @author Mathieu Fortin - June 2016
 */
public interface REpiceaClimateChangeGenerator<P extends GeographicalCoordinatesProvider> extends REpiceaClimateGenerator<P> {

	
	
	
	/**
	 * Returns Map of the annual changes for the different climate variable for a MonteCarloSimulationCompliantObject
	 * instance, typically a plot.
	 * @param obj a MonteCarloSimulationCompliantObject instance
	 * @return a Map of ClimateVariable (key) and Double (value)
	 */
	public REpiceaClimateChangeTrend getClimateTrendForThisStand(MonteCarloSimulationCompliantObject obj);
	
	
}
