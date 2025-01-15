/*
 * This file is part of the repicea-mathstats library.
 *
 * Copyright (C) 2025 His Majesty the King in right of Canada
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
package repicea.stats.sampling;

import repicea.math.Matrix;
import repicea.math.SymmetricMatrix;

/**
 * An interface to define finite population point estimates.
 * 
 * @author Mathieu Fortin - January 2025
 */
public interface FinitePopulationPointEstimate extends PointEstimate {

	/**
	 * Provide the size of the population.<p>
	 * The size is expressed in terms of number of population units
	 * under the assumption they are all the same size.
	 * @return a double
	 */
	public double getPopulationSize();
	
	/**
	 * Provide the estimated total of the population.
	 * @return a Matrix instance
	 */
	public default Matrix getTotal() {
		return getMean().scalarMultiply(getPopulationSize());
	}

	/**
	 * Provide the estimated variance for the total of the population.
	 * @return a SymmetricMatrix instance
	 */
	public default SymmetricMatrix getVarianceOfTotal() {
		return getVariance().scalarMultiply(getPopulationSize() * getPopulationSize());
	}
}
