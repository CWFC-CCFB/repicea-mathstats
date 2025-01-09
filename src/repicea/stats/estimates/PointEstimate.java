/*
 * This file is part of the repicea-mathstats library.
 *
 * Copyright (C) 2009-2018 Mathieu Fortin for Rouge-Epicea
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
package repicea.stats.estimates;

import repicea.math.Matrix;
import repicea.math.SymmetricMatrix;
import repicea.stats.distributions.GaussianDistribution;

/**
 * An abstract class for point estimates (total or mean).<p>
 * This class is more focused on sampling. 
 * @author Mathieu Fortin - March 2021, January 2025
 */
public interface PointEstimate extends Estimate<Matrix, SymmetricMatrix, GaussianDistribution> {

	/**
	 * Check if the population size is available.
	 * @return a boolean
	 */
	public default boolean isPopulationSizeKnown() {
		return getPopulationSize() != -1;
	}

	/**
	 * Provide the population size in terms of number of population units.
	 * @return a double
	 */
	public double getPopulationSize();

}