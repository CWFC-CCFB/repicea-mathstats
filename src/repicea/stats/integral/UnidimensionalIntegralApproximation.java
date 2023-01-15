/*
 * This file is part of the repicea library.
 *
 * Copyright (C) 2009-2022 Mathieu Fortin for Rouge-Epicea
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
package repicea.stats.integral;

import repicea.math.EvaluableFunction;

/**
 * Provides a simple method for numerical integration of unidimensional integrals.
 * @author Mathieu Fortin - July 2022
 */
public interface UnidimensionalIntegralApproximation<P extends EvaluableFunction<Double>> {

	/**
	 * Compute the numerical integration for one-dimension integrals.
	 * @param functionToEvaluate an AbstractMathematicalFunction function
	 * @param index the index of the parameter or the variable that is integrated
	 * @param isParameter a boolean true it is a parameter, false it is a variable
	 * @return a double
	 */
	public double getIntegralApproximation(P functionToEvaluate, 
			int index, 
			boolean isParameter);

	

}
