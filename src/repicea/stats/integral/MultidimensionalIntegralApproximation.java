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

import java.util.List;

import repicea.math.EvaluableFunction;

/**
 * Provides a simple method for numerical integration of multidimensional integrals.
 * @author Mathieu Fortin - July 2022
 */
public interface MultidimensionalIntegralApproximation<P extends EvaluableFunction<Double>> {

	
	/**
	 * This method returns the value of a multi-dimension integral
	 * @param functionToEvaluate a MathematicalFunction instance that returns Double 
	 * @param indices the indices of the parameters over which the integration is made
	 * @param isParameter a boolean to indicate that indices refer to parameters. If false, it is assumed that the
	 * indices refer to variables.
	 * @param lowerCholeskyTriangle the lower triangle of the Cholesky factorization of the variance-covariance matrix
	 * @return the approximation of the integral
	 */
	public double getMultiDimensionalIntegralApproximation(P functionToEvaluate,
			List<Integer> indices, 
			boolean isParameter);
}
