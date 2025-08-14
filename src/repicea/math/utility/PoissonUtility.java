/*
 * This file is part of the repicea-mathstats library.
 *
 * Copyright (C) 2025 His Majesty the King in Right of Canada
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
package repicea.math.utility;

import java.security.InvalidParameterException;

/**
 * An implementation of the Poisson distribution.
 * @author Mathieu Fortin - August 2025
 */
public class PoissonUtility {

	static final double SMALL_VALUE = 1E-8;

	/**
	 * Provide the probability mass of y given a Poisson distribution with parameter 
	 * lambda.
	 * @param y a positive integer (must be greater than or equal to 0)
	 * @param lambda the mean of the Poisson distribution
	 * @return the probability mass of observing y
	 */
	public static double getProbabilityMass(int y, double lambda) {
		if (y < 0) {
			throw new InvalidParameterException("The parameter y must be an integer greater to or equal to 0!");
		}
		return Math.pow(lambda, y) * Math.exp(-lambda) / MathUtility.Factorial(y);
	}

	/**
	 * Provide the probability mass of y given a Poisson distribution with parameter 
	 * lambda.
	 * @param y a positive integer (must be greater than or equal to 0)
	 * @param lambda the mean of the Poisson distribution
	 * @return the probability mass of observing y or smaller values
	 */
	public static double getCumulativeProbabilityMass(int y, double lambda) {
		if (y < 0) {
			throw new InvalidParameterException("The parameter y must be an integer greater to or equal to 0!");
		}
		double cumulativeMass = 0d;
		for (int x = 0; x <= y; x++) {
			double probabilityMass = getProbabilityMass(x, lambda);
			cumulativeMass += probabilityMass;
			if (x > lambda && probabilityMass < SMALL_VALUE) {
				break;
			} 
		}
		return cumulativeMass;
	}
	
	/**
	 * Provide a quantile of the distribution.
	 * @param cmfValue the cumulative mass (must greater than 0 and smaller than 1)
	 * @param lambda the mean of the distribution
	 * @return a quantile
	 */
	public static int getQuantile(double cmfValue, double lambda) {
		if (cmfValue <= 0 || cmfValue >= 1) 
			throw new InvalidParameterException("The cmfValue parameter should be a double between 0 and 1!");
		double cumulativeMass = 0d;
		int y = 0;
		while(cumulativeMass < cmfValue) 
			cumulativeMass += getProbabilityMass(y++, lambda); 
		return --y;
	}

}
