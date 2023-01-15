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
package repicea.stats.distributions.utility;

import java.security.InvalidParameterException;

import repicea.stats.StatisticalUtility;

public class WeibullUtility {

	
	/**
	 * Compute the probability density of a Weibull distribution. <br>
	 * <br>
	 * The Weibull parameterization is pdf(y) = k/lambda * ((y - theta)/lambda)^(k-1) * e^(-((y - theta)/lambda)^k).
	 * @param y the value
	 * @param k the shape parameter (must be greater than 0)
	 * @param lambda the scale parameter (must be greater than 0)
	 * @param theta the location parameter
	 * @return a probability density
	 */
	public static double getProbabilityDensity(double y, double k, double lambda, double theta) {
		WeibullUtility.checkShapeAndScale(k, lambda);
		if (y - theta < 0) {
			return 0d;
	 	} else {
			return k / lambda * Math.pow((y - theta)/lambda, k - 1) * Math.exp(-Math.pow((y - theta)/lambda, k));
	 	}
	}
	
	/**
	 * Compute the probability density of a Weibull distribution with location parameter theta = 0. <br>
	 * <br>
	 * The Weibull parameterization is pdf(y) = k/lambda * (y/lambda)^(k-1) * e^(-(y/lambda)^k).
	 * @param y the value
	 * @param k the shape parameter (must be greater than 0)
	 * @param lambda the scale parameter (must be greater than 0)
	 * @return a probability density
	 */
	public static double getProbabilityDensity(double y, double k, double lambda) {
		return getProbabilityDensity(y, k, lambda, 0);  // theta is set to 0
	}

	
	private static void checkShapeAndScale(double k, double lambda) {
		if (k <= 0d || lambda <= 0d) {
			throw new InvalidParameterException("The k and lambda argument must be strictly positive (i.e. > 0)!");
		}
	}

	
	/**
	 * Compute the cumulative probability of a Weibull distribution.<br>
	 * <br>
	 * The Weibull parameterization is cdf(y) = 1 - e^(-((y-theta)/lambda)^k).
	 * @param y the value
	 * @param k the shape parameter (must be greater than 0)
	 * @param lambda the scale parameter (must be greater than 0)
	 * @param theta the location parameter
	 * @return a cumulative probability
	 */
	public static double getCumulativeProbability(double y, double k, double lambda, double theta) {
		WeibullUtility.checkShapeAndScale(k, lambda);
		if (y - theta < 0) {
			return 0d;
	 	} else {
	 		return 1 - Math.exp(-Math.pow((y - theta)/lambda, k));
	 	}
	}

	/**
	 * Compute the cumulative probability of a Weibull distribution with location parameter theta = 0.<br>
	 * <br>
	 * The Weibull parameterization is cdf(y) = 1 - e^(-(y/lambda)^k).
	 * @param y the value
	 * @param k the shape parameter (must be greater than 0)
	 * @param lambda the scale parameter (must be greater than 0)
	 * @return a cumulative probability
	 */
	public static double getCumulativeProbability(double y, double k, double lambda) {
		return WeibullUtility.getCumulativeProbability(y, k, lambda, 0); // location parameter set to 0
	}

	/**
	 * Return the quantile of a Weibull distribution. <br>
	 * <br>
	 * The Weibull parameterization is cdf(y) = 1 - e^(-(y/lambda)^k).
	 * @param q the quantile (between 0 and 1)
	 * @param k the shape parameter (must be greater than 0)
	 * @param lambda the scale parameter (must be greater than 0)
	 * @param theta the location parameter
	 * @return the quantile of the distribution
	 */
	public static double getQuantile(double q, double k, double lambda, double theta) {
		WeibullUtility.checkShapeAndScale(k, lambda);
		if (q < 0 || q > 1) {
			throw new InvalidParameterException("The q argument must range between 0 and 1!");
		}
		return lambda * Math.pow(-Math.log(1 - q), 1/k) + theta;
	}
	
	/**
	 * Return the quantile of a Weibull distribution. <br>
	 * <br>
	 * The Weibull parameterization is cdf(y) = 1 - e^(-(y/lambda)^k).
	 * @param q the quantile (between 0 and 1)
	 * @param k the shape parameter (must be greater than 0)
	 * @param lambda the scale parameter (must be greater than 0)
	 * @return the quantile of the distribution
	 */
	public static double getQuantile(double q, double k, double lambda) {
		return getQuantile(q, k, lambda, 0d);
	}


}
