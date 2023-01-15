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

public class TruncatedGaussianUtility {

	/**
	 * Compute the probability density of a truncated normal distribution. 
	 * 
	 * @param x the value for which the density is to be calculated
	 * @param mu the mean of the original normal distribution
	 * @param sigma2 the variance of the original distribution
	 * @param a the lower bound
	 * @param b the upper bound
	 * @return the probability density
	 */
	public static double getProbabilityDensity(double x, double mu, double sigma2, double a, double b) {
		checkValue(sigma2, a, b);
		if (x < a || x > b) {
			return 0d;
		}
		double sigma = Math.sqrt(sigma2);
		double alpha = (a - mu) / sigma;
		double beta = (b - mu) / sigma;
		double epsilon = (x - mu) / sigma;
		double Z = GaussianUtility.getCumulativeProbability(beta) - GaussianUtility.getCumulativeProbability(alpha);
		double density = GaussianUtility.getProbabilityDensity(epsilon) / (sigma * Z);
		return density;
	}
	
	private static void checkValue(double sigma2, double a, double b) {
		if (b <= a) {
			throw new InvalidParameterException("The b argument must be larger than a!");
		}
		if (sigma2 <= 0) {
			throw new InvalidParameterException("The sigma2 argument must be larger than 0!");
		}
	}
	
	/**
	 * Compute the cumulative probability of a truncated normal distribution. 
	 * 
	 * @param x the value for which the cumulative density is to be calculated, i.e. F(X<x)
	 * @param mu the mean of the original normal distribution
	 * @param sigma2 the variance of the original distribution
	 * @param a the lower bound
	 * @param b the upper bound
	 * @return the cumulative probability 
	 */
	public static double getCumulativeProbability(double x, double mu, double sigma2, double a, double b) {
		checkValue(sigma2, a, b);
		if (x < a) {
			return 0;
		} else if (x > b) {
			return 1;
		} else {
			double sigma = Math.sqrt(sigma2);
			double alpha = (a - mu) / sigma;
			double beta = (b - mu) / sigma;
			double epsilon = (x - mu) / sigma;
			double Z = GaussianUtility.getCumulativeProbability(beta) - GaussianUtility.getCumulativeProbability(alpha);
			double cumulativeProbability = (GaussianUtility.getCumulativeProbability(epsilon) - GaussianUtility.getCumulativeProbability(alpha)) / Z;
			return cumulativeProbability;
		}
	}


	/**
	 * Compute the quantile of a truncated normal distribution corresponding to a given
	 * cumulative probability. 
	 * 
	 * @param x the cumulative probability
	 * @param mu the mean of the original normal distribution
	 * @param sigma2 the variance of the original distribution
	 * @param a the lower bound
	 * @param b the upper bound
	 * @return the quantile
	 */
	public static double getQuantile(double x, double mu, double sigma2, double a, double b) {
		checkValue(sigma2, a, b);
		if (x <= 0 || x >= 1) {
			throw new InvalidParameterException("The x argument must take a value between 0 and 1!");
		}
		double sigma = Math.sqrt(sigma2);
		double alpha = (a - mu) / sigma;
		double beta = (b - mu) / sigma;
		double Z = GaussianUtility.getCumulativeProbability(beta) - GaussianUtility.getCumulativeProbability(alpha);
		return GaussianUtility.getQuantile(x * Z + GaussianUtility.getCumulativeProbability(alpha)) * sigma + mu;
	}

}
