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

import org.junit.Assert;
import org.junit.Test;

import repicea.math.Matrix;
import repicea.stats.StatisticalUtility;
import repicea.stats.estimates.MonteCarloEstimate;


public class TruncatedGaussianUtilityTest {

	@Test
	public void simpleProbabilityDensity() {
		double value = TruncatedGaussianUtility.getProbabilityDensity(0, 0, 1, Double.NEGATIVE_INFINITY, 0);
		Assert.assertEquals("Testing probability density", 0.7978845608028654, value, 1E-8);
	}
	
	@Test
	public void simpleCumulativeProbability() {
		double value = TruncatedGaussianUtility.getCumulativeProbability(0, 0, 1, Double.NEGATIVE_INFINITY, 0);
		Assert.assertEquals("Testing cumulative probability", 1d, value, 1E-8);
	}

	@Test
	public void simpleQuantile() {
		double value = TruncatedGaussianUtility.getQuantile(0.99, 0, 1, Double.NEGATIVE_INFINITY, 0);
		Assert.assertEquals("Testing cumulative probability", -0.012533469522069105, value, 1E-8);
	}

	@Test
	public void simpleQuantile2() {
		double value = TruncatedGaussianUtility.getQuantile(0.01, 0, 1, Double.NEGATIVE_INFINITY, 0);
		Assert.assertEquals("Testing cumulative probability", -2.5758293064439264, value, 1E-8);
	}


	@Test
	public void simpleQuantile3() {
		double value = TruncatedGaussianUtility.getQuantile(0.005, 0, 1, Double.NEGATIVE_INFINITY, Double.POSITIVE_INFINITY);
		Assert.assertEquals("Testing cumulative probability", -2.5758293064439264, value, 1E-8);
	}

	@Test
	public void randomGeneration() {
		double mu = 6;
		double var = 10;
		double std = Math.sqrt(var);
		double a = 2;
		double b = Double.POSITIVE_INFINITY;
		double alpha = (2 - mu) / std;
		double Z = 1d - GaussianUtility.getCumulativeProbability(alpha);
		
		MonteCarloEstimate estimate = new MonteCarloEstimate();
		for (int i = 0; i < 1000000; i++)
			estimate.addRealization(new Matrix(1,1, StatisticalUtility.getRandom().nextTruncatedGaussian(mu, var, a, b), 0));
		
		double expectedMean = GaussianUtility.getProbabilityDensity(alpha) * std / Z + mu;
		double observedMean = estimate.getMean().getValueAt(0, 0);
		Assert.assertEquals("Comparing means", expectedMean, observedMean, 1E-2);
		double expectedVariance = var * (1 + alpha * GaussianUtility.getProbabilityDensity(alpha)/Z - (GaussianUtility.getProbabilityDensity(alpha)/Z) * (GaussianUtility.getProbabilityDensity(alpha)/Z));
		double observedVariance = estimate.getVariance().getValueAt(0, 0);
		Assert.assertEquals("Comparing variances", expectedVariance, observedVariance, 3E-2);
	}
	
	@Test
	public void checkingDensityCalculation() {
		double mu = 6;
		double var = 10;
		double std = Math.sqrt(var);
		double a = 2;
		double b = Double.POSITIVE_INFINITY;
	
		double sum = 0;
		for (double i = a; i < 10000; i += 0.1) {
			double x1 = TruncatedGaussianUtility.getProbabilityDensity(i, mu, var, a, b);
			double x2 = TruncatedGaussianUtility.getProbabilityDensity(i + .1, mu, var, a, b);
			sum += (x1 + x2) * .5 * .1;
		}
		Assert.assertEquals("Comparing cumulative probabilities", 1d, sum, 1E-3);
	}

	@Test
	public void checkingDensityCalculation2() {
		double mu = 4;
		double var = 12;
		double std = Math.sqrt(var);
		double a = 2;
		double b = Double.POSITIVE_INFINITY;
	
		double sum = 0;
		for (double i = a; i < 10000; i += 0.1) {
			double x1 = TruncatedGaussianUtility.getProbabilityDensity(i, mu, var, a, b);
			double x2 = TruncatedGaussianUtility.getProbabilityDensity(i + .1, mu, var, a, b);
			sum += (x1 + x2) * .5 * .1;
		}
		Assert.assertEquals("Comparing cumulative probabilities", 1d, sum, 1E-3);
	}


}
