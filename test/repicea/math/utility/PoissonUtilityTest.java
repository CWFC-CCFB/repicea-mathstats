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

import org.junit.Assert;
import org.junit.Test;

import repicea.math.Matrix;
import repicea.stats.StatisticalUtility;
import repicea.stats.estimates.MonteCarloEstimate;

public class PoissonUtilityTest {

	@Test
	public void test01SimpleValue() {
		double observed = PoissonUtility.getProbabilityMass(2, 1);
		Assert.assertEquals("Testing expected probability mass", 0.18393972058572117, observed, 1E-8);
	}

	
	@Test
	public void test02SimpleValue() {
		double observed = PoissonUtility.getProbabilityMass(22, 12);
		Assert.assertEquals("Testing expected probability mass", 0.0030177762602241563, observed, 1E-8);
	}

	@Test
	public void test03CumulativeMass() {
		double observed = PoissonUtility.getCumulativeProbabilityMass(22, 12);
		Assert.assertEquals("Testing expected cumulative probabality mass", 0.9969526288558276, observed, 1E-8);
	}

	@Test
	public void test04Quantile() {
		int observed = PoissonUtility.getQuantile(0.9, 2);
		Assert.assertEquals("Testing a quantile", 4, observed);
	}

	@Test
	public void test05Quantile() {
		int observed = PoissonUtility.getQuantile(0.6, 20);
		Assert.assertEquals("Testing a quantile", 21, observed);
	}

	
	@Test
	public void test05StochasticMean() {
		double lambda = 2d;
		MonteCarloEstimate est = new MonteCarloEstimate();
		int nbRealizations = 50000;
		for (int i = 0; i < nbRealizations; i++) {
			int observed = PoissonUtility.getQuantile(StatisticalUtility.getRandom().nextDouble(), lambda);
			est.addRealization(new Matrix(1,1,observed,0));
		}
		double mean = est.getMean().getValueAt(0, 0);
		double variance = est.getVariance().getValueAt(0, 0);
		System.out.println("Expected mean = " + lambda + "; Actual mean = " + mean);
		Assert.assertEquals("Testing the mean", lambda, mean, 2.5E-2);
		double expectedVariance = lambda;
		System.out.println("Expected variance = " + expectedVariance + "; Actual variance = " + variance);
		Assert.assertEquals("Testing the variance", expectedVariance, variance, 1E-1);
	}

	@Test
	public void test06StochasticMean() {
		double lambda = 5d;
		MonteCarloEstimate est = new MonteCarloEstimate();
		int nbRealizations = 50000;
		for (int i = 0; i < nbRealizations; i++) {
			int observed = PoissonUtility.getQuantile(StatisticalUtility.getRandom().nextDouble(), lambda);
			est.addRealization(new Matrix(1,1,observed,0));
		}
		double mean = est.getMean().getValueAt(0, 0);
		double variance = est.getVariance().getValueAt(0, 0);
		System.out.println("Expected mean = " + lambda + "; Actual mean = " + mean);
		Assert.assertEquals("Testing the mean", lambda, mean, 4E-2);
		double expectedVariance = lambda;
		System.out.println("Expected variance = " + expectedVariance + "; Actual variance = " + variance);
		Assert.assertEquals("Testing the variance", expectedVariance, variance, 1E-1);
	}
}
