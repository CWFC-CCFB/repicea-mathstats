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


public class NegativeBinomialUtilityTest {

	@Test
	public void simpleValueTest() {
		double expected = NegativeBinomialUtility.getMassProbabilityOLD(2, 1, .5);
		double observed = NegativeBinomialUtility.getMassProbability(2, 1, .5);
		Assert.assertEquals("Testing the two methods", expected, observed, 1E-8);
	}

	@Test
	public void performanceTest() {
		for (int i = 0; i < 100000; i++) 
			NegativeBinomialUtility.getMassProbabilityOLD(2, 1, .5);

		long startTime = System.currentTimeMillis();
		for (int i = 0; i < 100000; i++) 
			NegativeBinomialUtility.getMassProbabilityOLD(2, 1, .5);
		long elapsedTimeReference = System.currentTimeMillis() - startTime;
		
		startTime = System.currentTimeMillis();
		for (int i = 0; i < 100000; i++) 
			NegativeBinomialUtility.getMassProbability(2, 1, .5);
		long elapsedTimeNewMethod = System.currentTimeMillis() - startTime;
		System.out.println("Performance: Reference time = " + elapsedTimeReference + " - New method = " + elapsedTimeNewMethod);
		Assert.assertTrue("Testing that the new method is faster", elapsedTimeNewMethod < elapsedTimeReference);
	}
	
	
	@Test
	public void simpleValueTest2() {
		double observed = NegativeBinomialUtility.getMassProbability(22, 1, 1);
		double observed2 = NegativeBinomialUtility.getMassProbabilityOLD(22, 1, 1);
		Assert.assertEquals("Testing the two methods", observed, observed2, 1E-8);
	}

	@Test
	public void zeroValueFasterImplementationTest() {
		double observed = NegativeBinomialUtility.getMassProbability(0, 1, 1);
		double observed2 = NegativeBinomialUtility.getMassProbabilityOLD(0, 1, 1);
		Assert.assertEquals("Testing the two methods", observed, observed2, 1E-8);
	}
	
	@Test
	public void quantileTest1() {
		int observed = NegativeBinomialUtility.getQuantile(0.5, 1, .8);
		Assert.assertEquals("Testing a quantile", 1, observed);
	}

	@Test
	public void meanTest1() {
		double mu = 1d;
		double theta = .8;
		MonteCarloEstimate est = new MonteCarloEstimate();
		int nbRealizations = 50000;
		for (int i = 0; i < nbRealizations; i++) {
			int observed = NegativeBinomialUtility.getQuantile(StatisticalUtility.getRandom().nextDouble(),
					mu, 
					theta);
			est.addRealization(new Matrix(1,1,observed,0));
		}
		double mean = est.getMean().getValueAt(0, 0);
		double variance = est.getVariance().getValueAt(0, 0);
		Assert.assertEquals("Testing the mean", mu, mean, 1E-2);
		Assert.assertEquals("Testing the variance", mu + theta*mu*mu, variance, 1E-1);
	}


	@Test
	public void meanTest2() {
		double mu = 1.5;
		double theta = .5;
		MonteCarloEstimate est = new MonteCarloEstimate();
		int nbRealizations = 50000;
		for (int i = 0; i < nbRealizations; i++) {
			int observed = NegativeBinomialUtility.getQuantile(StatisticalUtility.getRandom().nextDouble(),
					mu, 
					theta);
			est.addRealization(new Matrix(1,1,observed,0));
		}
		double mean = est.getMean().getValueAt(0, 0);
		double variance = est.getVariance().getValueAt(0, 0);
		Assert.assertEquals("Testing the mean", mu, mean, 2E-2);
		Assert.assertEquals("Testing the variance", mu + theta*mu*mu, variance, 1E-1);
	}

}
