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
package repicea.stats;

import org.junit.Assert;
import org.junit.Test;

import repicea.math.Matrix;
import repicea.stats.estimates.MonteCarloEstimate;

public class REpiceaRandomTest {

	@Test
	public void testGammaMean() {
		double dDispersionGamma = 0.61766452969979;
		double fGammaMean = 1.813311388045523;
//		double scale = dDispersionGamma;
//		double shape = fGammaMean / dDispersionGamma;
		double shape = dDispersionGamma;
		double scale = fGammaMean / dDispersionGamma;
		REpiceaRandom randomGenerator = new REpiceaRandom();
		
		double mean = 0;
		int maxIter = 500000;
		double meanFactor = 1d / maxIter;
		for (int i = 0; i < maxIter; i++) {
			double randomDeviate = 0;
			try {
				randomDeviate = randomGenerator.nextGamma(shape, scale);
				mean += randomDeviate * meanFactor;
			} catch (Exception e) {
			}

		}
		
		System.out.println ("Simulated mean = " + mean + "; Expected mean = " + fGammaMean);
		Assert.assertEquals("Testing mean for gamma random values", fGammaMean, mean, 1E-2);
	}

	
	@Test
	public void testGammaVariance() {
		double dDispersionGamma = 0.61766452969979;
		double fGammaMean = 1.813311388045523;
//		double scale = dDispersionGamma;
//		double shape = fGammaMean / dDispersionGamma;
		double shape = dDispersionGamma;
		double scale = fGammaMean / dDispersionGamma;
		REpiceaRandom randomGenerator = new REpiceaRandom();
//		RandomDataImpl randomGenerator = new RandomDataImpl();
		double variance = fGammaMean * fGammaMean / dDispersionGamma;
		
		int maxIter = 500000;
		MonteCarloEstimate estimate = new MonteCarloEstimate();
		Matrix realization;
		for (int i = 0; i < maxIter; i++) {
			realization = new Matrix(1,1);
			realization.setValueAt(0, 0, randomGenerator.nextGamma(shape, scale));
			estimate.addRealization(realization);
		}
		double actual = estimate.getVariance().getValueAt(0, 0);
		@SuppressWarnings("unused")
		double mean = estimate.getMean().getValueAt(0, 0);
		System.out.println ("Simulated variance = " + actual + "; Expected variance = " + variance);
		Assert.assertEquals("Testing mean for gamma random values", variance, actual, 0.07);
	}

	@Test
	public void testNegativeBinomialMean() {
		double dispersion = 1.41209128903651;
		double expectedMean = 0.5172784824575968;	
		REpiceaRandom randomGenerator = new REpiceaRandom();
		
		double mean = 0;
		int maxIter = 500000;
		double meanFactor = 1d / maxIter;
		for (int i = 0; i < maxIter; i++) {
			double randomDeviate = 0;
			try {
				randomDeviate = randomGenerator.nextNegativeBinomial(expectedMean, dispersion);
				mean += randomDeviate * meanFactor;
			} catch (Exception e) {
			}

		}
		
		System.out.println ("Simulated mean = " + mean + "; Expected mean = " + expectedMean);
		Assert.assertEquals("Testing mean for gamma random values", expectedMean, mean, 1E-2);
	}

	
	@Test
	public void testNegativeBinomialVariance() {
		double dispersion = 1.41209128903651;
		double expectedMean = 0.5172784824575968;	
		REpiceaRandom randomGenerator = new REpiceaRandom();
		double expectedVariance = expectedMean + dispersion * expectedMean * expectedMean;
		int maxIter = 500000;
		MonteCarloEstimate estimate = new MonteCarloEstimate();
		Matrix realization;
		for (int i = 0; i < maxIter; i++) {
			realization = new Matrix(1,1);
			realization.setValueAt(0, 0, randomGenerator.nextNegativeBinomial(expectedMean, dispersion));
			estimate.addRealization(realization);
		}
		double actual = estimate.getVariance().getValueAt(0, 0);
		@SuppressWarnings("unused")
		double mean = estimate.getMean().getValueAt(0, 0);
		System.out.println ("Simulated variance = " + actual + "; Expected variance = " + expectedVariance);
		Assert.assertEquals("Testing mean for gamma random values", expectedVariance, actual, 1E-1);
	}

	@Test
	public void testBetaMeanAndVariance() {
		double scale1 = 1d;
		double scale2 = 2d;	
		double expectedMean = scale1 / (scale1 + scale2);
		double expectedVariance = scale1 * scale2 / ((scale1 + scale2) * (scale1 + scale2) * (scale1 + scale2 + 1));
		REpiceaRandom randomGenerator = new REpiceaRandom();
		int maxIter = 500000;
		MonteCarloEstimate estimate = new MonteCarloEstimate();
		Matrix realization;
		for (int i = 0; i < maxIter; i++) {
			realization = new Matrix(1,1);
			realization.setValueAt(0, 0, randomGenerator.nextBeta(scale1, scale2));
			estimate.addRealization(realization);
		}
		double actualMean = estimate.getMean().getValueAt(0, 0);
		double actualVariance = estimate.getVariance().getValueAt(0, 0);
		System.out.println ("Simulated mean = " + actualMean + "; Expected variance = " + expectedMean);
		Assert.assertEquals("Testing mean for gamma random values", expectedMean, actualMean, 5E-3);
		System.out.println ("Simulated variance = " + actualVariance + "; Expected variance = " + expectedVariance);
		Assert.assertEquals("Testing mean for gamma random values", expectedVariance, actualVariance, 1E-3);
	}


}
