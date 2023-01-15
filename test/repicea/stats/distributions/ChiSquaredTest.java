/*
 * This file is part of the repicea library.
 *
 * Copyright (C) 2009-2017 Mathieu Fortin for Rouge-Epicea
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
package repicea.stats.distributions;

import org.junit.Assert;
import org.junit.Test;

import repicea.math.Matrix;
import repicea.math.SymmetricMatrix;

public class ChiSquaredTest {

	@Test
	public void randomGenerationTestForUnivariateChiSquared() {
		int df = 100;
		double expectedMean = 5d;
		ChiSquaredDistribution dist = new ChiSquaredDistribution(df, expectedMean);
		EmpiricalDistribution receiver = new EmpiricalDistribution();
		for (int i = 0; i < 1000000; i++) {
			receiver.addRealization(dist.getRandomRealization());
		}
		
		Matrix mean = receiver.getMean();
		Assert.assertEquals("Testing the mean", dist.getMean().getValueAt(0, 0), mean.getValueAt(0, 0), 1E-2);
		Matrix variance = receiver.getVariance();
		Assert.assertEquals("Testing the variance", dist.getVariance().getValueAt(0, 0), variance.getValueAt(0, 0), 1E-2);
	}
	
	
	@Test
	public void randomGenerationTestForMultivariateChiSquared() {
		int df = 100;
		SymmetricMatrix expectedMean = new SymmetricMatrix(2);
		expectedMean.setValueAt(0, 0, 5d);
		expectedMean.setValueAt(1, 1, 5d);
		expectedMean.setValueAt(0, 1, 0d);
		
		ChiSquaredDistribution dist = new ChiSquaredDistribution(df, expectedMean);
		EmpiricalDistribution receiver = new EmpiricalDistribution();
		EmpiricalDistribution receiver00 = new EmpiricalDistribution();
		EmpiricalDistribution receiver01 = new EmpiricalDistribution();
		EmpiricalDistribution receiver10 = new EmpiricalDistribution();
		EmpiricalDistribution receiver11 = new EmpiricalDistribution();
		for (int i = 0; i < 1000000; i++) {
			SymmetricMatrix realization = dist.getRandomRealization();
			receiver.addRealization(realization);
			receiver00.addRealization(realization.getSubMatrix(0, 0, 0, 0));
			receiver01.addRealization(realization.getSubMatrix(0, 0, 1, 1));
			receiver10.addRealization(realization.getSubMatrix(1, 1, 0, 0));
			receiver11.addRealization(realization.getSubMatrix(1, 1, 1, 1));
		}
		
		Matrix obsMean = receiver.getMean();
		Matrix meanDiff = expectedMean.subtract(obsMean);
		Assert.assertTrue("Testing the mean", !meanDiff.getAbsoluteValue().anyElementLargerThan(5E-3));	// test if no absolute difference is larger than 5E-3

		Matrix obsVariance = new Matrix(2,2);
		obsVariance.setValueAt(0, 0, receiver00.getVariance().getValueAt(0, 0));
		obsVariance.setValueAt(0, 1, receiver01.getVariance().getValueAt(0, 0));
		obsVariance.setValueAt(1, 0, receiver10.getVariance().getValueAt(0, 0));
		obsVariance.setValueAt(1, 1, receiver11.getVariance().getValueAt(0, 0));
		Matrix expectedVariance = dist.getVariance();
		Matrix varianceDiff = expectedVariance.subtract(obsVariance);
		
		Assert.assertTrue("Testing the variance", !varianceDiff.getAbsoluteValue().anyElementLargerThan(5E-3)); 	// test if no absolute difference is larger than 5E-3
	}
	
}
