/*
 * This file is part of the repicea-mathstats library.
 *
 * Copyright (C) 2023 His Majesty the King in Right of Canada 
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
package repicea.stats.distributions;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

import org.junit.Assert;
import org.junit.Test;

import repicea.math.ComplexMatrix;
import repicea.math.ComplexNumber;
import repicea.math.HermitianMatrix;
import repicea.math.Matrix;
import repicea.stats.StatisticalUtility;

public class ComplexEmpiricalDistributionTest {

	@Test
	public void simpleMeanAndVarianceTest() {
		int nbReals = 100000;
		ComplexEmpiricalDistribution ced = new ComplexEmpiricalDistribution();
		List<ComplexNumber> l = new ArrayList<ComplexNumber>();
		for (int i = 0; i < nbReals; i++) {
			l.clear();
			l.add(new ComplexNumber(2 + StatisticalUtility.getRandom().nextGaussian(), StatisticalUtility.getRandom().nextGaussian()));
			ced.addRealization(new ComplexMatrix(l));
		}
		
		ComplexMatrix mean = ced.getMean();
		double actualMeanRealPart = mean.getValueAt(0, 0).realPart;
		System.out.println("Expected mean (real part) = 2; Actual mean = " + actualMeanRealPart);
		Assert.assertEquals("Testing mean (real part)", 2, actualMeanRealPart, 1E-2);
		double actualMeanImagPart = mean.getValueAt(0, 0).imaginaryPart;
		System.out.println("Expected mean (imaginary part) = 0; Actual mean = " + actualMeanImagPart);
		Assert.assertEquals("Testing mean (imaginary part)", 0, actualMeanImagPart, 1E-2);
		HermitianMatrix variance = ced.getVariance();
		System.out.println("Expected variance (total) = 2; Actual variance = " + variance.getValueAt(0, 0).realPart);
		Assert.assertEquals("Testing variance", 2d, variance.getValueAt(0, 0).realPart, 1.5E-2);
		System.out.println("Expected variance (imaginary part) = 0; Actual variance = " + variance.getValueAt(0, 0).imaginaryPart);
		Assert.assertEquals("Testing variance", 0d, variance.getValueAt(0, 0).imaginaryPart, 1E-15);
	}
	
	
	
	@Test
	public void multivariateMeanAndVarianceTest() {
		int nbReals = 100000;
		ComplexEmpiricalDistribution ced = new ComplexEmpiricalDistribution();
		List<ComplexNumber> l = new ArrayList<ComplexNumber>();
		for (int i = 0; i < nbReals; i++) {
			l.clear();
			l.add(new ComplexNumber(2 + StatisticalUtility.getRandom().nextGaussian(), StatisticalUtility.getRandom().nextGaussian()));
			l.add(new ComplexNumber(4 + StatisticalUtility.getRandom().nextGaussian() * 2, StatisticalUtility.getRandom().nextGaussian()*.5));
			ced.addRealization(new ComplexMatrix(l));
		}
		List<Double> expectedMeanReal = Arrays.asList(new Double[] {2d,4d});
		List<Double> expectedMeanImag = Arrays.asList(new Double[] {0d,0d});
		ComplexMatrix mean = ced.getMean();
		for (int i = 0; i < mean.m_iRows; i++) {
			double actualMeanRealPart = mean.getValueAt(i, 0).realPart;
			System.out.println("Expected mean (real part) = " + expectedMeanReal.get(i) + "; Actual mean = " + actualMeanRealPart);
			Assert.assertEquals("Testing mean (real part)", expectedMeanReal.get(i), actualMeanRealPart, 5E-2);
			double actualMeanImagPart = mean.getValueAt(i, 0).imaginaryPart;
			System.out.println("Expected mean (imaginary part) = " + expectedMeanImag.get(i) + "; Actual mean = " + actualMeanImagPart);
			Assert.assertEquals("Testing mean (imaginary part)", expectedMeanImag.get(i), actualMeanImagPart, 1E-2);
		}
		
		HermitianMatrix variance = ced.getVariance();
		System.out.println("Expected variance (total) = 2; Actual variance = " + variance.getValueAt(0, 0).realPart);
		Assert.assertEquals("Testing variance", 2d, variance.getValueAt(0, 0).realPart, 2E-2);
		System.out.println("Expected variance (imaginary part) = 0; Actual variance = " + variance.getValueAt(0, 0).imaginaryPart);
		Assert.assertEquals("Testing variance", 0d, variance.getValueAt(0, 0).imaginaryPart, 1E-15);
		
		System.out.println("Expected variance (total) = 4.25; Actual variance = " + variance.getValueAt(1, 1).realPart);
		Assert.assertEquals("Testing variance", 4.25, variance.getValueAt(1, 1).realPart, 0.1);
		System.out.println("Expected variance (imaginary part) = 0; Actual variance = " + variance.getValueAt(1, 1).imaginaryPart);
		Assert.assertEquals("Testing variance", 0d, variance.getValueAt(1, 1).imaginaryPart, 1E-15);

		System.out.println("Expected variance (total) = 0; Actual variance = " + variance.getValueAt(0, 1).realPart);
		Assert.assertEquals("Testing variance", 0, variance.getValueAt(0, 1).realPart, 2E-2);
		System.out.println("Expected variance (imaginary part) = 0; Actual variance = " + variance.getValueAt(0, 1).imaginaryPart);
		Assert.assertEquals("Testing variance", 0d, Math.abs(variance.getValueAt(0, 1).imaginaryPart), 0.1);
	}

	
	@Test
	public void varianceDecompositionTest() {
		int nbReals = 100000;
		ComplexEmpiricalDistribution ced = new ComplexEmpiricalDistribution();
		List<ComplexNumber> l = new ArrayList<ComplexNumber>();
		for (int i = 0; i < nbReals; i++) {
			l.clear();
			l.add(new ComplexNumber(2 + StatisticalUtility.getRandom().nextGaussian(), StatisticalUtility.getRandom().nextGaussian()));
			l.add(new ComplexNumber(4 + StatisticalUtility.getRandom().nextGaussian() * 2, StatisticalUtility.getRandom().nextGaussian()*.5));
			ced.addRealization(new ComplexMatrix(l));
		}
		HermitianMatrix variance = ced.getVariance();
		Matrix realPart = ced.getVarianceRealPart();
		Matrix imagPart = ced.getVarianceImaginaryPart();
		for (int i = 0; i < variance.m_iRows; i++) {
			double totalVariance = variance.getValueAt(i, i).realPart;
			double rPart = realPart.getValueAt(i, 0);
			double iPart = imagPart.getValueAt(i, 0);
			Assert.assertEquals("Testing variance decomposition", totalVariance, rPart + iPart, 1E-8);
		}
		
	}

}
