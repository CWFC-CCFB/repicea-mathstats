/*
 * This file is part of the repicea-mathstats library.
 *
 * Copyright (C) 2024 His Majesty the King in Right of Canada
 * Author: Mathieu Fortin, Canadian Wood Fibre Centre, Canadian Forest Service
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
package repicea.math;

import org.junit.Assert;
import org.junit.Test;

import repicea.stats.Distribution;
import repicea.stats.StatisticalUtility;

public class ComplexRandomNumberTest {

		@Test
		public void randomChiSquareIntegrationTest() {
			int upsilon = 25;
			ComplexNumber g = new ComplexNumber(0,-1);
			ComplexNumber mean = new ComplexNumber(0,0);
			int nbRealizations = 100000;
			for (int i = 0; i < nbRealizations; i++) {
				double Q = StatisticalUtility.getRandom().nextChiSquare(upsilon);
				mean = mean.add(g.multiply(Q/upsilon).exp());
			}
			mean = mean.multiply(1d/nbRealizations);
			ComplexNumber actual = g.multiply(-2d/upsilon).add(1d).log().multiply(-0.5 * upsilon).exp();
			System.out.println("Monte Carlo estimate = " + mean.toString());
			System.out.println("Analytical estimate = " + actual.toString());
			Assert.assertEquals("Testing real part", mean.realPart, actual.realPart, 0.004);
			Assert.assertEquals("Testing imaginary part", mean.imaginaryPart, actual.imaginaryPart, 0.002);
		}
	
		/*
		 * Test the limit of the complex Monte Carlo point estimator 
		 */
		@Test
		public void limitExpectationTest() {
			int upsilon = 25;
//			int nbRealizations = 100000000;
			int nbRealizations = 10000000;
			double sigma2Hat = 1;
			Matrix betaHat = new Matrix(2,1);
			betaHat.setValueAt(0, 0, 2);
			betaHat.setValueAt(1, 0, .25);
			
			SymmetricMatrix invXtX = new SymmetricMatrix(2);
			invXtX.setValueAt(0, 0, 0.285);
			invXtX.setValueAt(1, 1, 0.005);
			invXtX.setValueAt(1, 0, -0.035);
			
			Matrix invXtXChol = invXtX.getLowerCholTriangle();
			
			Matrix xMat = new Matrix(1,2);
			xMat.setValueAt(0, 0, 1d);
			xMat.setValueAt(0, 1, 3);
			double xBetaHat = xMat.multiply(betaHat).getValueAt(0, 0);

			ComplexNumber mean = new ComplexNumber(0,0);
			for (int i = 0; i < nbRealizations; i++) {
				double chiSquareDeviate = StatisticalUtility.getRandom().nextChiSquare(upsilon);
				ComplexNumber sigma2Hat_b = new ComplexNumber(sigma2Hat, sigma2Hat * (chiSquareDeviate/upsilon - 1));
				ComplexNumber sigmaHat_b = sigma2Hat_b.sqrt();
				Matrix betaDeviatesMat = invXtXChol.multiply(StatisticalUtility.drawRandomVector(betaHat.m_iRows, Distribution.Type.GAUSSIAN)); // here we get C * epsilon
				ComplexNumber xBetaDeviates = new ComplexNumber(0, xMat.multiply(betaDeviatesMat).getValueAt(0, 0)); // here we get 0 + x C epsilon i 
				xBetaDeviates = xBetaDeviates.multiply(sigmaHat_b); // here we get 0 + sqrt{sigma2Hat_b} x C epsilon i
				ComplexNumber meanTmp = xBetaDeviates.add(xBetaHat); // here we get x beta + sqrt{sigma2Hat_b} x C epsilon i
				ComplexNumber meanPlusCorrectionFactor = meanTmp.add(sigma2Hat_b.multiply(0.5)); // here we get x beta + + sqrt{sigma2Hat_b} x C epsilon i + 0.5 * sigma2_b  				
				mean = mean.add(meanPlusCorrectionFactor.exp());
			}
			mean = mean.multiply(1d/nbRealizations);

			double omegaValue = xMat.multiply(invXtX).multiply(xMat.transpose()).getValueAt(0, 0);
			
			ComplexNumber oneMinusI = new ComplexNumber(1,-1);
			ComplexNumber term2 = oneMinusI.multiply((1-omegaValue) * 0.5 * sigma2Hat);
			ComplexNumber term3 = new ComplexNumber(1,-(1-omegaValue) * sigma2Hat / upsilon);
			term3 = term3.log().multiply(-upsilon * .5);
			ComplexNumber limitMean = term2.add(term3).add(xBetaHat).exp();
			System.out.println("Mean of realizations = " + mean.toString());
			System.out.println("Limit = " + limitMean.toString());
			Assert.assertEquals("Testing real part", mean.realPart, limitMean.realPart, 0.01);
			Assert.assertEquals("Testing imaginary part", mean.imaginaryPart, limitMean.imaginaryPart, 0.01);
		}

}
