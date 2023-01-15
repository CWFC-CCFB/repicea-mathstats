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
package repicea.math.functions;

import org.junit.Assert;
import org.junit.Test;

import repicea.math.Matrix;

/**
 * The GaussianFunctionTest class tests the density, the gradient and the hessian of
 * the GaussianFunction class.
 * @author Mathieu Fortin - 2022
 */
public class GaussianFunctionTest {

	@Test
	public void valueTest() {
		double mu = 2;
		double sigma2 = 3;
		double x = .5;
		GaussianFunction gf = new GaussianFunction(mu, sigma2);
		gf.setVariableValue(0, x);
		double observedValue = gf.getValue();
		double expectedValue = Math.exp(-(x - mu) * (x - mu) / (2 * sigma2))/ Math.sqrt(2 * Math.PI * sigma2);
		Assert.assertEquals("Comparing densities", expectedValue, observedValue, 1E-8);
	}
	
	
	@Test
	public void gradientTest() {
		double mu = 2;
		double sigma2 = 3;
		double x = .5;
		GaussianFunction gf = new GaussianFunction(mu, sigma2);
		gf.setVariableValue(0, x);
		Matrix gradient = gf.getGradient();
		double df_dmu = (x - mu) * Math.exp(-(x - mu) * (x - mu) / (2 * sigma2))/ (Math.sqrt(2 * Math.PI) * Math.pow(sigma2, 1.5));
		double df_dsigma2 = ((x - mu) * (x - mu) - sigma2) * Math.exp(-(x - mu) * (x - mu) / (2 * sigma2)) / (2 * Math.sqrt(2 * Math.PI) * Math.pow(sigma2, 2.5));
		
		Assert.assertEquals("Comparing mu derivatives", df_dmu, gradient.getValueAt(0, 0), 1E-8);
		Assert.assertEquals("Comparing sigma2 derivatives", df_dsigma2, gradient.getValueAt(1, 0), 1E-8);
	}

	
	@Test
	public void hessianTest() {
		double mu = 2;
		double sigma2 = 3;
		double x = .5;
		GaussianFunction gf = new GaussianFunction(mu, sigma2);
		gf.setVariableValue(0, x);
		Matrix hessian = gf.getHessian();
		double d2f_d2mu = (x - mu) * (x - mu) * Math.exp(-(x - mu) * (x - mu) / (2 * sigma2))/ (Math.sqrt(2 * Math.PI) * Math.pow(sigma2, 2.5)) -
				Math.exp(-(x - mu) * (x - mu) / (2 * sigma2))/ (Math.sqrt(2 * Math.PI) * Math.pow(sigma2, 1.5));
		
		double df_dmu_dsigma2 = (x - mu) * Math.exp(-(x - mu) * (x - mu) / (2 * sigma2)) * ((x - mu) * (x-mu) - 3 *sigma2)/ (2 * Math.sqrt(2 * Math.PI) * Math.pow(sigma2, 3.5));

		double df_d2sigma2 = (3 * sigma2 * sigma2 - 6 * sigma2 * (x - mu) * (x -mu) + (x - mu) * (x -mu) * (x - mu) * (x -mu)) *
				Math.exp(-(x - mu) * (x - mu) / (2 * sigma2)) / (4 * Math.sqrt(2 * Math.PI) * Math.pow(sigma2, 4.5));

		Assert.assertEquals("Comparing mu second derivatives", d2f_d2mu, hessian.getValueAt(0, 0), 1E-8);
		Assert.assertEquals("Comparing mu x sigma2 second derivatives", df_dmu_dsigma2, hessian.getValueAt(1, 0), 1E-8);
		Assert.assertEquals("Comparing mu x sigma2 second derivatives", df_dmu_dsigma2, hessian.getValueAt(0, 1), 1E-8);
		Assert.assertEquals("Comparing sigma2 second derivatives", df_d2sigma2, hessian.getValueAt(1, 1), 1E-8);
	}

}
