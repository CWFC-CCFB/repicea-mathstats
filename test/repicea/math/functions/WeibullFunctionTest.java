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
 * The WeibullFunctionTest class tests the density, the gradient and the hessian of
 * the WeibullFunction class.
 * @author Mathieu Fortin - 2022
 */
public class WeibullFunctionTest {

	@Test
	public void valueTest() {
		double k = 2;
		double lambda = 1.5;
		double theta = 3;
		double x = 5;
		WeibullFunction gf = new WeibullFunction(k, lambda, theta);
		gf.setVariableValue(0, x);
		double observedValue = gf.getValue();
		double expectedValue = k / lambda * Math.pow((x - theta)/lambda,k-1) * Math.exp(-Math.pow((x - theta)/lambda, k));
		Assert.assertEquals("Comparing densities", expectedValue, observedValue, 1E-8);
	}
	
	
	@Test
	public void gradientTest() {
		double k = 2;
		double lambda = 1.5;
		double theta = 3;
		double x = 5;
		WeibullFunction gf = new WeibullFunction(k, lambda, theta);
		gf.setVariableValue(0, x);
		Matrix gradient = gf.getGradient();
		double df_dk = Math.exp(-Math.pow((x-theta)/lambda, k)) * Math.pow((x-theta)/lambda, k) * (k * (Math.pow((x-theta)/lambda, k) - 1) * Math.log((x-theta)/lambda) - 1) / (theta - x);
		double df_dlambda = k * k * Math.exp(-Math.pow((x-theta)/lambda, k)) * (Math.pow((x-theta)/lambda, k)   - 1)  * Math.pow((x-theta)/lambda, k) / (lambda * (x - theta));
		double df_dt = k * Math.exp(-Math.pow((x-theta)/lambda, k)) * (k * (Math.pow((x-theta)/lambda, k)   - 1) + 1) * Math.pow((x-theta)/lambda, k) / ((x - theta) * (x - theta));
				
		Assert.assertEquals("Comparing k derivatives", df_dk, gradient.getValueAt(0, 0), 1E-8);
		Assert.assertEquals("Comparing lambda derivatives", df_dlambda, gradient.getValueAt(1, 0), 1E-8);
		Assert.assertEquals("Comparing theta derivatives", df_dt, gradient.getValueAt(2, 0), 1E-8);
	}

	
	@Test
	public void hessianTest() {
		double k = 2;
		double lambda = 3;
		double x = 2.5;
		double theta = 2;
		WeibullFunction gf = new WeibullFunction(k, lambda, theta);
		gf.setVariableValue(0, x);
		Matrix hessian = gf.getHessian();
		double xMod = x - theta;
		double xMod_lambda = xMod / lambda;
		double xMod_lambda_powK = Math.pow(xMod_lambda, k);
		
		double d2f_d2k = Math.exp(-xMod_lambda_powK) * xMod_lambda_powK * Math.log(xMod_lambda) * (-2 * xMod_lambda_powK + k *(-3 * xMod_lambda_powK + xMod_lambda_powK*xMod_lambda_powK + 1) * Math.log(xMod_lambda) + 2) / xMod;
		double d2f_dk_dlambda = - k * Math.exp(-xMod_lambda_powK) * Math.pow(xMod_lambda, k - 1) * (-2 * xMod_lambda_powK + k*(-3*xMod_lambda_powK + xMod_lambda_powK*xMod_lambda_powK + 1)* Math.log(xMod_lambda) + 2) / (lambda * lambda);
		double d2f_d2lambda = k * k * Math.exp(-xMod_lambda_powK) * Math.pow(xMod_lambda, k - 1) * (- xMod_lambda_powK + k*(-3*xMod_lambda_powK + xMod_lambda_powK*xMod_lambda_powK + 1) + 1) / (lambda * lambda* lambda);

		double d2f_dk_dTheta = Math.exp(-xMod_lambda_powK) * xMod_lambda_powK * (2*k*(xMod_lambda_powK - 1) - k*(xMod_lambda_powK + k*(-3*xMod_lambda_powK + xMod_lambda_powK*xMod_lambda_powK +1)-1)*Math.log(xMod_lambda) + 1)/(xMod * xMod);
		double d2f_dLambda_dTheta = k * k * Math.exp(-xMod_lambda_powK) * Math.pow(xMod_lambda, k + 1) * (xMod_lambda_powK + k*(-3*xMod_lambda_powK + xMod_lambda_powK*xMod_lambda_powK + 1) - 1) / (xMod * xMod * xMod);

		double d2f_d2theta = k * Math.exp(-xMod_lambda_powK) * xMod_lambda_powK * (k*k*(-3*xMod_lambda_powK + xMod_lambda_powK*xMod_lambda_powK+1 ) + 3*k*(xMod_lambda_powK -1)+2) / (xMod * xMod * xMod);
		
		Assert.assertEquals("Comparing k second derivatives", d2f_d2k, hessian.getValueAt(0, 0), 1E-8);
		Assert.assertEquals("Comparing k x lambda second derivatives", d2f_dk_dlambda, hessian.getValueAt(1, 0), 1E-8);
		Assert.assertEquals("Comparing k x lambda second derivatives", d2f_dk_dlambda, hessian.getValueAt(0, 1), 1E-8);
		Assert.assertEquals("Comparing lambda second derivatives", d2f_d2lambda, hessian.getValueAt(1, 1), 1E-8);
		
		Assert.assertEquals("Comparing k x theta derivatives", d2f_dk_dTheta, hessian.getValueAt(0, 2), 1E-8);
		Assert.assertEquals("Comparing k x theta derivatives", d2f_dk_dTheta, hessian.getValueAt(2, 0), 1E-8);
		
		Assert.assertEquals("Comparing lambda  x theta derivatives", d2f_dLambda_dTheta, hessian.getValueAt(1, 2), 1E-8);
		Assert.assertEquals("Comparing lambda  x theta derivatives", d2f_dLambda_dTheta, hessian.getValueAt(2, 1), 1E-8);

		Assert.assertEquals("Comparing theta second derivatives", d2f_d2theta, hessian.getValueAt(2, 2), 1E-8);
	}

}
