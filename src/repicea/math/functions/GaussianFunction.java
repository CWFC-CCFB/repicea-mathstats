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

import java.security.InvalidParameterException;

import repicea.math.AbstractMathematicalFunction;
import repicea.math.Matrix;
import repicea.math.ParameterBound;
import repicea.math.SymmetricMatrix;
import repicea.stats.distributions.utility.GaussianUtility;

/**
 * The GaussianFunction class implements the pdf of the normal distribution. <br>
 * <br>
 * The function has two parameters, namely mu (the mean) and sigma2 (the variance).
 * 
 * @author Mathieu Fortin - July 2022
 */
@SuppressWarnings("serial")
public class GaussianFunction extends AbstractMathematicalFunction {

	private static final int MU_INDEX = 0;
	private static final int SIGMA2_INDEX = 1;
	
	/**
	 * Constructor 1.
	 * @param mu the mean of the function
	 * @param sigma2 the variance of the function
	 */
	public GaussianFunction(double mu, double sigma2) {
		if (sigma2 <= 0) {
			throw new InvalidParameterException("The sigma2 argument must be strictly positive (ie. > 0)!");
		}
		setParameterValue(MU_INDEX, mu);
		setParameterValue(SIGMA2_INDEX, sigma2);
		setVariableValue(0, 0d);
		setBounds(SIGMA2_INDEX, new ParameterBound(MINIMUM_ACCEPTABLE_POSITIVE_VALUE, null));	// sigma2 must be strictly positive
	}

	/**
	 * A default contructor with mu and sigma2 set to 0 and 1 respectively.
	 */
	public GaussianFunction() {
		this(0,1);
	}

	@Override
	public void setParameterValue(int index, double value) {
		if (index > 1) {
			throw new InvalidParameterException("The Gaussian function only has two parameters!");
		} else {
			super.setParameterValue(index, value);
		}
	}

	@Override
	public void setVariableValue(int index, double value) {
		if (index > 0) {
			throw new InvalidParameterException("The Gaussian function only has one variable (namely the observation x)!");
		} else {
			super.setVariableValue(index, value);
		}
	}
	
	@Override
	public Double getValue() {
		double x = getVariableValue(0);
		double mu = getParameterValue(MU_INDEX);
		double sigma2 = getParameterValue(SIGMA2_INDEX);
		return GaussianUtility.getProbabilityDensity(x, mu, sigma2);
	}

	@Override
	public Matrix getGradient() {
		double x = getVariableValue(0);
		double mu = getParameterValue(MU_INDEX);
		double sigma2 = getParameterValue(SIGMA2_INDEX);
		
		double f = GaussianUtility.getProbabilityDensity(x, mu, sigma2);
		
		Matrix gradient = new Matrix(2,1);
		double df_dMu =  f * (x - mu) / sigma2;
		gradient.setValueAt(MU_INDEX, 0, df_dMu);
		double df_dSigma2 = f * ((x-mu) * (x-mu)/(2 * sigma2 * sigma2) - 1d / (2 * sigma2));
		gradient.setValueAt(SIGMA2_INDEX, 0, df_dSigma2);

		return gradient;
	}

	@Override
	public SymmetricMatrix getHessian() {
		double x = getVariableValue(0);
		double mu = getParameterValue(MU_INDEX);
		double sigma2 = getParameterValue(SIGMA2_INDEX);

		double f = GaussianUtility.getProbabilityDensity(x, mu, sigma2);

		Matrix gradient = getGradient();
		
		SymmetricMatrix hessian = new SymmetricMatrix(2);
		double d2f_d2Mu = gradient.getValueAt(MU_INDEX, 0) * (x - mu) / sigma2 - f / sigma2;
		hessian.setValueAt(MU_INDEX, MU_INDEX, d2f_d2Mu);
		double d2f_dMu_dSigma2 = gradient.getValueAt(SIGMA2_INDEX, 0) * (x - mu) / sigma2 - f * (x - mu) / (sigma2 * sigma2);
		hessian.setValueAt(MU_INDEX, SIGMA2_INDEX, d2f_dMu_dSigma2);
//		hessian.setValueAt(SIGMA2_INDEX, MU_INDEX, d2f_dMu_dSigma2);
		double d2f_d2Sigma2 = gradient.getValueAt(SIGMA2_INDEX, 0) * ((x-mu) * (x-mu)/(2 * sigma2 * sigma2) - 1d / (2* sigma2)) +
				f * (-(x-mu) * (x-mu)/(sigma2 * sigma2 * sigma2) + 1d / (2 * sigma2 * sigma2));
		hessian.setValueAt(SIGMA2_INDEX, SIGMA2_INDEX, d2f_d2Sigma2);
		return hessian;
	}

	
}
