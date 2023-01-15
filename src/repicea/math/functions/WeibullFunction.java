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
import repicea.stats.distributions.utility.WeibullUtility;

public class WeibullFunction extends AbstractMathematicalFunction {

	private static final int K_INDEX = 0;
	private static final int LAMBDA_INDEX = 1;
	public static final int THETA_INDEX = 2;
	
	private final boolean isLocationParameterEnabled;
	
	/**
	 * Private constructor.
	 * @param k the shape parameter (must be greater than 0)
	 * @param lambda the scale parameter (must be greater than 0)
	 * @param isLocationParameterEnabled a boolean to inform on whether the function has two or three parameters.
	 */
	protected WeibullFunction(double k, double lambda, boolean isLocationParameterEnabled) {
		if (k < 0 || lambda < 0) {
			throw new InvalidParameterException("The k and lambda arguments must be strictly positive (ie. > 0)!");
		}
		setParameterValue(K_INDEX, k);
		setParameterValue(LAMBDA_INDEX, lambda);
		setVariableValue(0, 0d);
		setBounds(K_INDEX, new ParameterBound(MINIMUM_ACCEPTABLE_POSITIVE_VALUE, null));	// sigma2 must be strictly positive
		setBounds(LAMBDA_INDEX, new ParameterBound(MINIMUM_ACCEPTABLE_POSITIVE_VALUE, null));	// sigma2 must be strictly positive
		this.isLocationParameterEnabled = isLocationParameterEnabled;
	}

	/**
	 * Private constructor for two-parameter Weibull function.
	 * @param k the shape parameter (must be greater than 0)
	 * @param lambda the scale parameter (must be greater than 0)
	 */
	public WeibullFunction(double k, double lambda) {
		this(k, lambda, false);
	}
	
	/**
	 * Constructor for three-parameter Weibull function.
	 * @param k the shape parameter (must be greater than 0)
	 * @param lambda the scale parameter (must be greater than 0)
	 * @param theta the location parameter
	 */
	public WeibullFunction(double k, double lambda, double theta) {
		this(k, lambda, true);
		setParameterValue(THETA_INDEX, theta);
	}		

	@Override
	public void setParameterValue(int index, double value) {
		int nbParameters = isLocationParameterEnabled ? 3 : 2; 
		if (index > nbParameters - 1) {
			throw new InvalidParameterException("This Weibull function only has " + nbParameters + " parameters!");
		} else {
			super.setParameterValue(index, value);
		}
	}

	@Override
	public void setVariableValue(int index, double value) {
		if (index > 0) {
			throw new InvalidParameterException("The Weibull function only has one variable (namely the observation x)!");
		} else {
			double y = isLocationParameterEnabled ? value - getParameterValue(THETA_INDEX) : value;
			if (y < 0) {
				throw new InvalidParameterException("The Weibull function does not support null or negative values (when considering the location parameter)!");
			}
			super.setVariableValue(index, value);
		}
	}
	
	@Override
	public Double getValue() {
		double x = getVariableValue(0);
		double k = getParameterValue(K_INDEX);
		double lambda = getParameterValue(LAMBDA_INDEX);
		double y = isLocationParameterEnabled ? x - getParameterValue(THETA_INDEX) : x;
		return y < 0 ? 0 : WeibullUtility.getProbabilityDensity(y, k, lambda);
	}

	@Override
	public Matrix getGradient() {
		double x = getVariableValue(0);
		double k = getParameterValue(K_INDEX);
		double lambda = getParameterValue(LAMBDA_INDEX);
		double theta = isLocationParameterEnabled ?  getParameterValue(THETA_INDEX) : 0;
		double xMod = x - theta;
		double xMod_lambda = xMod/lambda;
		double xMod_lambda_powK = Math.pow(xMod_lambda, k);
		double df_dk = Math.exp(-xMod_lambda_powK) * xMod_lambda_powK * (1 - k * (xMod_lambda_powK - 1) * Math.log(xMod_lambda)) / xMod;
		double df_dLambda = k * k * Math.exp(-xMod_lambda_powK) * (xMod_lambda_powK - 1) * Math.pow(xMod_lambda, k-1) / (lambda * lambda);
		Matrix gradient;
		if (isLocationParameterEnabled) {
			gradient = new Matrix(3,1);
			gradient.setValueAt(K_INDEX, 0, df_dk);
			gradient.setValueAt(LAMBDA_INDEX, 0, df_dLambda);
			double df_dtheta = k * Math.exp(-xMod_lambda_powK) * (k * (xMod_lambda_powK - 1) + 1) * xMod_lambda_powK / (xMod * xMod);
			gradient.setValueAt(THETA_INDEX, 0, df_dtheta);
			return gradient;
			
		} else {
			gradient = new Matrix(2,1);
			gradient.setValueAt(K_INDEX, 0, df_dk);
			gradient.setValueAt(LAMBDA_INDEX, 0, df_dLambda);
			return gradient;
		}
		
	}

	@Override
	public SymmetricMatrix getHessian() {
		double x = getVariableValue(0);
		double k = getParameterValue(K_INDEX);
		double lambda = getParameterValue(LAMBDA_INDEX);
		double theta = isLocationParameterEnabled ?  getParameterValue(THETA_INDEX) : 0;
		double xMod = x - theta;
		double xMod_lambda = xMod/lambda;
		double xMod_lambda_powK = Math.pow(xMod_lambda, k);

		double d2_d2k = Math.exp(-xMod_lambda_powK) * xMod_lambda_powK * Math.log(xMod_lambda) * (-2 * xMod_lambda_powK + k*(-3 * xMod_lambda_powK + xMod_lambda_powK*xMod_lambda_powK + 1) * Math.log(xMod_lambda) + 2) / xMod;
		
		
		double d2_d2lambda = k * k * Math.exp(-xMod_lambda_powK) * (-xMod_lambda_powK + k*(-3*xMod_lambda_powK + xMod_lambda_powK*xMod_lambda_powK + 1) + 1) * Math.pow(xMod_lambda, k - 1) / (lambda * lambda * lambda);
		double d2_dkdLambda = - k * Math.exp(-xMod_lambda_powK) * Math.pow(xMod_lambda, k - 1) * (-2 * xMod_lambda_powK + k * (-3 * xMod_lambda_powK + xMod_lambda_powK * xMod_lambda_powK + 1) * Math.log(xMod_lambda)  + 2) / (lambda * lambda);
		
		SymmetricMatrix hessian;
		if (isLocationParameterEnabled) {
			
			double d2_d2theta = k * Math.exp(-xMod_lambda_powK) * (k *k *(-3 * xMod_lambda_powK + xMod_lambda_powK*xMod_lambda_powK + 1) + 3 * k * (xMod_lambda_powK - 1) + 2) *xMod_lambda_powK / (xMod * xMod * xMod);
			double d2_dkdTheta = Math.exp(-xMod_lambda_powK) * xMod_lambda_powK * (2 * k * (xMod_lambda_powK - 1) - k * (xMod_lambda_powK + k*(-3*xMod_lambda_powK + xMod_lambda_powK*xMod_lambda_powK + 1) - 1)*Math.log(xMod_lambda) + 1) / (xMod * xMod);
			double d2_dLambdadTheta = k*k* Math.exp(-xMod_lambda_powK) * (xMod_lambda_powK + k * (-3*xMod_lambda_powK + xMod_lambda_powK*xMod_lambda_powK + 1) -1) *Math.pow(xMod_lambda, k+1)  / (xMod * xMod * xMod);
			
			hessian = new SymmetricMatrix(3);
			hessian.setValueAt(K_INDEX, K_INDEX, d2_d2k);
			hessian.setValueAt(LAMBDA_INDEX, LAMBDA_INDEX, d2_d2lambda);
			hessian.setValueAt(K_INDEX, LAMBDA_INDEX, d2_dkdLambda);
//			hessian.setValueAt(LAMBDA_INDEX, K_INDEX, d2_dkdLambda);

			hessian.setValueAt(LAMBDA_INDEX, THETA_INDEX, d2_dLambdadTheta);
//			hessian.setValueAt(THETA_INDEX, LAMBDA_INDEX, d2_dLambdadTheta);

			hessian.setValueAt(K_INDEX, THETA_INDEX, d2_dkdTheta);
//			hessian.setValueAt(THETA_INDEX, K_INDEX, d2_dkdTheta);
			hessian.setValueAt(THETA_INDEX, THETA_INDEX, d2_d2theta);
			
			return hessian;
			
		} else {
			hessian = new SymmetricMatrix(2);
			hessian.setValueAt(K_INDEX, K_INDEX, d2_d2k);
			hessian.setValueAt(LAMBDA_INDEX, LAMBDA_INDEX, d2_d2lambda);
			hessian.setValueAt(K_INDEX, LAMBDA_INDEX, d2_dkdLambda);
//			hessian.setValueAt(LAMBDA_INDEX, K_INDEX, d2_dkdLambda);
			return hessian;
		}
	}

	
}

