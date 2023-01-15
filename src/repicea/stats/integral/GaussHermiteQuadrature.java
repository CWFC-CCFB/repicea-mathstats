/*
 * This file is part of the repicea-statistics library.
 *
 * Copyright (C) 2009-2012 Mathieu Fortin for Rouge-Epicea
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
package repicea.stats.integral;

import java.security.InvalidParameterException;
import java.util.Arrays;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import repicea.math.EvaluableFunction;
import repicea.math.Matrix;
import repicea.stats.integral.GaussHermiteQuadrature.GaussHermiteQuadratureCompatibleFunction;

/**
 * The GaussHermiteQuadrature class provides the x values and their weights for numerical integration.
 * The class implements the 5-point, 10-point and 15-point quadrature.
 * @author Mathieu Fortin - July 2012
 */
@SuppressWarnings("serial")
public class GaussHermiteQuadrature extends AbstractGaussHermiteQuadrature implements UnidimensionalIntegralApproximation<GaussHermiteQuadratureCompatibleFunction<Double>>,
																				UnidimensionalIntegralApproximationForMatrix<GaussHermiteQuadratureCompatibleFunction<Matrix>>,
																				MultidimensionalIntegralApproximation<GaussHermiteQuadratureCompatibleFunction<Double>> {
	
	
	/**
	 * An interface for a better management of function rescaling.
	 * @author Mathieu Fortin - July 2022
	 * @see <a href=https://en.wikipedia.org/wiki/Gauss%E2%80%93Hermite_quadrature /a> Gaussian-Hermite quadrature 
	 */
	public static interface GaussHermiteQuadratureCompatibleFunction<P> extends EvaluableFunction<P> {
		
		/**
		 * Convert the value to the original scale.<br>
		 * <br>
		 * In most cases, the function is not e^-(x^2) and consequently, the function
		 * has to be rescaled. For instance, with the pdf of the normal distribution, 
		 * x = (y - mu)/((2 * sigma2)^(1/2)). From x, we can calculate the y on the 
		 * original scale as y = (2 * sigma2)^(1/2) * x + mu.
		 * 
		 * @param x the value used in the Gauss-Hermite quadrature.
		 * @param mu the mu value on the original scale
		 * @param covarianceIndexI the row index where to find the std in the Cholesky matrix.
		 * @param covarianceIndexJ the column index where to find the std in the Cholesky matrix.
		 * @return the value on the original scale.
		 * @see <a href=https://en.wikipedia.org/wiki/Gauss%E2%80%93Hermite_quadrature /a> Gaussian-Hermite quadrature 
		 */
		public double convertFromGaussToOriginal(double x, double mu, int covarianceIndexI, int covarianceIndexJ);
		
		/**
		 * Calculate the adjustment to the integral. <br>
		 * <br>
		 * For instance in the context of a Gaussian distribution, the relationship between x and y 
		 * is y = (2 * sigma2)^(1/2) * x + mu. Consequently, dy = (2 * sigma2)^(1/2) * dx and the integral
		 * adjustment reductes to to 1/PI^(-N/2) where N is the number of dimensions. This is the
		 * default implementation.
		 * 
		 * @param dimensions the number of dimensions
		 * @return a double
		 * @see <a href=https://en.wikipedia.org/wiki/Gauss%E2%80%93Hermite_quadrature /a> Gaussian-Hermite quadrature 
		 */
		public default double getIntegralAdjustment(int dimensions) {
			return Math.pow(Math.PI, -dimensions/2d);
		}
		
	}

	
	/**
	 * Constructor.
	 * @param numberOfPoints a NumberOfPoints enum variable (either NumberOfPoints.N5, NumberOfPoints.N10, or NumberOfPoints.N15) 
	 */
	public GaussHermiteQuadrature(NumberOfPoints numberOfPoints) {
		super(numberOfPoints);
	}
	
	/**
	 * Constructor. Default Gauss-Hermite quadrature with 5 points.
	 */
	public GaussHermiteQuadrature() {
		this(NumberOfPoints.N5);
	}
	

	@Override
	protected double getOneDimensionIntegral(GaussHermiteQuadratureCompatibleFunction<Double> functionToEvaluate,
												List<Integer> indices, 
												boolean isParameter,
												int startingIndex) {
		Integer thisIndex = indices.get(startingIndex);
		double originalValue = isParameter ? functionToEvaluate.getParameterValue(thisIndex) : functionToEvaluate.getVariableValue(thisIndex);
		double sum = 0;
		double value;
		for (int i = 0; i < getXValues().size(); i++) {
			double thisValueOnTheOriginalScale = functionToEvaluate.convertFromGaussToOriginal(getXValues().get(i), originalValue, startingIndex, startingIndex); // the first variance in the vcov matrix
			if (isParameter) {
				functionToEvaluate.setParameterValue(thisIndex, thisValueOnTheOriginalScale);
			} else {
				functionToEvaluate.setVariableValue(thisIndex, thisValueOnTheOriginalScale);
			}
			value = functionToEvaluate.getValue() * getWeights().get(i);
			sum += value;
		}
		functionToEvaluate.setParameterValue(thisIndex, originalValue);
		return sum;
	}
	
	@Override
	public double getMultiDimensionalIntegralApproximation(GaussHermiteQuadratureCompatibleFunction<Double> functionToEvaluate,
											List<Integer> indices, 
											boolean isParameter) {
		if (indices == null || indices.isEmpty()) {
			throw new InvalidParameterException("The indices argument must be a non empty list of integers!");
		} else {
			int maxIndex = isParameter ? functionToEvaluate.getNumberOfParameters() : functionToEvaluate.getNumberOfVariables();
			for (Integer index : indices) {
				if (index < 0 || index >= maxIndex) {
					throw new InvalidParameterException("One index is either negative or it exceeds the number of " + (isParameter ? "parameters" : "variables") + " in the function!");
				}
			}
			return functionToEvaluate.getIntegralAdjustment(indices.size()) * getMultiDimensionIntegral(functionToEvaluate, indices, isParameter, 0);
		}
	}

	@Override
	public double getIntegralApproximation(GaussHermiteQuadratureCompatibleFunction<Double> functionToEvaluate, 
			int index,
			boolean isParameter) {
		return functionToEvaluate.getIntegralAdjustment(1) * getOneDimensionIntegral(functionToEvaluate, 
				Arrays.asList(new Integer[] {index}), 
				isParameter, 
				0);
	}

	@Override
	public Matrix getIntegralApproximationForMatrixFunction(GaussHermiteQuadratureCompatibleFunction<Matrix> functionToEvaluate, 
			int index,
			boolean isParameter) {
		return getOneDimensionIntegralForMatrix(functionToEvaluate, 
				Arrays.asList(new Integer[] {index}), 
				isParameter, 
				0).scalarMultiply(functionToEvaluate.getIntegralAdjustment(1));

	}

	private Matrix getOneDimensionIntegralForMatrix(GaussHermiteQuadratureCompatibleFunction<Matrix> functionToEvaluate,
			List<Integer> indices, 
			boolean isParameter, 
			int startingIndex) {
		Integer thisIndex = indices.get(startingIndex);
		double originalValue = isParameter ? functionToEvaluate.getParameterValue(thisIndex) : functionToEvaluate.getVariableValue(thisIndex);
		Matrix sum = null;
		Matrix value;
		for (int i = 0; i < getXValues().size(); i++) {
			double thisValueOnTheOriginalScale = functionToEvaluate.convertFromGaussToOriginal(getXValues().get(i), originalValue, startingIndex, startingIndex); // the first variance in the vcov matrix
			if (isParameter) {
				functionToEvaluate.setParameterValue(thisIndex, thisValueOnTheOriginalScale);
			} else {
				functionToEvaluate.setVariableValue(thisIndex, thisValueOnTheOriginalScale);
			}
			value = functionToEvaluate.getValue().scalarMultiply(getWeights().get(i));
			sum =  sum == null ? value : sum.add(value);
		}
		functionToEvaluate.setParameterValue(thisIndex, originalValue);
		return sum;
	}

	
	
	
	
}
