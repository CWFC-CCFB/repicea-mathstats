/*
 * This file is part of the repicea library.
 *
 * Copyright (C) 2009-2015 Mathieu Fortin for Rouge-Epicea
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

/**
 * The AdaptativeGaussHermiteQuadrature class implements the adaptative algorithm, i.e.
 * the maximum of the function is first found before integrating.
 * @author Mathieu Fortin - December 2015, July 2022
 *
 */
@SuppressWarnings("serial")
abstract class AdaptativeGaussHermiteQuadrature { //extends AbstractGaussHermiteQuadrature implements UnidimensionalIntegralApproximation<AdaptativeGaussHermiteQuadratureCompatibleFunction>,
													//											MultidimensionalIntegralApproximation<AdaptativeGaussHermiteQuadratureCompatibleFunction>{
	

//	/**
//	 * An interface for a better management of function rescaling.
//	 * @author Mathieu Fortin - July 2022
//	 * @see <a href=https://en.wikipedia.org/wiki/Gauss%E2%80%93Hermite_quadrature /a> Gaussian-Hermite quadrature 
//	 */
//	public static interface AdaptativeGaussHermiteQuadratureCompatibleFunction extends MathematicalFunction, GaussHermiteQuadratureCompatibleFunction<Double> {
//		
//		/**
//		 * Return the variance-covariance matrix.
//		 * @return a positive definite square matrix
//		 */
//		public Matrix getVariance();
//	}
//	
//	private static class InternalExponentialFunctionWrapper extends ExponentialFunctionWrapper implements AdaptativeGaussHermiteQuadratureCompatibleFunction {
//
//		final Matrix variance;
//		
//		private InternalExponentialFunctionWrapper(AdaptativeGaussHermiteQuadratureCompatibleFunction originalFunction, Matrix variance) {
//			super(originalFunction);
//			this.variance = variance;
//		}
//
//		@Override
//		public AdaptativeGaussHermiteQuadratureCompatibleFunction getOriginalFunction() {
//			return (AdaptativeGaussHermiteQuadratureCompatibleFunction) super.getOriginalFunction();
//		}
//		
//		@Override
//		public double convertFromGaussToOriginal(double x, double mu, int covarianceIndexI, int covarianceIndexJ) {
//			return getOriginalFunction().convertFromGaussToOriginal(x, mu, covarianceIndexI, covarianceIndexJ);
//		}
//
//		@Override
//		public double getIntegralAdjustment(int dimensions) {
//			// TODO Auto-generated method stub
//			return 0;
//		}
//
//		@Override
//		public Matrix getVariance() {
//			return variance;
//		}
//	}
//
//	/**
//	 * Constructor.
//	 * @param numberOfPoints the number of points the integral is based on.
//	 */
//	public AdaptativeGaussHermiteQuadrature(NumberOfPoints numberOfPoints) {
//		super(numberOfPoints);
//	}
//	
//	/**
//	 * Constructor for the LaplaceApproximation class.
//	 */
//	AdaptativeGaussHermiteQuadrature() {
//		super();
//	}
//	
//	
//	@Override
//	public double getMultiDimensionalIntegralApproximation(AdaptativeGaussHermiteQuadratureCompatibleFunction functionToEvaluate,
//											List<Integer> indices, 
//											boolean isParameter) {
//		if (!isParameter) {
//			throw new UnsupportedOperationException("The AdaptativeGaussianHermiteQuadrature class has not been implemented for integral over the variables yet!");
//		}
//		if (indices == null || indices.isEmpty()) {
//			throw new InvalidParameterException("The indices argument must be a non empty list of integers!");
//		} else {
//			int maxIndex = functionToEvaluate.getNumberOfParameters();
//			for (Integer index : indices) {
//				if (index < 0 || index >= maxIndex) {
//					throw new InvalidParameterException("One index is either negative or it exceeds the number of parameters in the function!");
//				}
//			}
//			InternalLogWrapperFunction functionToBeOptimized = new InternalLogWrapperFunction(functionToEvaluate, indices, functionToEvaluate.getVariance());
//			
//			NewtonRaphsonOptimizer nro = new NewtonRaphsonOptimizer();
//			try {
//				nro.optimize(functionToBeOptimized, indices);
//			} catch (OptimizationException e) {
//				e.printStackTrace();
//			}
//
//			Matrix newHessian = nro.getHessianAtMaximum();
//			Matrix varCov = newHessian.getInverseMatrix().scalarMultiply(-1d);
////			Matrix newLowerCholeskyTriangle = varCov.getLowerCholTriangle();
//
//			ExponentialFunctionWrapper efw = new InternalExponentialFunctionWrapper(functionToBeOptimized, varCov);
//			double approximation = getMultiDimensionIntegral(efw, indices, isParameter, newLowerCholeskyTriangle);
//			int dimensions = newHessian.m_iRows;
//			approximation *= Math.pow(2d, (2 * dimensions - 1) * .5) * Math.sqrt(varCov.getDeterminant());
//			return approximation;
//		}
//	}
//
//	@Override
//	protected double getOneDimensionIntegral(GaussHermiteQuadratureCompatibleFunction<Double> functionToEvaluate,
//			List<Integer> indices, 
//			boolean isParameter,
//			int startingIndex) {
//		int index = indices.get(startingIndex);
//		double originalValue = functionToEvaluate.getParameterValue(index);
//		double sum = 0;
//		double value;
//		for (int i = 0; i < getXValues().size(); i++) {
//			double z = getXValues().get(i);
//			double valueOnOriginalScale = functionToEvaluate.convertFromGaussToOriginal(z, originalValue, index, index);
////			double tmp =  z * standardDeviation * Math.sqrt(2d);
////			functionToEvaluate.setParameterValue(index, originalValue + tmp);
//			functionToEvaluate.setParameterValue(index, valueOnOriginalScale);
//			value = functionToEvaluate.getValue() * Math.exp(z*z) * getWeights().get(i);
//			sum += value;
//		}
//		functionToEvaluate.setParameterValue(index, originalValue);
//		return sum;
//	}
//	
//	
//	
//	@Override
//	public double getIntegralApproximation(AdaptativeGaussHermiteQuadratureCompatibleFunction functionToEvaluate, 
//			int index,
//			boolean isParameter) {
//		return functionToEvaluate.getIntegralAdjustment(1) * getOneDimensionIntegral(functionToEvaluate, 
//				Arrays.asList(new Integer[] {index}), 
//				isParameter, 
//				0);
//	}
//

	
}
