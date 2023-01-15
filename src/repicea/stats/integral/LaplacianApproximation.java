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

import java.security.InvalidParameterException;
import java.util.List;

import repicea.math.AbstractMathematicalFunction;
import repicea.math.Matrix;
import repicea.math.optimizer.AbstractOptimizer.OptimizationException;
import repicea.math.optimizer.NewtonRaphsonOptimizer;

/**
 * The LaplacianApproximation class implements the Laplace approximation for integrals.
 * 
 * @author Mathieu Fortin - December 2015
 *
 */
@SuppressWarnings("serial")
abstract class LaplacianApproximation extends AdaptativeGaussHermiteQuadrature {

//	/**
//	 * Constructor.
//	 */
//	public LaplacianApproximation() {
//		super();
//		weights.add(Math.sqrt(2d * Math.PI));
//		xValues.add(0d);
//	}
//
//	@Override
//	public List<Double> getWeights() {return weights;}
//
//	@Override
//	public List<Double> getXValues() {return xValues;}
//	
//	
//	@Override
//	public double getMultiDimensionalIntegralApproximation(AbstractMathematicalFunction functionToEvaluate,
//			List<Integer> parameterIndices, 
//			boolean isParameter,
//			Matrix lowerCholeskyTriangle) {
//		if (!isParameter) {
//			throw new UnsupportedOperationException("The LaplacianApproximation class has not been implemented for integral over the variables yet!");
//		}
//		if (!lowerCholeskyTriangle.isSquare() || parameterIndices.size() != lowerCholeskyTriangle.m_iRows) {
//			throw new InvalidParameterException("The indices are not compatible with the lower Cholesky triangle!");
//		} else {
//			for (Integer index : parameterIndices) {
//				if (index < 0 || index >= functionToEvaluate.getNumberOfParameters()) {
//					throw new InvalidParameterException("One index is either negative or it exceeds the number of parameters in the function!");
//				}
//			}
//			Matrix matrixG = lowerCholeskyTriangle.multiply(lowerCholeskyTriangle.transpose());
//			InternalLogWrapperFunction functionToBeOptimized = new InternalLogWrapperFunction(functionToEvaluate, parameterIndices, matrixG);
//			NewtonRaphsonOptimizer nro = new NewtonRaphsonOptimizer();
//			try {
//				nro.optimize(functionToBeOptimized, parameterIndices);
//			} catch (OptimizationException e) {
//				e.printStackTrace();
//			}
//			Matrix newHessian = nro.getHessianAtMaximum();
//			int dimensions = parameterIndices.size();
//			double fOptimal = functionToBeOptimized.getValue();
//			double approximation = Math.pow(getWeights().get(0), dimensions) * Math.pow(newHessian.scalarMultiply(-1d).getDeterminant(), -.5) * Math.exp(fOptimal); 
//			return approximation;
//		}
//	}
//
//
////	public static void main(String[] args) {
////		Random random = new Random();
////		LinkFunction logit = new LinkFunction(LinkFunction.Type.Logit);
////		double xBeta = -1.5;
////		logit.setParameterValue(0, xBeta);
////		logit.setVariableValue(0, 1d);
////		double mean = 0;
////		int nbIter = 1000000;
////		double factor = 1d / nbIter;
////		double stdDev = 1d;
////		for (int i = 0; i < nbIter; i++) {
////			logit.setParameterValue(0, xBeta + random.nextGaussian() * stdDev);
////			mean += logit.getValue() * factor;
////		}
////
////		Matrix lowerCholeskyTriangle = new Matrix(1,1);
////		lowerCholeskyTriangle.m_afData[0][0] = 1d;
////		
////		System.out.println("Simulated mean =  " + mean);
////
////		logit.setParameterValue(0, xBeta);
////		
////		
////		List<Integer> parameterIndices = new ArrayList<Integer>();
////		parameterIndices.add(0);
////
////		LaplaceApproximation la = new LaplaceApproximation();
////		double sum = la.getIntegralApproximation(logit, parameterIndices, lowerCholeskyTriangle);
////		int u = 0;
////	}
//

}
