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
package repicea.math.optimizer;

import java.security.InvalidParameterException;
import java.util.List;
import java.util.logging.Level;

import repicea.math.MathematicalFunction;
import repicea.math.Matrix;
import repicea.math.SymmetricMatrix;
import repicea.math.utility.MatrixUtility;
import repicea.util.REpiceaLogManager;

/**
 * The NewtonRaphsonOptimizer class implements the Optimizer interface. It optimizes a log-likelihood function using the
 * Newton-Raphson algorithm. The vector of parameter changes is estimated as - inv(d2LLK/d2Beta) * dLLK/dBeta.
 * @author Mathieu Fortin - June 2011
 */
public class NewtonRaphsonOptimizer extends AbstractOptimizer {
	
	public final static String InnerIterationStarted = "InnerIterationStarted";

	public static String LOGGER_NAME = "NewtonRaphsonOptimizer";
	
	protected int maxNumberOfIterations = 50;
	protected double gradientCriterion = 1E-3;
	private int iterationID;
	private LineSearchMethod lineSearchMethod;
	
	public NewtonRaphsonOptimizer(LineSearchMethod lineSearchMethod) {
		this.convergenceCriterion = 1E-8; // default value
		setLineSearchMethod(lineSearchMethod);
	}
	
	public NewtonRaphsonOptimizer() {
		this(null);
	}
	
	/**
	 * Sets the line search method. <br>
	 * <br>
	 *  
	 * If the lsm parameter is null,
	 * the line search method is set to LineSearchMethod.TEN_EQUAL 
	 * by default.
	 * @param lsm a LineSearchMethod enum
	 */
	public void setLineSearchMethod(LineSearchMethod lsm) {
		lineSearchMethod = lsm == null ? LineSearchMethod.TEN_EQUAL : lsm;
	}

	/**
	 * @param model the model that provides the log-likelihood function
	 */
	
	
	/**
	 * Optimize the log-likelihood function using the Newton-Raphson algorithm.
	 * @param function an MathematicalFunction instance
	 * @param indicesOfParametersToOptimize a list of the indices of the parameters to be optimized
	 * @param originalBeta the vector that contains the parameters of the previous outer optimization
	 * @param optimisationStep the optimization step from the Newton-Raphson algorithm
	 * @param previousLogLikelihood the value of the log-likelihood function in the last outer optimization
	 * @return the value of the function
	 * @throws OptimizationException if the inner optimization is not successful
	 */
	protected double runInnerOptimisation(MathematicalFunction function, 
			List<Integer> indicesOfParametersToOptimize,
			Matrix originalBeta, 
			Matrix optimisationStep,
			double previousLogLikelihood) throws OptimizationException {
		
		int numberSubIter = 0;
		int maxNumberOfSubiterations = (Integer) lineSearchMethod.getParameter();
		
		double scalingFactor;
		double currentLlkValue;
		
		Matrix newBeta;
		boolean areParametersWithinBounds;
		do {
			fireOptimizerEvent(NewtonRaphsonOptimizer.InnerIterationStarted);
			scalingFactor = getScalingFactor(numberSubIter);
			newBeta = originalBeta.add(optimisationStep.scalarMultiply(scalingFactor));
			areParametersWithinBounds = areParametersWithinBounds(function, newBeta);
			if (areParametersWithinBounds) {
				setParameters(function, indicesOfParametersToOptimize, newBeta);
				currentLlkValue = function.getValue();
				REpiceaLogManager.logMessage(LOGGER_NAME, Level.FINEST, LOGGER_NAME, "Subiteration : " + numberSubIter + "; Parms = " + newBeta.toString() + "; Log-likelihood : " + currentLlkValue);
			} else {
				REpiceaLogManager.logMessage(LOGGER_NAME, Level.FINEST, LOGGER_NAME, "Subiteration : " + numberSubIter + "; Parms = " + newBeta.toString() + "; Some parameters exceed bounds!");
				currentLlkValue = Double.NaN;
			}
			numberSubIter++;
		} while ((Double.isNaN(currentLlkValue) || currentLlkValue < previousLogLikelihood) && numberSubIter < maxNumberOfSubiterations); // loop if the number of iterations is not over the maximum number and either the likelihood is still higher or non defined
		
		if (Double.isNaN(currentLlkValue) ||  currentLlkValue < previousLogLikelihood) {
			throw new OptimizationException("Failed to improve the log-likelihood function !");
		} else {
			return currentLlkValue;
		}
	}
	
	private boolean areParametersWithinBounds(MathematicalFunction function, Matrix beta) {
		for (int i = 0; i < beta.m_iRows; i++) {
			if (!function.isThisParameterValueWithinBounds(i, beta.getValueAt(i, 0))) {
				return false;
			}
		}
		return true;
	}
	
	private double getScalingFactor(int numberSubIter) {
		switch(lineSearchMethod) {
		case TEN_EQUAL:
			return 1d - numberSubIter * .1;
		case SINGLE_TRIAL:
			return 1d;
		case HALF_STEP:
			return 1d * Math.pow(0.5, numberSubIter);
		default:
			throw new InvalidParameterException("The line search method " + lineSearchMethod.name() + " has not been implemented yet!");
		}
	}

	@Override
	public boolean optimize(MathematicalFunction function, List<Integer> indicesOfParametersToOptimize) throws OptimizationException {
		convergenceAchieved = false;

		if (function instanceof OptimizerListener) {
			addOptimizerListener((OptimizerListener) function);
		}
		
		fireOptimizerEvent(OptimizerListener.optimizationStarted);
		
		double value0 = function.getValue();
		REpiceaLogManager.logMessage(LOGGER_NAME, Level.FINE, LOGGER_NAME, "Initial parameter estimates = " + function.getParameters().toString() + "; Initial LLK = " + value0);
		Matrix gradient = function.getGradient();
		REpiceaLogManager.logMessage(LOGGER_NAME, Level.FINER, LOGGER_NAME, "Gradient = " + gradient.toString());
		SymmetricMatrix hessian = function.getHessian();
				
		iterationID = 0;
		
		double gconv = calculateConvergence(gradient, hessian, value0);
		if (Math.abs(gconv) < convergenceCriterion) {
			convergenceAchieved = true;
		}
		Matrix currentBeta = function.getParameters();
		double previousLLK = 0;
		double fconv;
		try {
			while (!convergenceAchieved && iterationID <= maxNumberOfIterations) {
				iterationID++;
				Matrix inverseHessian = hessian.getInverseMatrix();
				Matrix optimisationStep = inverseHessian.multiply(gradient).scalarMultiply(-1d);
				REpiceaLogManager.logMessage(LOGGER_NAME, Level.FINEST, LOGGER_NAME, "Optimization step at iteration " + iterationID + " = " + optimisationStep.toString());

				Matrix originalBeta = extractParameters(function,indicesOfParametersToOptimize);
				previousLLK = value0;	
				value0 = runInnerOptimisation(function, 
						indicesOfParametersToOptimize, 
						originalBeta, 
						optimisationStep, 
						value0);		// if it does not throw an Exception, it means the inner optimisation was successful. 
				REpiceaLogManager.logMessage(LOGGER_NAME, Level.FINE, LOGGER_NAME, "LLK at iteration " + iterationID + " = " + value0);

				gradient = function.getGradient();
				hessian = function.getHessian();
				currentBeta = extractParameters(function, indicesOfParametersToOptimize);
				boolean originalHessian = true;
				gconv = calculateConvergence(gradient, hessian, value0);
				fconv = Math.abs(value0 - previousLLK) / Math.abs(previousLLK);
				int iterForHessianCorrection = 0;
				while (gconv < 0 && iterForHessianCorrection < 1000) {
					originalHessian = false;
					hessian = (SymmetricMatrix) hessian.add(Matrix.getIdentityMatrix(hessian.m_iRows));
					gconv = calculateConvergence(gradient, hessian, value0);
					iterForHessianCorrection++;
				}
				if (!originalHessian) {
					REpiceaLogManager.logMessage(LOGGER_NAME, Level.FINER, LOGGER_NAME, "Hessian was not positive-definite!");
				}
				REpiceaLogManager.logMessage(LOGGER_NAME, Level.FINER, LOGGER_NAME, "Iteration : " + iterationID + "; Log-likelihood : " + value0 + "; gconv : " + gconv + "; parms : " + currentBeta.toString() + "; gradient : " + gradient.toString());
				
				if (gconv < 0) {
					REpiceaLogManager.logMessage(LOGGER_NAME, Level.FINE, LOGGER_NAME, "GConv is negative!");
					convergenceAchieved = !gradient.getAbsoluteValue().anyElementLargerThan(1E-5) || fconv < 1E-8;
					break;
				} else if (gconv < convergenceCriterion) {
					convergenceAchieved = true;
				}
			}  

			if (iterationID > maxNumberOfIterations && !convergenceAchieved) {
				throw new OptimizationException("Convergence could not be achieved after " + ((Integer) iterationID).toString() + " iterations!");
			}

			if (gradient.getAbsoluteValue().anyElementLargerThan(gradientCriterion)) {
				REpiceaLogManager.logMessage(LOGGER_NAME, Level.FINE, LOGGER_NAME, "A least one element of the gradient vector is larger than 1E-3.");
			}

			optimalValue = value0;

			fireOptimizerEvent(OptimizerListener.optimizationEnded);

			return convergenceAchieved;
		} catch (OptimizationException e) {
			e.printStackTrace();
			throw e;
		} finally {
			betaVector = currentBeta;
			hessianMatrix = hessian;
			if (function instanceof OptimizerListener) {
				removeOptimizerListener((OptimizerListener) function);
			}
		}
	}
	
	private Matrix extractParameters(MathematicalFunction function, List<Integer> indicesOfParametersToOptimize) {
		Matrix beta = new Matrix(indicesOfParametersToOptimize.size(),1);
		for (int i = 0; i < indicesOfParametersToOptimize.size(); i++) {
			int parameterIndex = indicesOfParametersToOptimize.get(i);
			beta.setValueAt(i, 0, function.getParameterValue(parameterIndex));
		}
		return beta;
	}

	private void setParameters(MathematicalFunction function, List<Integer> indicesOfParametersToOptimize, Matrix newBeta) {
		for (int i = 0; i < indicesOfParametersToOptimize.size(); i++) {
			int index = indicesOfParametersToOptimize.get(i);
			function.setParameterValue(index, newBeta.getValueAt(i, 0));
		}
	}
	
	
	/**
	 * This method returns a double that is the convergence indicator based on the gradient, the hessian and the value of the log-likelihood function.
	 * @param gradient a Matrix instance
	 * @param hessian a Matrix instance
	 * @param llk a double that is the value of the log-likelihood function
	 * @return the indicator (a Double instance)
	 */
	protected double calculateConvergence(Matrix gradient, Matrix hessian, double llk) {
		return gradient.transpose().multiply(hessian.getInverseMatrix()).multiply(gradient).scalarMultiply(1d / llk).getValueAt(0, 0);
	}



	/**
	 * This method returns the number of iterations to reach convergence. It returns -1 if 
	 * convergence could not be reached.
	 * @return an integer
	 */
	public int getNumberOfIterations() {
		if (convergenceAchieved) {
			return iterationID;
		} else {
			return -1;
		}
	}
	
}
