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

import java.util.List;
import java.util.concurrent.CopyOnWriteArrayList;

import repicea.math.MathematicalFunction;
import repicea.math.Matrix;
import repicea.math.SymmetricMatrix;

/**
 * This class is the main class for all optimizers.
 * @author Mathieu Fortin - August 2011
 */
public abstract class AbstractOptimizer {

	
	public static enum LineSearchMethod {
		TEN_EQUAL(10),	
		SINGLE_TRIAL(0),
		HALF_STEP(20);
		private Number numberOfTrial;

		private LineSearchMethod(int numberOfTrial) {
			this.numberOfTrial = numberOfTrial;
		}

		public Number getParameter() {return numberOfTrial;}
	}

	/**
	 * The OptimizationException class encompasses all the exception that can be thrown when the
	 * optimizer fails to reach convergence.
	 * @author Mathieu Fortin - June 2011
	 */
	public static class OptimizationException extends Exception {
		private static final long serialVersionUID = 20110614L;

		public OptimizationException(String message) {
			super(message);
		}
	}
	
	private final CopyOnWriteArrayList<OptimizerListener> listeners;

	protected double convergenceCriterion;
	protected boolean convergenceAchieved;
	protected double optimalValue;
	
	protected Matrix betaVector;
	protected SymmetricMatrix hessianMatrix;

	
	/**
	 * Default constructor.
	 */
	public AbstractOptimizer() {
		listeners = new CopyOnWriteArrayList<OptimizerListener>();
	}
	
	protected void addOptimizerListener(OptimizerListener listener) {
		if (!listeners.contains(listener)) {
			listeners.add(listener);
		}
	}

	protected void removeOptimizerListener(OptimizerListener listener) {
		listeners.remove(listener);
	}
	
	public boolean isConvergenceAchieved() {return convergenceAchieved;}

	public Matrix getParametersAtMaximum() {return betaVector;}

	public SymmetricMatrix getHessianAtMaximum() {return hessianMatrix;}

	
	@Override
	public String toString() {
		return this.getClass().getSimpleName();
	}
	
	/**
	 * This method returns the optimal value resulting from the optimization.
	 * @return a double
	 */
	public double getOptimalValue() {
		if (convergenceAchieved) {
			return optimalValue;
		} else {
			return Double.NaN;
		}
	}

	public double getConvergenceCriterion() {return convergenceCriterion;}
	
	public void setConvergenceCriterion(double conv) {convergenceCriterion = conv;}
	
	/**
	 * This method optimizes the likelihood of the StatisticalModel instance.
	 * @param function a AbstractMathematicalFunction object
	 * @param indicesOfParametersToOptimize	a List instance that contains the indices of the parameters to be optimized
	 */
	public abstract boolean optimize(MathematicalFunction function, List<Integer> indicesOfParametersToOptimize) throws OptimizationException;

	protected void fireOptimizerEvent(String actionString) {
		for (OptimizerListener listener : listeners) {
			listener.optimizerDidThis(actionString);
		}
	}
	
}
