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
package repicea.stats.estimators;

import java.security.InvalidParameterException;
import java.text.DecimalFormat;
import java.text.NumberFormat;
import java.util.ArrayList;
import java.util.Collections;
import java.util.List;
import java.util.logging.Level;

import repicea.math.Matrix;
import repicea.math.SymmetricMatrix;
import repicea.math.optimizer.AbstractOptimizer.LineSearchMethod;
import repicea.math.optimizer.AbstractOptimizer.OptimizationException;
import repicea.stats.data.DataSet;
import repicea.stats.estimates.Estimate;
import repicea.stats.estimates.GaussianEstimate;
import repicea.stats.estimators.AbstractEstimator.EstimatorCompatibleModel;
import repicea.stats.estimators.MaximumLikelihoodEstimator.MaximumLikelihoodCompatibleModel;
import repicea.stats.model.CompositeLogLikelihood;
import repicea.util.REpiceaLogManager;

/**
 * Implement a maximum likelihood estimator based on 
 * the Newton-Raphson algorithm.
 * @author Mathieu Fortin - August 2011
 */
public class MaximumLikelihoodEstimator extends AbstractEstimator<MaximumLikelihoodCompatibleModel> {

	
	public interface MaximumLikelihoodCompatibleModel extends EstimatorCompatibleModel {
	
		public double getConvergenceCriterion();

		/**
		 * Return the model log-likelihood function.
		 * @return a CompositeLogLikelihood instance
		 */
		public CompositeLogLikelihood getCompleteLogLikelihood();
		
		/**
		 * Set the parameters in the log-likelihood function.
		 * @param beta a Matrix instance
		 */
		public default void setParameters(Matrix beta) {
			for (int i = 0; i < beta.m_iRows; i++) {
				getCompleteLogLikelihood().setParameterValue(i, beta.getValueAt(i, 0));
			}
		}

		
		
	}
	
	
	protected static class LikelihoodValue implements Comparable<LikelihoodValue> {

		private double llk;
		private Matrix beta;
		
		protected LikelihoodValue(Matrix beta, double llk) {
			this.beta = beta.getDeepClone();
			this.llk = llk;
		}
		
		@Override
		public int compareTo(LikelihoodValue arg0) {
			double reference = ((LikelihoodValue) arg0).llk;
			if (this.llk < reference) {
				return 1;
			} else if (this.llk == reference) {
				return 0;
			} else {
				return -1;
			}
		}
		
		protected Matrix getParameters() {return beta;}
	}

	public static String LOGGER_NAME = "MLEstimator";
	protected GaussianEstimate parameterEstimate;
	protected final repicea.math.optimizer.NewtonRaphsonOptimizer nro;
	
	/**
	 * Constructor. 
	 */
	public MaximumLikelihoodEstimator(MaximumLikelihoodCompatibleModel model) {
		super(model);
		nro = new repicea.math.optimizer.NewtonRaphsonOptimizer();
	}
	
	
	/**
	 * Return the model parameters. <br>
	 * <br>
	 * A new Matrix instance is created each time this method is called.
	 * @return a Matrix instance
	 */
	private Matrix getParameters(MaximumLikelihoodCompatibleModel model) {
		int nbParms = model.getCompleteLogLikelihood().getNumberOfParameters();
		Matrix beta = new Matrix(nbParms, 1);
		for (int i = 0; i < nbParms; i++) {
			beta.setValueAt(i, 0, model.getCompleteLogLikelihood().getParameterValue(i));
		}
		return beta;
	}
	
	
	/**
	 * This method scans the log likelihood function within a range of values for a particular parameter.
	 * @param parameterName the index of the parameter
	 * @param start the starting value
	 * @param end the ending value
	 * @param step the step between these two values.
	 */
	public void gridSearch(int parameterName, double start, double end, double step) {
		REpiceaLogManager.logMessage(MaximumLikelihoodEstimator.LOGGER_NAME, Level.FINER, MaximumLikelihoodEstimator.LOGGER_NAME, "Initializing grid search...");
		ArrayList<LikelihoodValue> likelihoodValues = new ArrayList<LikelihoodValue>();
		Matrix originalParameters = getParameters(model);
		double llk;
		for (double value = start; value < end + step; value+=step) {
			Matrix beta = originalParameters.getDeepClone();
			beta.setValueAt(parameterName, 0, value);
			model.setParameters(beta);
			model.getCompleteLogLikelihood().reset();
			llk = model.getCompleteLogLikelihood().getValue();
			likelihoodValues.add(new LikelihoodValue(beta, llk));
			REpiceaLogManager.logMessage(MaximumLikelihoodEstimator.LOGGER_NAME, Level.FINER, MaximumLikelihoodEstimator.LOGGER_NAME, "Parameters : " + beta.toString() + "; Log-likelihood : " + llk);
		}
		
		Collections.sort(likelihoodValues);
		LikelihoodValue lk;
		Matrix bestFittingParameters = null;
		for (int i = 0; i < likelihoodValues.size(); i++) {
			lk = likelihoodValues.get(i);
			if (!Double.isNaN(lk.llk)) {
				bestFittingParameters = lk.getParameters();
				break;
			}
		}
		if (bestFittingParameters == null) {
			throw new InvalidParameterException("All the likelihoods of the grid are NaN!");
		} else {
			model.setParameters(bestFittingParameters);
		}
	}

	
	@Override
	public boolean doEstimation() throws EstimatorException {
		nro.setConvergenceCriterion(model.getConvergenceCriterion());
		
		CompositeLogLikelihood llk = model.getCompleteLogLikelihood();
		List<Integer> indices = new ArrayList<Integer>();
		int nbParameters = model.getCompleteLogLikelihood().getNumberOfParameters();
		for (int i = 0; i < nbParameters; i++) {
			indices.add(i);
		}
		try {
			REpiceaLogManager.logMessage(LOGGER_NAME, Level.INFO, LOGGER_NAME, "Starting optimization");
			nro.optimize(llk, indices);
		} catch (OptimizationException e) {
			REpiceaLogManager.logMessage(LOGGER_NAME, Level.SEVERE, LOGGER_NAME, e.getMessage());
			parameterEstimate = null;
			return false;
		}
		if (nro.isConvergenceAchieved()) {
			SymmetricMatrix varCov = nro.getHessianAtMaximum().getInverseMatrix().scalarMultiply(-1d);
			parameterEstimate = new GaussianEstimate(nro.getParametersAtMaximum(), varCov);
			return true;
		} else {
			parameterEstimate = null;
			return false;
		}
		
	}

	/**
	 * Returns the maximum log likelihood value after convergence.
	 * @return a double
	 */
	public double getMaximumLogLikelihood() {return nro.getOptimalValue();}


	@Override
	public boolean isConvergenceAchieved() {return nro.isConvergenceAchieved();}

	@Override
	public Estimate<?> getParameterEstimates() {
		return parameterEstimate;
	}
	
	@Override
	public String toString() {return "Maximum likelihood estimator";}

	/**
	 * Sets the line search method. <br>
	 * <br>
	 *  
	 * If the lsm parameter is null,
	 * the line search method is set to LineSearchMethod.TEN_EQUAL 
	 * by default.
	 * @param lsm a LineSearchMethod enum
	 */
	public void setLineSearchMethod(LineSearchMethod lineSearchMethod) {
		nro.setLineSearchMethod(lineSearchMethod);
	}
	
	@Override
	public DataSet getConvergenceStatusReport() {
		List<String> fieldNames = new ArrayList<String>();
		fieldNames.add("Element");
		fieldNames.add("Value");
		DataSet dataSet = new DataSet(fieldNames);
		Object[] record = new Object[2];
		record[0] = "Converged";
		record[1] = isConvergenceAchieved();
		dataSet.addObservation(record);
		record[0] = "Maximum log-likelihood";
		record[1] = getMaximumLogLikelihood();
		dataSet.addObservation(record);
		record[0] = "AIC";
		record[1] = - 2 * getMaximumLogLikelihood() + 2 * nro.getParametersAtMaximum().getNumberOfElements();
		dataSet.addObservation(record);
		record[0] = "BIC";
		record[1] = - 2 * getMaximumLogLikelihood() + nro.getParametersAtMaximum().getNumberOfElements() * Math.log(model.getNumberOfObservations());
		dataSet.addObservation(record);
		return dataSet;
	}


	@Override
	public String getReport() {
		if (!isConvergenceAchieved()) {
			return "The log-likelihood function has not been or cannot be optimized.";
		} else {
			StringBuilder sb = new StringBuilder();
			DataSet convergenceDataset = getConvergenceStatusReport();
			DecimalFormat decFormat = new DecimalFormat();
			decFormat.setMaximumFractionDigits(4);
			decFormat.setMinimumFractionDigits(4);
			convergenceDataset.setFormatter(1, decFormat);
			sb.append(convergenceDataset.toString() + System.lineSeparator());
			DataSet parameterDataset = getParameterEstimatesReport();
			decFormat = new DecimalFormat();
			decFormat.setMaximumFractionDigits(6);
			decFormat.setMinimumFractionDigits(6);
			parameterDataset.setFormatter(1, decFormat);
			parameterDataset.setFormatter(2, decFormat);
			parameterDataset.setFormatter(4, decFormat);
			decFormat = new DecimalFormat();
			decFormat.setMaximumFractionDigits(3);
			decFormat.setMinimumFractionDigits(3);
			parameterDataset.setFormatter(3, decFormat);
			sb.append(parameterDataset.toString() + System.lineSeparator());
			return sb.toString();
		}

	}
}
