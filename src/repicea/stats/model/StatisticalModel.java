/*
 * This file is part of the repicea library.
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
package repicea.stats.model;

import repicea.math.Matrix;
import repicea.stats.estimators.Estimator;

/**
 * This interface defines the services provided by a statistical model.
 * @author Mathieu Fortin - October 2011
 */
public interface StatisticalModel { 

	/**
	 * Return the parameter estimates produced by the estimator.
	 * @return a Matrix instance
	 */
	public default Matrix getParameters() {
		if (getEstimator().isConvergenceAchieved()) {
			return getEstimator().getParameterEstimates().getMean();
		} else {
			throw new UnsupportedOperationException("The model parameters have not been estimated yet!");
		}
	}

	/**
	 * This method returns the results of the fit on screen.
	 * @return a string containing the summary of the model fit
	 */
	public default String getSummary() {
		StringBuilder sb = new StringBuilder();
		sb.append(toString() + System.lineSeparator());
		sb.append("Model definition : " + getModelDefinition() + System.lineSeparator() + System.lineSeparator());
		sb.append(getEstimator().getReport() + System.lineSeparator());
		return sb.toString();
	}

	
	/**
	 * This method computes the parameter estimates.
	 */
	public void doEstimation();
	
	/**
	 * This method returns the model definition as entered by the user.
	 * @return a String
	 */
	public String getModelDefinition();
	
	
	/**
	 * This method returns the optimizer of the log-likelihood function.
	 * @return an Optimizer instance
	 */
	public Estimator getEstimator();
	
}
