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

import java.util.logging.Level;

import repicea.stats.data.StatisticalDataException;
import repicea.stats.estimators.Estimator;
import repicea.stats.estimators.Estimator.EstimatorException;
import repicea.util.REpiceaLogManager;

/**
 * The AbstractStatisticalModel class implements the StatisticalModel interface. It contains the
 * basic features for a StatisticalModel, namely a data structure, a log-likelihood function and an
 * optimizer.
 * @author Mathieu Fortin - August 2011
 */
public abstract class AbstractStatisticalModel implements StatisticalModel {

	public static String LOGGER_NAME = "AbstractStatisticalModel";
	
	protected Estimator estimator;
	
	private Object optimizerParameters; 
	
	private String modelDefinition;

	/**
	 * Default constructor.
	 */
	protected AbstractStatisticalModel() {}
	
	/**
	 * Set the Estimator for the model.
	 * @param estimator an Estimator instance
	 */
	public void setEstimator(Estimator estimator) {this.estimator = estimator;}
	
	@Override
	public Estimator getEstimator() {
		if (estimator == null) {
			setEstimator(instantiateDefaultEstimator());
		}
		return estimator;
	}

	
	/**
	 * This method defines the default optimizer which is to be specific to the derived classes.
	 * @return an Optimizer instance 
	 */
	protected abstract Estimator instantiateDefaultEstimator();

	
	/**
	 * Set the parameter for the optimizer.
	 * @param optimizerParameters an object containing the parameters
	 */
	public void setOptimizerParameters(Object optimizerParameters) {this.optimizerParameters = optimizerParameters;}
	
	protected Object getOptimizerParameters() {return optimizerParameters;}
	
	@Override
	public void doEstimation() {
		REpiceaLogManager.logMessage(LOGGER_NAME, Level.FINE, LOGGER_NAME, "Optimization using " + getEstimator().toString() + ".");
		try {
			if (getEstimator().doEstimation()) {
				REpiceaLogManager.logMessage(LOGGER_NAME, Level.FINE, LOGGER_NAME,"Convergence achieved!");
			} else {
				REpiceaLogManager.logMessage(LOGGER_NAME, Level.WARNING, LOGGER_NAME, "Unable to reach convergence.");
			}
		} catch (EstimatorException e) {
			REpiceaLogManager.logMessage(LOGGER_NAME, Level.SEVERE, LOGGER_NAME,"An error occured while optimizing the log likelihood function.");
			e.printStackTrace();
		}
	}


	@Override
	public String getModelDefinition() {return modelDefinition;}

	/**
	 * This method sets the model definition and computes the appropriate matrix from the data.
	 * @param modelDefinition a String that defines the model
	 * @param additionalParm an object containing additional parameters if need be
	 * @throws StatisticalDataException if the model definition cannot be parsed
	 */
	protected void setModelDefinition(String modelDefinition, Object additionalParm) throws StatisticalDataException {
		this.modelDefinition = modelDefinition;
	}
	
	protected void setModelDefinition(String modelDefinition) throws StatisticalDataException {
		setModelDefinition(modelDefinition, null);
	}

}
