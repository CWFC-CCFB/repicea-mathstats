/*
 * This file is part of the repicea library.
 *
 * Copyright (C) 2009-2016 Mathieu Fortin for Rouge-Epicea
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
package repicea.simulation;

import java.io.Serializable;
import java.util.HashMap;
import java.util.Map;

import repicea.math.Matrix;
import repicea.simulation.covariateproviders.StochasticImplementation;
import repicea.stats.estimates.Estimate;

/**
 * The SensitivityAnalysisParameter class is an abstract class for the implementation of 
 * varying parameters in a sensitivity analysis.
 * @author Mathieu Fortin - April 2016
 *
 * @param <E> an Estimate-derived class
 */
@SuppressWarnings({ "serial", "rawtypes" })
public abstract class SensitivityAnalysisParameter<E extends Estimate> implements Serializable, StochasticImplementation {

	final Map<Integer, Matrix> simulatedParameters;		// refers to the realization id only
	private E parameterEstimates;
	protected boolean isParametersVariabilityEnabled;

	protected SensitivityAnalysisParameter(boolean isParametersVariabilityEnabled) {
		this.isParametersVariabilityEnabled = isParametersVariabilityEnabled;
		simulatedParameters = new HashMap<Integer, Matrix>();
	}
	
	protected void setParameterEstimates(E estimate) {
		this.parameterEstimates = estimate;
	}

	protected E getParameterEstimates() {return parameterEstimates;}

	
	/**
	 * This method calls the setSpecificParametersDeviateForThisRealization method if the parameter variability is enabled and returns 
	 * a realization-specific simulated vector of model parameters. Otherwise it returns a default vector (beta). Note that the simulated
	 * parameters are related to the Monte Carlo realization. For instance, all subject in a given Monte Carlo realization will have the
	 * same simulation parameters. 
	 * @param subject a subject that implements the MonteCarloSimulationCompliantObject interface
	 * @return a vector of parameters
	 */
	protected synchronized Matrix getParametersForThisRealization(MonteCarloSimulationCompliantObject subject) {
		if (isParametersVariabilityEnabled) {
			String subjectPlusMonteCarloId = REpiceaPredictor.getSubjectPlusMonteCarloSpecificId(subject.getSubjectId(), subject.getMonteCarloRealizationId());
			int hashCodeSubjectId = subjectPlusMonteCarloId.hashCode();
			if (!simulatedParameters.containsKey(hashCodeSubjectId)) {		// the simulated parameters remain constant within the same Monte Carlo iteration
				Matrix randomDeviates = getParameterEstimates().getRandomDeviate();
				simulatedParameters.put(hashCodeSubjectId, randomDeviates);
			}
			return simulatedParameters.get(hashCodeSubjectId);
		} else {
			return getParameterEstimates().getMean();
		}
	}

	@Override
	public boolean isStochastic() {return isParametersVariabilityEnabled;}

}
