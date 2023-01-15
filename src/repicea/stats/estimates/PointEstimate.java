/*
 * This file is part of the repicea-statistics library.
 *
 * Copyright (C) 2009-2018 Mathieu Fortin for Rouge-Epicea
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
package repicea.stats.estimates;

import java.security.InvalidParameterException;
import java.util.ArrayList;
import java.util.Collections;
import java.util.List;
import java.util.Map;
import java.util.concurrent.ConcurrentHashMap;

import repicea.math.Matrix;
import repicea.stats.distributions.GaussianDistribution;
import repicea.stats.distributions.utility.GaussianUtility;
import repicea.stats.sampling.PopulationUnit;

/**
 * An abstract class for point estimates (total or mean). 
 * @author Mathieu Fortin - March 2021
 *
 * @param <O> a PopulationEstimate-derived class
 */
@SuppressWarnings("serial")
public abstract class PointEstimate<O extends PopulationUnit> extends Estimate<GaussianDistribution> {

	private final Map<String, O> observations;
	protected int nRows;
	protected int nCols;
	private final double populationSize;

	
	/**
	 * Basic constructor without population size.
	 */
	protected PointEstimate() {
		super(new GaussianDistribution(0d, 1d));
		observations = new ConcurrentHashMap<String, O>();
		populationSize = -1d;
		estimatorType = EstimatorType.LeastSquares;
	}

	/**
	 * Constructor with population size.
	 * @param populationSize the number of units in the population.
	 */
	protected PointEstimate(double populationSize) {
		super(new GaussianDistribution(0d, 1d));
		if (populationSize <= 0) {
			throw new InvalidParameterException("The population size must be greater than 0!");
		}
		observations = new ConcurrentHashMap<String, O>();
		this.populationSize = populationSize;
		estimatorType = EstimatorType.LeastSquares;
	}
	

	/**
	 * Create a Matrix instance with each row representing one observation. The order is ensured by
	 * the list of sample Ids.
	 * @return a Matrix
	 */
	protected Matrix getObservationMatrix() {
		Matrix outputMatrix = null;
		int nbObservations = getObservations().size();
		int nbElementsPerObs = 0;
		List<String> sampleIds = getSampleIds();
		for (int i = 0; i < sampleIds.size(); i++) {
			String sampleId = sampleIds.get(i);
			O obs = getObservations().get(sampleId);
			if (outputMatrix == null) {
				nbElementsPerObs = obs.getData().m_iRows;
				outputMatrix = new Matrix(nbObservations, nbElementsPerObs);
			}
			outputMatrix.setSubMatrix(obs.getData().transpose(), i, 0);
		}
		return outputMatrix;
	}

	/**
	 * Create a List with the ordered sample ids 
	 * @return a List instance
	 */
	protected final List<String> getSampleIds() {
		List<String> sampleIds = new ArrayList<String>();
		sampleIds.addAll(observations.keySet());
		Collections.sort(sampleIds);
		return sampleIds;
	}
	
	protected int getNumberOfElementsPerObservation() {
		if (!getObservations().isEmpty()) {
			return getObservations().values().iterator().next().getData().m_iRows;
		} else {
			return -1;
		}
		
	}
	
	

	/**
	 * Add an observation to the sample.
	 * 
	 * @param sampleId a string that represents the sample id
	 * @param obs a PopulationUnitObservation instance
	 */
	public void addObservation(O obs) {
		if (obs == null) {
			throw new InvalidParameterException("The obs argument must be non null!");
		}
		String sampleId = obs.getSampleId();
		if (observations.containsKey(sampleId)) {
			throw new InvalidParameterException("The sample id " + sampleId + " is already contained in the observation map!");
		}
		if (nCols == 0) {
			nCols = obs.getData().m_iCols;
		}
		if (nRows == 0) {
			nRows = obs.getData().m_iRows;
		}
		if (obs.getData().m_iCols != nCols || obs.getData().m_iRows != nRows) {
			throw new InvalidParameterException("The observation is incompatible with what was already observed!");
		} else {
			observations.put(sampleId, obs);
		}
	}

	@SuppressWarnings("rawtypes")
	@Override
	protected boolean isMergeableEstimate(Estimate<?> estimate) {
		if (estimate.getClass().equals(getClass())) {
			PointEstimate pe = (PointEstimate) estimate;
			if (getSampleIds().equals(pe.getSampleIds())) {	// make sure we have the same sample ids
				if (nRows == pe.nRows) {
					if (nCols == pe.nCols) {
						return true;
					}
				}
			}
		}
		return false;
	}
	
	protected Map<String, O> getObservations() {return observations;}

	public boolean isPopulationSizeKnown() {return populationSize != -1;}
	
	public double getPopulationSize() {return populationSize;}
	
	protected Matrix getQuantileForProbability(double probability) {
		Matrix stdDev = getVariance().diagonalVector().elementWisePower(.5); 
		double quantile = GaussianUtility.getQuantile(probability);
		return getMean().add(stdDev.scalarMultiply(quantile));
	}
	
	@Override
	public ConfidenceInterval getConfidenceIntervalBounds(double oneMinusAlpha) {
		Matrix lowerBoundValue = getQuantileForProbability(.5 * (1d - oneMinusAlpha));
		Matrix upperBoundValue = getQuantileForProbability(1d - .5 * (1d - oneMinusAlpha));
		return new ConfidenceInterval(lowerBoundValue, upperBoundValue, oneMinusAlpha);
	}

	
	protected abstract PointEstimate<?> add(PointEstimate<?> pointEstimate);

	protected abstract PointEstimate<?> subtract(PointEstimate<?> pointEstimate);

	protected abstract PointEstimate<?> multiply(double scalar);


}