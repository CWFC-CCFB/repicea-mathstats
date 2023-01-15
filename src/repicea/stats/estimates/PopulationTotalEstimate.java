/*
 * This file is part of the repicea-statistics library.
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
package repicea.stats.estimates;

import java.security.InvalidParameterException;
import java.util.List;

import repicea.math.Matrix;
import repicea.math.SymmetricMatrix;
import repicea.stats.sampling.PopulationUnitWithUnequalInclusionProbability;

/**
 * This class implements the classical Horvitz-Thompson estimator of the total (tau) in a
 * context of random sampling WITHOUT replacement. 
 * @author Mathieu Fortin - September 2016
 */
@SuppressWarnings("serial")
public class PopulationTotalEstimate extends PointEstimate<PopulationUnitWithUnequalInclusionProbability> {

	/**
	 * Constructor.
	 */
	public PopulationTotalEstimate() {
		super();
	}

	/**
	 * Constructor with population size.
	 * @param populationSize the number of units in the population
	 */
	public PopulationTotalEstimate(double populationSize) {
		super(populationSize);
	}

	@Override
	protected Matrix getMeanFromDistribution() {
		Matrix total = new Matrix(nRows, nCols);
		int sampleSize = getObservations().size();
		for (PopulationUnitWithUnequalInclusionProbability observation : getObservations().values()) {
			total = total.add(observation.getData().scalarMultiply(1d/(sampleSize * observation.getInclusionProbability())));
		}
		return total;
	}
	
	
	/**
	 * This method returns the variance of the tau estimate.
	 * @return a Matrix
	 */
	@Override
	protected SymmetricMatrix getVarianceFromDistribution() {
		int n = getObservations().size();
		PopulationUnitWithUnequalInclusionProbability obs_i;
		PopulationUnitWithUnequalInclusionProbability obs_j;
		double pi_i;
		double pi_j;
		double pi_ij;
		Matrix varianceContribution;
		Matrix variance = null;
		List<String> sampleIds = getSampleIds();
		for (int i = 0; i < getObservations().size(); i++) {
			for (int j = i; j < getObservations().size(); j++) {
				obs_i = getObservations().get(sampleIds.get(i));
				obs_j = getObservations().get(sampleIds.get(j));
				pi_i = n * obs_i.getInclusionProbability();
				pi_j = n * obs_j.getInclusionProbability();
				if (i == j) {
					varianceContribution = obs_i.getData().multiply(obs_i.getData().transpose()).scalarMultiply((1 - pi_i)/(pi_i*pi_i));
				} else {
					pi_ij = pi_i * (n-1) * obs_j.getInclusionProbability() / (1 - obs_i.getInclusionProbability());
					double factor = (pi_ij - pi_i * pi_j)/(pi_i * pi_j * pi_ij);
					varianceContribution = obs_i.getData().multiply(obs_j.getData().transpose()).scalarMultiply(2 * factor);
				}
				if (variance == null) {
					variance = varianceContribution;
				} else {
					variance = variance.add(varianceContribution);
				}
			}
		}
		return SymmetricMatrix.convertToSymmetricIfPossible(variance);
	}

	
	@Override
	protected boolean isMergeableEstimate(Estimate<?> estimate) {
		boolean isMergeable = super.isMergeableEstimate(estimate);
		if (isMergeable) {
			PopulationTotalEstimate est = (PopulationTotalEstimate) estimate;
			for (String sampleId : getSampleIds()) {
				PopulationUnitWithUnequalInclusionProbability thisUnit = getObservations().get(sampleId);
				PopulationUnitWithUnequalInclusionProbability thatUnit = est.getObservations().get(sampleId);
				if (thisUnit.getInclusionProbability() != thatUnit.getInclusionProbability()) {
					return false;
				}
			}
		}
		return isMergeable;
	}

	
	@Override
	protected PopulationTotalEstimate add(PointEstimate<?> pointEstimate) {
		if (isMergeableEstimate(pointEstimate)) {
			PopulationTotalEstimate newEstimate = new PopulationTotalEstimate();
			PopulationTotalEstimate totalEstimate = (PopulationTotalEstimate) pointEstimate;
			for (String sampleId : getSampleIds()) {
				PopulationUnitWithUnequalInclusionProbability thisUnit = getObservations().get(sampleId);
				PopulationUnitWithUnequalInclusionProbability thatUnit = totalEstimate.getObservations().get(sampleId);
				PopulationUnitWithUnequalInclusionProbability newUnit = new PopulationUnitWithUnequalInclusionProbability(sampleId, 
						thisUnit.getData().add(thatUnit.getData()), 
						thisUnit.getInclusionProbability());
				newEstimate.addObservation(newUnit);
			}
			return newEstimate;
		} else {
			throw new InvalidParameterException("Incompatible point estimates!");
		}
	}

	@Override
	protected PopulationTotalEstimate subtract(PointEstimate<?> pointEstimate) {
		if (isMergeableEstimate(pointEstimate)) {
			PopulationTotalEstimate newEstimate = new PopulationTotalEstimate();
			PopulationTotalEstimate totalEstimate = (PopulationTotalEstimate) pointEstimate;
			for (String sampleId : getSampleIds()) {
				PopulationUnitWithUnequalInclusionProbability thisUnit = getObservations().get(sampleId);
				PopulationUnitWithUnequalInclusionProbability thatUnit = totalEstimate.getObservations().get(sampleId);
				PopulationUnitWithUnequalInclusionProbability newUnit = new PopulationUnitWithUnequalInclusionProbability(sampleId, 
						thisUnit.getData().subtract(thatUnit.getData()), 
						thisUnit.getInclusionProbability());
				newEstimate.addObservation(newUnit);
						
			}
			return newEstimate;
		} else {
			throw new InvalidParameterException("Incompatible point estimates!");
		}
	}

	@Override
	protected PopulationTotalEstimate multiply(double scalar) {
		PopulationTotalEstimate newEstimate = new PopulationTotalEstimate();
		for (String sampleId : getSampleIds()) {
			PopulationUnitWithUnequalInclusionProbability thisUnit = getObservations().get(sampleId);
			PopulationUnitWithUnequalInclusionProbability newUnit =	new PopulationUnitWithUnequalInclusionProbability(sampleId, 
					thisUnit.getData().scalarMultiply(scalar), 
					thisUnit.getInclusionProbability());
			newEstimate.addObservation(newUnit);
		}
		return newEstimate;
	}

}
