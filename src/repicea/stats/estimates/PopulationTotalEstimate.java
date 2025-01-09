/*
 * This file is part of the repicea-statistics library.
 *
 * Copyright (C) 2009-2016 Mathieu Fortin for Rouge-Epicea
 * Copyright (C) 2025 His Majesty the King in right of Canada
 * Author: Mathieu Fortin, Canadian Forest Service
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
import repicea.stats.sampling.PopulationUnit;

/**
 * Implement a classical Horvitz-Thompson estimator of the total (tau). <p>
 * 
 * The estimator assumes an even inclusion probability for each population unit and
 * random sampling WITHOUT replacement. 
 * @author Mathieu Fortin - September 2016, January 2025
 */
@SuppressWarnings("serial")
public class PopulationTotalEstimate extends AbstractSimplePointEstimate {

	private final double individualInclusionProbability;
	
	/**
	 * Constructor with population size.
	 * @param populationSize the number of units in the population
	 */
	public PopulationTotalEstimate(double populationSize) {
		super(populationSize);
		individualInclusionProbability = 1d/getPopulationSize();
	}

	@Override
	protected Matrix getMeanFromDistribution() {
		Matrix total = new Matrix(nRows, nCols);
		int sampleSize = getObservations().size();
		
		for (PopulationUnit observation : getObservations().values()) {
			total = total.add(observation.getData().scalarMultiply(1d/(sampleSize * individualInclusionProbability)));
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
		PopulationUnit obs_i;
		PopulationUnit obs_j;
		double pi_i = n * individualInclusionProbability;
		double pi_j = n * individualInclusionProbability;
		double pi_ij = pi_i * (n-1d) / (getPopulationSize() - 1d);
		Matrix varianceContribution;
		Matrix variance = null;
		List<String> sampleIds = getSampleIds();
		for (int i = 0; i < getObservations().size(); i++) {
			for (int j = i; j < getObservations().size(); j++) {
				obs_i = getObservations().get(sampleIds.get(i));
				obs_j = getObservations().get(sampleIds.get(j));
				if (i == j) {
					varianceContribution = obs_i.getData().multiply(obs_i.getData().transpose()).scalarMultiply((1 - pi_i)/(pi_i*pi_i));
				} else {
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
	protected boolean isMergeableEstimate(Estimate<?,?,?> estimate) {
		boolean isMergeable = super.isMergeableEstimate(estimate);
		if (isMergeable) {
			PopulationTotalEstimate thatPopulationTotalEstimate = (PopulationTotalEstimate) estimate;
			if (this.individualInclusionProbability == thatPopulationTotalEstimate.individualInclusionProbability) {
				return true;
			}
		}
		return false;
	}

	
	@Override
	protected PopulationTotalEstimate add(PointEstimate pointEstimate) {
		if (isMergeableEstimate(pointEstimate)) {
			PopulationTotalEstimate newEstimate = new PopulationTotalEstimate(getPopulationSize());
			PopulationTotalEstimate totalEstimate = (PopulationTotalEstimate) pointEstimate;
			for (String sampleId : getSampleIds()) {
				PopulationUnit thisUnit = getObservations().get(sampleId);
				PopulationUnit thatUnit = totalEstimate.getObservations().get(sampleId);
				PopulationUnit newUnit = new PopulationUnit(sampleId, thisUnit.getData().add(thatUnit.getData()));
				newEstimate.addObservation(newUnit);
			}
			return newEstimate;
		} else {
			throw new InvalidParameterException("Incompatible point estimates!");
		}
	}

	@Override
	protected PopulationTotalEstimate subtract(PointEstimate pointEstimate) {
		if (isMergeableEstimate(pointEstimate)) {
			PopulationTotalEstimate newEstimate = new PopulationTotalEstimate(getPopulationSize());
			PopulationTotalEstimate totalEstimate = (PopulationTotalEstimate) pointEstimate;
			for (String sampleId : getSampleIds()) {
				PopulationUnit thisUnit = getObservations().get(sampleId);
				PopulationUnit thatUnit = totalEstimate.getObservations().get(sampleId);
				PopulationUnit newUnit = new PopulationUnit(sampleId, thisUnit.getData().subtract(thatUnit.getData()));
				newEstimate.addObservation(newUnit);
			}
			return newEstimate;
		} else {
			throw new InvalidParameterException("Incompatible point estimates!");
		}
	}

	@Override
	protected PopulationTotalEstimate multiply(double scalar) {
		PopulationTotalEstimate newEstimate = new PopulationTotalEstimate(getPopulationSize());
		for (String sampleId : getSampleIds()) {
			PopulationUnit thisUnit = getObservations().get(sampleId);
			PopulationUnit newUnit =	new PopulationUnit(sampleId, thisUnit.getData().scalarMultiply(scalar));
			newEstimate.addObservation(newUnit);
		}
		return newEstimate;
	}

}
