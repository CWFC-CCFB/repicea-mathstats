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
package repicea.stats.sampling;

import java.security.InvalidParameterException;

import repicea.math.Matrix;
import repicea.math.SymmetricMatrix;
import repicea.stats.distributions.EmpiricalDistribution;

/**
 * An implementation of the estimator of the mean.<p>
 * The population is assumed to be infinite.
 * @author Mathieu Fortin - April 2016
 */
@SuppressWarnings("serial")
public class PopulationMeanEstimate extends AbstractSimplePointEstimate {
		
	protected final EmpiricalDistribution sample;
	protected Matrix mean;
	protected SymmetricMatrix variance;
	
	/**
	 * Basic constructor.
	 */
	public PopulationMeanEstimate() {
		super();
		sample = new EmpiricalDistribution();
	}
	
	protected void recalculate() {
		mean = sample.getMean();
		variance = sample.getVariance().scalarMultiply(1d/getSampleSize());// * finitePopulationCorrectionFactor);
		getDistribution().setMean(mean);		// the mean and variance and not tied to the the distribution
		getDistribution().setVariance(variance);	// consequently, they have to be specified before drawing the random deviates
		needsToBeRecalculated = false;
	}
	
	@Override
	protected final Matrix getMeanFromDistribution() {
		if (needsToBeRecalculated) {
			recalculate();
		}
		return mean;
	}
	
	@Override
	protected final SymmetricMatrix getVarianceFromDistribution() {
		if (needsToBeRecalculated) {
			recalculate();
		}
		return variance;
	}

	@Override
	public void addObservation(Matrix obs, String obsId, String stratumName) {
		super.addObservation(obs, obsId, stratumName);
		sample.addRealization(obs);
	}

	@Override
	protected final PopulationMeanEstimate add(PointEstimate pointEstimate) {
		if (isMergeableEstimate(pointEstimate)) {
			PopulationMeanEstimate newEstimate = getEmptyEstimate();
			PopulationMeanEstimate meanEstimate = (PopulationMeanEstimate) pointEstimate;
			for (String obsId : getPopulationUnitIds()) {
				Matrix thisUnit = getObservations().get(obsId);
				Matrix thatUnit = meanEstimate.getObservations().get(obsId);
				newEstimate.addObservation(thisUnit.add(thatUnit), obsId);
			}
			return newEstimate;
		} else {
			throw new InvalidParameterException("Incompatible point estimates!");
		}
	}

	@Override
	protected final PopulationMeanEstimate subtract(PointEstimate pointEstimate) {
		if (isMergeableEstimate(pointEstimate)) {
			PopulationMeanEstimate newEstimate = getEmptyEstimate();
			PopulationMeanEstimate meanEstimate = (PopulationMeanEstimate) pointEstimate;
			for (String obsId : getPopulationUnitIds()) {
				Matrix thisUnit = getObservations().get(obsId);
				Matrix thatUnit = meanEstimate.getObservations().get(obsId);
				newEstimate.addObservation(thisUnit.subtract(thatUnit), obsId);
			}
			return newEstimate;
		} else {
			throw new InvalidParameterException("Incompatible point estimates!");
		}
	}

	@Override
	protected final PopulationMeanEstimate multiply(double scalar) {
		PopulationMeanEstimate newEstimate = getEmptyEstimate();
		for (String obsId : getPopulationUnitIds()) {
			Matrix thisUnit = getObservations().get(obsId);
			addObservation(thisUnit.scalarMultiply(scalar), obsId);
		}
		return newEstimate;
	}

	@Override
	PopulationMeanEstimate getEmptyEstimate() {
		return new PopulationMeanEstimate();
	}

}
