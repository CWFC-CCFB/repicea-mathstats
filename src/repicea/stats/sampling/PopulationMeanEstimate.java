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

/**
 * An implementation of the estimator of the mean.<p>
 * The population is assumed to be infinite. The variance is 
 * that of the estimate (ie. s^2 / (n * (n-1))) and not the
 * residual variance.
 * @author Mathieu Fortin - April 2016
 */
@SuppressWarnings("serial")
public class PopulationMeanEstimate extends AbstractSimplePointEstimate {
		
//	protected final EmpiricalDistribution sample;
	protected Matrix mean;
	protected SymmetricMatrix variance;
	
	/**
	 * Basic constructor.
	 */
	public PopulationMeanEstimate() {
		super();
//		sample = new EmpiricalDistribution();
	}
	
	protected void recalculate() {
		mean = computeMeanInternally();
		variance = computeVarianceInternally();
		getDistribution().setMean(mean);		// the mean and variance and not tied to the the distribution
		getDistribution().setVariance(variance);	// consequently, they have to be specified before drawing the random deviates
		needsToBeRecalculated = false;
	}
	
	
	protected Matrix computeMeanInternally() {
		if (observations.size() > 0) {
			Matrix mean = new Matrix(nRows, nCols);
			for (Matrix m : observations.values()) {
				mean = mean.add(m);
			}
			return mean.scalarMultiply(1d / getSampleSize());
		} else {
			return null;
		}
	}

	protected SymmetricMatrix computeVarianceInternally() {
		if (observations.size() > 1) {
			Matrix sse = null;
			Matrix error;
			for (Matrix mat : observations.values()) {
				error = mat.subtract(mean);
				if (sse == null) {
					sse = error.multiply(error.transpose());
				} else {
					sse = sse.add(error.multiply(error.transpose()));
				}
			}
			int n = observations.size();
			sse = sse.scalarMultiply(1d / (n * (n-1)));
			return  SymmetricMatrix.convertToSymmetricIfPossible(sse);
		}
		return null;
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

//	@Override
//	public void addObservation(Matrix obs, String obsId, String stratumName) {
//		super.addObservation(obs, obsId, stratumName);
////		sample.addRealization(obs);
//	}

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
			newEstimate.addObservation(thisUnit.scalarMultiply(scalar), obsId);
		}
		return newEstimate;
	}

	@Override
	PopulationMeanEstimate getEmptyEstimate() {
		return new PopulationMeanEstimate();
	}

}
