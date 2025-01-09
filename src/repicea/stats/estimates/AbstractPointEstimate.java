/*
 * This file is part of the repicea-mathstats library.
 *
 * Copyright (C) 2009-2018 Mathieu Fortin for Rouge-Epicea
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

import java.util.List;
import java.util.Map;

import repicea.math.Matrix;
import repicea.math.SymmetricMatrix;
import repicea.math.utility.GaussianUtility;
import repicea.stats.distributions.GaussianDistribution;
import repicea.stats.sampling.PopulationUnit;

/**
 * An abstract class for point estimates (total or mean).<p>
 * This class is more focused on sampling. 
 * @author Mathieu Fortin - March 2021, January 2025
 */
@SuppressWarnings("serial")
public abstract class AbstractPointEstimate extends AbstractEstimate<Matrix, SymmetricMatrix, GaussianDistribution> 
											implements PointEstimate {

	protected int nRows;
	protected int nCols;

	
	/**
	 * Basic constructor.
	 */
	protected AbstractPointEstimate() {
		super(new GaussianDistribution(0d, 1d));
		estimatorType = EstimatorType.LeastSquares;
	}

	@Override
	protected abstract boolean isMergeableEstimate(Estimate<?,?,?> estimate);	
	
	protected final Matrix getQuantileForProbability(double probability) {
		Matrix stdDev = getVariance().diagonalVector().elementWisePower(.5); 
		double quantile = GaussianUtility.getQuantile(probability);
		return getMean().add(stdDev.scalarMultiply(quantile));
	}
	
	@Override
	public final ConfidenceInterval getConfidenceIntervalBounds(double oneMinusAlpha) {
		Matrix lowerBoundValue = getQuantileForProbability(.5 * (1d - oneMinusAlpha));
		Matrix upperBoundValue = getQuantileForProbability(1d - .5 * (1d - oneMinusAlpha));
		return new ConfidenceInterval(lowerBoundValue, upperBoundValue, oneMinusAlpha);
	}

	/**
	 * Map the observations with their id as key and the population units
	 * as values.
	 * @return a Map instance
	 */
	protected abstract Map<String, PopulationUnit> getObservations();

	protected final int getNumberOfElementsPerObservation() {return nRows;}

	protected abstract AbstractPointEstimate add(PointEstimate pointEstimate);

	protected abstract AbstractPointEstimate subtract(PointEstimate pointEstimate);

	protected abstract AbstractPointEstimate multiply(double scalar);

	/**
	 * Provide a list of the population unit ids.
	 * @return a List of String
	 */
	protected abstract List<String> getPopulationUnitIds();
	
}