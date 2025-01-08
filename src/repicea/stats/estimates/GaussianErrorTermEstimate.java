/*
 * This file is part of the repicea-statistics library.
 *
 * Copyright (C) 2009-2014 Mathieu Fortin for Rouge-Epicea
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

import repicea.math.Matrix;
import repicea.math.SymmetricMatrix;
import repicea.math.utility.GaussianUtility;
import repicea.stats.StatisticalUtility.TypeMatrixR;
import repicea.stats.distributions.CenteredGaussianDistribution;
import repicea.stats.distributions.GaussianErrorTermList;

/**
 * The GaussianErrorTermEstimate class handles the complex covariance structure in linear and nonlinear models. It assumes the mean of the
 * distribution is 0 and the variance-covariance matrix adapts to the size of the GaussianErrorTermList instance that is passed to the object 
 * through the setErrorTermList method. IMPORTANT: this method should be called before any other regular method. Otherwise, a null object is 
 * returned. The class is not thread safe either.
 * @author Mathieu Fortin - August 2014
 */
@SuppressWarnings("serial")
public final class GaussianErrorTermEstimate extends AbstractEstimate<Matrix, SymmetricMatrix, CenteredGaussianDistribution> {

	/**
	 * General constructor.
	 * @param variance a double
	 * @param correlationParameter a double 
	 * @param type a TypeMatrixR enum
	 */
	public GaussianErrorTermEstimate(SymmetricMatrix variance, double correlationParameter, TypeMatrixR type) {
		super(new CenteredGaussianDistribution(variance, correlationParameter, type));
		estimatorType = EstimatorType.LikelihoodBased;
	}

	/**
	 * Constructor for univariate distribution.
	 * @param variance a double
	 */
	public GaussianErrorTermEstimate(SymmetricMatrix variance) {
		this(variance, 0, null);
	}

	public Matrix getMean(GaussianErrorTermList errorTermList) {
		return getDistribution().getMean(errorTermList);
	}
	
	public Matrix getVariance(GaussianErrorTermList errorTermList) {
		return getDistribution().getVariance(errorTermList);
	}
	
	public Matrix getRandomDeviate(GaussianErrorTermList errorTermList) {
		return getDistribution().getRandomRealization(errorTermList);
	}
	
	
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

}
