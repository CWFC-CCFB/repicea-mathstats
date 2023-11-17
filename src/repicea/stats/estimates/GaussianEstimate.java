/*
 * This file is part of the repicea-statistics library.
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
package repicea.stats.estimates;

import java.io.Serializable;

import repicea.math.Matrix;
import repicea.math.SymmetricMatrix;
import repicea.math.utility.GaussianUtility;
import repicea.stats.CentralMomentsSettable;
import repicea.stats.distributions.GaussianDistribution;

/**
 * This class contains the elements related to the random effects, i.e. the best linear unbiased predictors (blups) as well as their variances.
 * @author Mathieu Fortin - October 2011
 */
public class GaussianEstimate extends Estimate<Matrix, GaussianDistribution> implements CentralMomentsSettable, Serializable {

	private static final long serialVersionUID = 20120725L;
	
	/**
	 * Common constructor. By default the Gaussian distribution that supports this estimate has a mean 0 and a variance 1.
	 */
	public GaussianEstimate() {
		super(new GaussianDistribution());
		estimatorType = EstimatorType.LikelihoodBased;
	}
	
	/**
	 * Constructor with the mean and variance.
	 * @param mean a Matrix instance that contains the mean 
	 * @param variance a SymmetricMatrix instance that contains the variance-covariance
	 */
	public GaussianEstimate(Matrix mean, SymmetricMatrix variance) {
		this();
		setMean(mean);
		setVariance(variance);
	}

	/**
	 * Constructor for univariate distribution.
	 * @param mean a double that stands for the mean
	 * @param variance a double that stands for the variance
	 */
	public GaussianEstimate(double mean, double variance) {
		this();
		Matrix meanMat = new Matrix(1,1);
		meanMat.setValueAt(0, 0, mean);
		SymmetricMatrix varianceMat = new SymmetricMatrix(1);
		varianceMat.setValueAt(0, 0, variance);
		setMean(meanMat);
		setVariance(varianceMat);
	}
	
	

	@Override
	public void setVariance(SymmetricMatrix variance) {
		getDistribution().setVariance(variance);
	}

	@Override
	public void setMean(Matrix mean) {
		getDistribution().setMean(mean);
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
	
