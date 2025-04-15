/*
 * This file is part of the repicea-statistics library.
 *
 * Copyright (C) 2009-2021 Mathieu Fortin for Rouge-Epicea
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
import java.security.InvalidParameterException;

import repicea.math.Matrix;
import repicea.math.SymmetricMatrix;
import repicea.math.utility.GaussianUtility;
import repicea.stats.MomentSettable;
import repicea.stats.distributions.GaussianDistribution;

/**
 * The LogNormalEstimate class implements an estimate that follows
 * a log-normal distribution. At the moment, it only implements univariate
 * distributed estimates.
 * @author Mathieu Fortin - January 2021
 */
@SuppressWarnings("serial")
public class LogNormalEstimate extends AbstractEstimate<Matrix, SymmetricMatrix, GaussianDistribution> implements MomentSettable, Serializable{

	/**
	 * Constructor.
	 * @param mean the mean of the distribution
	 * @param variance the variance of the distribution
	 * @param onLogScale true if the mean and variance are specified on the log scale
	 */
	public LogNormalEstimate(double mean, double variance, boolean onLogScale) {
		super(new GaussianDistribution(null, null));
		if (variance <= 0d) {
			throw new InvalidParameterException("The variance must be positive!");
		}
		Matrix m = new Matrix(1,1);
		SymmetricMatrix v = new SymmetricMatrix(1);
		if (onLogScale) {
			m.setValueAt(0, 0, mean);
			setMeanOnLogScale(m);
			v.setValueAt(0, 0, variance);
			setVarianceOnLogScale(v);
		} else {
			double varianceOnLogScale = Math.log(variance / (mean * mean) + 1);
			v.setValueAt(0, 0, varianceOnLogScale);
			setVarianceOnLogScale(v);
			double meanOnLogScale = Math.log(mean) - varianceOnLogScale * .5;
			m.setValueAt(0, 0, meanOnLogScale);
			setMeanOnLogScale(m);
		}
	}

	/**
	 * Constructor with mean and variance specified on the log scale.
	 * @param mean the mean of the distribution
	 * @param variance the variance of the distribution
	 */
	public LogNormalEstimate(double mean, double variance) {
		this(mean, variance, true);
	}


	/**
	 * This method does not work with the LogNormalEstimate class. The mean is set 
	 * through the constructor.
	 */
	@Override
	public void setMean(Matrix mean) {
		throw new InvalidParameterException("The setMean(Matrix) method cannot be used with the LogNormalEstimate class!");
	}
	
	private void setMeanOnLogScale(Matrix mean) {
		getDistribution().setMean(mean);
	}

	private void setVarianceOnLogScale(SymmetricMatrix variance) {
		getDistribution().setVariance(variance);
	}

	/**
	 * This method does not work with the LogNormalEstimate class. The variance is set 
	 * through the constructor.
	 */
	@Override
	public void setVariance(SymmetricMatrix variance) {
		throw new InvalidParameterException("The setVariance(Matrix) method cannot be used with the LogNormalEstimate class!");
	}
	
	private Matrix getQuantileForProbability(double probability) {
		Matrix stdDev = getDistribution().getVariance().diagonalVector().elementWisePower(.5); 
		double quantile = GaussianUtility.getQuantile(probability);
		return getDistribution().getMean().add(stdDev.scalarMultiply(quantile));
	}

	@Override
	public ConfidenceInterval getConfidenceIntervalBounds(double oneMinusAlpha) {
		Matrix lowerBoundValue = getQuantileForProbability(.5 * (1d - oneMinusAlpha)).expMatrix();
		Matrix upperBoundValue = getQuantileForProbability(1d - .5 * (1d - oneMinusAlpha)).expMatrix();
		return new ConfidenceInterval(lowerBoundValue, upperBoundValue, oneMinusAlpha);
	}

	@Override
	protected Matrix getMeanFromDistribution() {
		Matrix mean = getDistribution().getMean();
		Matrix variance = getDistribution().getVariance();
		return mean.add(variance.scalarMultiply(0.5)).expMatrix();
	}
	
	@Override
	protected SymmetricMatrix getVarianceFromDistribution() {
		SymmetricMatrix variance = getDistribution().getVariance();
		return SymmetricMatrix.convertToSymmetricIfPossible(variance.expMatrix().scalarAdd(-1d).elementWiseMultiply(getMeanFromDistribution().elementWisePower(2)));
	}

	@Override
	public Matrix getRandomDeviate() {
		return super.getRandomDeviate().expMatrix();
	}



}
