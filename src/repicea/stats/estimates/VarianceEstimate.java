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

import repicea.math.Matrix;
import repicea.stats.distributions.ChiSquaredDistribution;

/**
 * The VarianceEstimate is a class of Estimate that suits the variance estimation. The distribution
 * of the estimate is assumed to be a Chi Square.
 * @author Mathieu Fortin - November 2012
 */
public class VarianceEstimate extends Estimate<Matrix, ChiSquaredDistribution> {

	private static final long serialVersionUID = 20121114L;

	/**
	 * Constructor with mean only.
	 * @param degreesOfFreedom the degrees of freedom
	 * @param mean the estimated mean
	 */
	public VarianceEstimate(int degreesOfFreedom, double mean) {
		super(new ChiSquaredDistribution(degreesOfFreedom, mean));
		estimatorType = EstimatorType.LikelihoodBased;
	}
	
	/**
	 * Provide the degrees of freedom associated with the variance estimate.
	 * @return an integer
	 */
	public int getDegreesOfFreedom() {
		return getDistribution().getDegreesOfFreedom();
	}
	

	@Override
	public ConfidenceInterval getConfidenceIntervalBounds(double oneMinusAlpha) {
		return null;
	}

	
}
