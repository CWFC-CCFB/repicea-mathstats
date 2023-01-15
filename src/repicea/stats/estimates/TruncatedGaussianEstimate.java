/*
 * This file is part of the repicea library.
 *
 * Copyright (C) 2009-2015 Mathieu Fortin for Rouge-Epicea
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
import repicea.stats.distributions.BoundedDistribution;
import repicea.stats.distributions.TruncatedGaussianDistribution;

/**
 * The TruncatedGaussianEstimate class allows to generate random deviates from a 
 * truncated Gaussian distribution. The bound of the distribution can set through 
 * the setLowerBound and setUpperBound methods. In its current form, the class is only
 * available for univariate distributions.
 * @author Mathieu Fortin - August 2015
 */
public class TruncatedGaussianEstimate extends Estimate<TruncatedGaussianDistribution> implements BoundedDistribution {

	private static final long serialVersionUID = -1426033536030992926L;

	/**
	 * Basic constructor with mu set to 0 and sigma2 set to 1.
	 */
	public TruncatedGaussianEstimate() {
		super(new TruncatedGaussianDistribution());
	}

	/**
	 * Constructor 2 with user specified mu and sigma2.
	 * @param mu a double
	 * @param sigma2 a double
	 */
	public TruncatedGaussianEstimate(double mu, double sigma2) {
		super(new TruncatedGaussianDistribution(mu, sigma2));
	}
	
	@Override
	public void setLowerBoundValue(Matrix lowerBoundValue) {
		getDistribution().setLowerBoundValue(lowerBoundValue);
	}

	@Override
	public void setUpperBoundValue(Matrix upperBoundValue) {
		getDistribution().setUpperBoundValue(upperBoundValue);
	}

	@Override
	public ConfidenceInterval getConfidenceIntervalBounds(double oneMinusAlpha) {
		return null;
	}
	

	
	
}
