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
package repicea.stats.distributions;

import repicea.math.Matrix;
import repicea.math.SymmetricMatrix;
import repicea.stats.CentralMomentsSettable;

/**
 * This class implements the Gaussian probability density function.
 * @author Mathieu Fortin - August 2012
 */
@SuppressWarnings("serial")
public class GaussianDistribution extends StandardGaussianDistribution implements CentralMomentsSettable {
		
	/**
	 * Constructor. <br>
	 * <br>
	 * Creates a Gaussian distribution with mean mu and variance sigma2. NOTE: Matrix sigma2 must be 
	 * positive definite.
	 * 
	 * @param mu the mean of the function
	 * @param sigma2 the variance of the function
	 * @throws UnsupportedOperationException if the matrix sigma2 is not positive definite
	 */
	public GaussianDistribution(Matrix mu, SymmetricMatrix sigma2) {
		setMean(mu);
		setVariance(sigma2);
	}
	
	/**
	 * Constructor for univariate Gaussian distribution.
	 * @param mean
	 * @param variance
	 */
	public GaussianDistribution(double mean, double variance) {
		Matrix mu = new Matrix(1,1);
		mu.setValueAt(0, 0, mean);
		setMean(mu);
		SymmetricMatrix sigma2 = new SymmetricMatrix(1);
		sigma2.setValueAt(0, 0, variance);
		setVariance(sigma2);
//		this(new Matrix(1,1,mean,0d), new Matrix(1,1,variance,0d));
	}
	
	/**
	 * Constructor for univariate Gaussian distribution centered on 0 with variance 1.
	 */
	public GaussianDistribution() {
		this(0d, 1d);
	}

	@Override
	public void setMean(Matrix mean) {
		super.setMean(mean);
	}

	@Override
	public void setVariance(SymmetricMatrix variance) {
		super.setVariance(variance);
	}

}
