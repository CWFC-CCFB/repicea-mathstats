/*
 * This file is part of the repicea library.
 *
 * Copyright (C) 2009-2019 Mathieu Fortin for Rouge-Epicea
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

import java.security.InvalidParameterException;

import repicea.math.Matrix;
import repicea.math.SymmetricMatrix;
import repicea.stats.Distribution;
import repicea.stats.StatisticalUtility;

/**
 * This class implements the standard Student's t distribution. The algorithm behind 
 * the random deviate generation is that of Bailey (1994) based on polar generation.
 * @author Mathieu Fortin - October 2019
 */
@SuppressWarnings("serial")
public class StudentTDistribution implements Distribution {

	
	private Matrix mu;
	private SymmetricMatrix sigma2;
	private Matrix lowerCholTriangle;
	private final int degreesOfFreedom;
	
	/**
	 * This constructor creates a Student't distribution with mean mu and variance sigma2. NOTE: It 
	 * does not handle multivariate distribution. If the mean and the variance are of dimensions > 1,
	 * they are assumed to be independent, i.e. the non zero off diagonal elements in the variance 
	 * matrix will not be considered.
	 * @param mu the mean of the function
	 * @param variance the variance of the distribution
	 * @param degreesOfFreedom the degrees of freedom of the distribution
	 */
	public StudentTDistribution(Matrix mu, SymmetricMatrix variance, int degreesOfFreedom) {
		this.degreesOfFreedom = degreesOfFreedom;
		setMean(mu);
		setVariance(variance);
	}
	
	/**
	 * Constructor for a standard Student's distribution, i.e. with mean 0 and sigma^2 = 1.
	 * @param degreesOfFreedom
	 */
	public StudentTDistribution(int degreesOfFreedom) {
		this.degreesOfFreedom = degreesOfFreedom;
		Matrix mu = new Matrix(1,1);
		SymmetricMatrix variance = new SymmetricMatrix(1);
		variance.setValueAt(0, 0, 1d);
		setMean(mu);
		setVariance(variance);
	}
	
	@Override
	public boolean isMultivariate() {
		return getMu().m_iRows > 1;
	}
	
	@Override
	public Matrix getRandomRealization() {
		int nbRows = getMean().m_iRows;
		Matrix realization = new Matrix(nbRows,1);
		for (int i = 0; i < nbRows; i++) {
			double deviate = StatisticalUtility.getRandom().nextStudentT(degreesOfFreedom);
			double scaledDeviate = getStandardDeviation().getValueAt(i, i) * deviate + getMean().getValueAt(i, 0); 
			realization.setValueAt(i, 0, scaledDeviate);
		}
		return realization;
	}

	/**
	 * This method returns the lower triangle of the Cholesky decomposition of the variance-covariance matrix.
	 * @return a Matrix instance
	 */
	public Matrix getStandardDeviation() {
		if (lowerCholTriangle == null) {
			lowerCholTriangle = getSigma2().getLowerCholTriangle();
		}
		return lowerCholTriangle;
	}
	
	@Override
	public Matrix getMean() {return getMu();}

	@Override
	public SymmetricMatrix getVariance() {return getSigma2().scalarMultiply(((double) degreesOfFreedom)/(degreesOfFreedom - 2));}

	@Override
	public Type getType() {return Distribution.Type.STUDENT;}

	private void setMean(Matrix mu) {
		this.mu = mu;
	}

	private void setVariance(SymmetricMatrix sigma2) {
		if (sigma2 != null && !sigma2.isSymmetric()) {
			throw new InvalidParameterException("The variance-covariance matrix must be symmetric!");
		}
		this.sigma2 = sigma2;
		lowerCholTriangle = null;
	}

	private Matrix getMu() {return mu;}
	
	private SymmetricMatrix getSigma2() {return sigma2;}
	
	
	@Override
	public boolean isParametric() {return true;}

}
