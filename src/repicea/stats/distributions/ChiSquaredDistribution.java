/*
 * This file is part of the repicea library.
 *
 * Copyright (C) 2009-2017 Mathieu Fortin for Rouge-Epicea
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
 * The ChiSquareDistribution class represents a univariate Chi Square distribution with a given degrees of freedom.
 * @author Mathieu Fortin - November 2012
 */
public final class ChiSquaredDistribution implements Distribution {

	private static final long serialVersionUID = 20121114L;

	private final int degreesOfFreedom;

	private SymmetricMatrix mean;
	private SymmetricMatrix variance;
	
	private Matrix lowerTriangleChol;
	private Matrix upperTriangleChol;
	

	/**
	 * Constructor for univariate Chi-squared distribution. 
	 * @param degreesOfFreedom the degrees of freedom
	 * @param meanValue the mean value
	 */
	public ChiSquaredDistribution(int degreesOfFreedom, double meanValue) {
		if (degreesOfFreedom < 1) {
			throw new InvalidParameterException("The number of degrees of freedom must be equal to or larger than 1!");
		}
		this.degreesOfFreedom = degreesOfFreedom;
		if (meanValue < 0) {
			throw new InvalidParameterException("The variance estimate must be larger than 0!");
		}
		SymmetricMatrix mean = new SymmetricMatrix(1);
		mean.setValueAt(0, 0, meanValue);
		this.mean = mean;
	}

	
	/**
	 * Constructor for multivariate Chi-squared distribution (Wishart distribution). 
	 * @param degreesOfFreedom the degrees of freedom
	 * @param meanValues the mean value
	 */
	public ChiSquaredDistribution(int degreesOfFreedom, SymmetricMatrix meanValues) {
		if (degreesOfFreedom < 1) {
			throw new InvalidParameterException("The number of degrees of freedom must be equal to or larger than 1!");
		}
		this.degreesOfFreedom = degreesOfFreedom;
		if (!meanValues.isSymmetric()) {
			throw new InvalidParameterException("The variance-covariance matrix must be symmetric!");
		}
		if (meanValues.diagonalVector().anyElementSmallerOrEqualTo(0d)) {
			throw new InvalidParameterException("The variance estimate must be larger than 0!");
		}
		this.mean = meanValues.getDeepClone();
	}


	/**
	 * This method returns the degrees of freedom.
	 * @return an integer
	 */
	public int getDegreesOfFreedom() {
		return degreesOfFreedom;
	}
	
	@Override
	public Matrix getMean() {
		return mean.getDeepClone();
	}

	@Override
	public SymmetricMatrix getVariance() {
		if (variance == null) {
			variance = calculateVarianceMatrix();
		}
		return variance;
	}

	/**
	 * This method is based on Wishart's distribution.
	 */
	private SymmetricMatrix calculateVarianceMatrix() {
		int nRows = mean.m_iRows;
		SymmetricMatrix varianceMat = new SymmetricMatrix(nRows);
		double result;
		for (int i = 0; i < nRows; i++) {
			for (int j = i; j < nRows; j++) {
				result = mean.getValueAt(i, j) * mean.getValueAt(i, j) + 
						mean.getValueAt(i, i) * mean.getValueAt(j, j);  
				varianceMat.setValueAt(i, j, result);
//				if (j != i) {
//					varianceMat.setValueAt(j, i, result);
//				}
			}
		}
		return varianceMat.scalarMultiply(1d / getDegreesOfFreedom());
	}
	
	

	@Override
	public boolean isParametric() {return true;}

	@Override
	public boolean isMultivariate() {
		return mean.m_iRows > 1;
	}

	@Override
	public Type getType() {
		if (isMultivariate()) {
			return Type.WISHART;
		} else {
			return Type.CHI_SQUARE;
		}
	}

	@Override
	public SymmetricMatrix getRandomRealization() {
		if (isMultivariate()) {	// then we use Bartlett decomposition for a Wishart distribution
			if (lowerTriangleChol == null) {
				lowerTriangleChol = mean.getLowerCholTriangle();
				upperTriangleChol = lowerTriangleChol.transpose();
			}
			Matrix aMat = StatisticalUtility.getRandom().nextBartlettDecompositionMatrix(getDegreesOfFreedom(), mean.m_iRows);
			Matrix randomMat = lowerTriangleChol.multiply(aMat).multiply(aMat.transpose()).multiply(upperTriangleChol).scalarMultiply(1d/getDegreesOfFreedom());
			return SymmetricMatrix.convertToSymmetricIfPossible(randomMat);
		} else {
			double factor = StatisticalUtility.getRandom().nextChiSquare(getDegreesOfFreedom()) / getDegreesOfFreedom();
			return mean.scalarMultiply(factor);
		}
	}


}
