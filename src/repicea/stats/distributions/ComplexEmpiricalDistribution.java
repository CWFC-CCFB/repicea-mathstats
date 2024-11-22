/*
 * This file is part of the repicea-mathstats library.
 *
 * Copyright (C) 2023 His Majesty the King in Right of Canada
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
package repicea.stats.distributions;

import java.util.ArrayList;
import java.util.List;
import java.util.stream.Collectors;

import repicea.math.ComplexMatrix;
import repicea.math.ComplexNumber;
import repicea.math.ComplexSymmetricMatrix;
import repicea.math.HermitianMatrix;
import repicea.math.Matrix;

/**
 * The ComplexEmpiricalDistribution class supports Monte Carlo estimation based
 * on complex numbers.
 * @author Mathieu Fortin - November 2023
 */
@SuppressWarnings("serial")
public class ComplexEmpiricalDistribution extends AbstractEmpiricalDistribution<ComplexMatrix, HermitianMatrix> {

	@Override
	public ComplexMatrix getMean() {
		ComplexMatrix mean = null;
		for (ComplexMatrix m : getRealizations()) {
			mean = mean == null ? m : mean.add(m);
		}
		return mean.scalarMultiply(1d / getNumberOfRealizations());
	}

	/**
	 * For complex random variable, the variance-covariance matrix is 
	 * a Hermitian matrix.
	 * 
	 * @see AbstractEmpiricalDistribution#getVariance() 
	 * @return a HermitianMatrix instance
	 */
	@Override
	public HermitianMatrix getVariance() {
		ComplexMatrix mean = getMean();
		ComplexMatrix variance = null;
		for (ComplexMatrix m : getRealizations()) {
			ComplexMatrix res = m.subtract(mean);
			ComplexMatrix product = res.multiply(res.getComplexConjugate().transpose());
			variance = variance == null ? product : variance.add(product);
		}
		variance = variance.scalarMultiply(1d/(getNumberOfRealizations() - 1));
		HermitianMatrix formattedMatrix = new HermitianMatrix(variance);
		return formattedMatrix;
	}

	/**
	 * Compute a vector of values that are the real part of the diagonal 
	 * of the Hermitian matrix.
	 * @see ComplexEmpiricalDistribution#getVariance()	
	 * @return a Matrix instance
	 */
	public Matrix getVarianceRealPart() {
		ComplexMatrix mean = getMean();
		double[] sumRes = new double[mean.m_iRows];
		for (ComplexMatrix m : getRealizations()) {
			for (int i = 0; i < mean.m_iRows; i++) {
				ComplexNumber n = m.getValueAt(i, 0);
				double diff = n.realPart - mean.getValueAt(i, 0).realPart;
				sumRes[i] += diff * diff;
			}
		}
		for (int i = 0; i < sumRes.length; i++) {
			sumRes[i] /= getNumberOfRealizations() - 1;
		}
		return new Matrix(sumRes);
	}
	
	/**
	 * Compute a vector of values that are the imaginary part of the diagonal 
	 * of the Hermitian matrix.
	 * @see ComplexEmpiricalDistribution#getVariance()	
	 * @return a Matrix instance
	 */
	public Matrix getVarianceImaginaryPart() {
		ComplexMatrix mean = getMean();
		double[] sumRes = new double[mean.m_iRows];
		for (ComplexMatrix m : getRealizations()) {
			for (int i = 0; i < mean.m_iRows; i++) {
				ComplexNumber n = m.getValueAt(i, 0);
				double diff = n.imaginaryPart - mean.getValueAt(i, 0).imaginaryPart;
				sumRes[i] += diff * diff;
			}
		}
		for (int i = 0; i < sumRes.length; i++) {
			sumRes[i] /= getNumberOfRealizations() - 1;
		}
		return new Matrix(sumRes);
	}

	/**
	 * Provide the pseudovariance of the complex
	 * random matrices.
	 * 
	 * @return a ComplexSymmetricMatrix instance
	 */
	public ComplexSymmetricMatrix getPseudoVariance() {
		ComplexMatrix variance = null;
		for (ComplexMatrix residual : getResiduals()) {
			ComplexMatrix product = residual.multiply(residual.transpose());
			variance = variance == null ? product : variance.add(product);
		}
		variance = variance.scalarMultiply(1d/(getNumberOfRealizations() - 1));
		ComplexSymmetricMatrix formattedMatrix = new ComplexSymmetricMatrix(variance);
		return formattedMatrix;
	}

	/**
	 * Provide imaginary realizations.<p>
	 * The imaginary realizations are calculated as the 
	 * residuals times i plus the mean.
	 * @return a List of ComplexMatrix instances
	 */
	public List<ComplexMatrix> getImaginaryRealizations() {
		ComplexMatrix mean = getMean();
		List<ComplexMatrix> realizations = getResiduals().
				stream().
				map(e -> e.scalarMultiply(ComplexNumber.COMPLEX_I).add(mean)).
				collect(Collectors.toList());
		return realizations;
	}
	
	/**
	 * Provide the residuals of the realizations.<p>
	 * Those are the realizations minus the mean.
	 * @return a List of ComplexMatrix instances
	 */
	private List<ComplexMatrix> getResiduals() {
		List<ComplexMatrix> residuals = new ArrayList<ComplexMatrix>();
		ComplexMatrix mean = getMean();
		for (ComplexMatrix m : getRealizations()) {
			residuals.add(m.subtract(mean));
		}
		return residuals;
	}
	
}
