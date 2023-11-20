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
package repicea.stats.estimates;

import repicea.math.ComplexMatrix;
import repicea.math.HermitianMatrix;
import repicea.math.Matrix;
import repicea.stats.distributions.ComplexEmpiricalDistribution;

/**
 * An implementation of the Monte Carlo estimator for complex numbers.
 * @author Mathieu Fortin - November 2023
 */
@SuppressWarnings("serial")
public class ComplexMonteCarloEstimate extends ResamplingBasedEstimate<ComplexMatrix, HermitianMatrix>{

	/**
	 * Constructor.
	 */
	public ComplexMonteCarloEstimate() {
		super(new ComplexEmpiricalDistribution());
	}

	@Override
	public ComplexEmpiricalDistribution getDistribution() {
		return (ComplexEmpiricalDistribution) super.getDistribution();
	}
	
	@Override
	public ConfidenceInterval getConfidenceIntervalBounds(double oneMinusAlpha) {
		// TODO Auto-generated method stub
		return null;
	}
	
	/**
	 * Produce the real part of variance of diagonal element of the Hermitian matrix.
	 * @see ComplexEmpiricalDistribution#getVarianceRealPart()
	 * @return a Matrix instance
	 */
	public Matrix getVarianceRealPart() {
		return getDistribution().getVarianceRealPart();
	}

	/**
	 * Produce the imaginary part of variance of diagonal element of the Hermitian matrix.
	 * @see ComplexEmpiricalDistribution#getVarianceImaginaryPart()
	 * @return a Matrix instance
	 */
	public Matrix getVarianceImaginaryPart() {
		return getDistribution().getVarianceImaginaryPart();
	}

}
