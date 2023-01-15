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
package repicea.stats.estimates;

import java.security.InvalidParameterException;

import repicea.math.Matrix;
import repicea.math.SymmetricMatrix;
import repicea.stats.StatisticalUtility;

/**
 * The JackknifeEstimate class implements the variance calculation typical of
 * delete-d Jackknife variance estimator.
 * @author Mathieu Fortin - February 2019
 */
@SuppressWarnings("serial")
public class JackknifeEstimate extends ResamplingBasedEstimate {

	private final int n;
	private final int d;		// delete d observations here

	/**
	 * Constructor.
	 * @param n the number of observations in the original sample
	 * @param d the number of left-out observations
	 */
	public JackknifeEstimate(int n, int d) {
		super();
		this.d = d;
		this.n = n;
	}

	@Override
	protected SymmetricMatrix getVarianceFromDistribution() {
		long nCombinations = StatisticalUtility.getCombinations(n, d);
		long nRealizations = getNumberOfRealizations();
		if (nCombinations != nRealizations) {
			throw new InvalidParameterException("The number of realizations is inconsistent with the n and d parameters of the constructor!");
		}
		SymmetricMatrix ss = getDistribution().getVariance().scalarMultiply(getNumberOfRealizations() - 1);	// sum of squared difference
		double scalingFactor = ((double) n - d) / (d * nCombinations);
		return ss.scalarMultiply(scalingFactor);
	}
	

	@Override
	public ConfidenceInterval getConfidenceIntervalBounds(double oneMinusAlpha) {
		// TODO To be implemented
		return null;
	}

}
