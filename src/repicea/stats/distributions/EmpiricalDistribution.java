/*
 * This file is part of the repicea library.
 *
 * Copyright (C) 2009-2016 Mathieu Fortin for Rouge-Epicea
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

/**
 * The EmpiricalDistribution class handles empirical distribution from observations
 * or realizations.
 * @author Mathieu Fortin - May 2018
 *
 */
public class EmpiricalDistribution extends AbstractEmpiricalDistribution<Matrix, SymmetricMatrix> {

	private static final long serialVersionUID = 20120826L;
	
	
	@Override
	public Matrix getMean() {
		if (observations == null || observations.isEmpty()) {
			return null;
		} else {
			Matrix sum = null;
			for (Matrix mat : observations) {
				if (sum == null) {
					sum = mat.getDeepClone();
				} else {
					sum = sum.add(mat);
				}
			}
			return sum.scalarMultiply(1d / observations.size());
		}
	}

	@Override
	public SymmetricMatrix getVariance() {
		Matrix mean = getMean();
		Matrix sse = null;
		Matrix error;
		for (Matrix mat : observations) {
			error = mat.subtract(mean);
			if (sse == null) {
				sse = error.multiply(error.transpose());
			} else {
				sse = sse.add(error.multiply(error.transpose()));
			}
		}
		SymmetricMatrix convertedSse = SymmetricMatrix.convertToSymmetricIfPossible(sse);
		return convertedSse.scalarMultiply(1d / (observations.size()-1));
	}


}
