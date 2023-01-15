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
import repicea.stats.Distribution;

/**
 * The UnknownDistribution class is a basic distribution for which the first two moments are known (and can be set) but the
 * underlying distribution remains unknown. NOTE the method getRandomObservation() always returns null for this class.
 * @author Mathieu Fortin - August 2012
 */
@SuppressWarnings("serial")
public class UnknownDistribution implements Distribution, CentralMomentsSettable {

	private Matrix mean;
	private SymmetricMatrix variance;
	
	
	@Override
	public boolean isMultivariate() {
		return ((Matrix) mean).m_iRows > 1;
	}

	@Override
	public Matrix getMean() {
		return mean;
	}

	@Override
	public SymmetricMatrix getVariance() {
		return variance;
	}

	@Override
	public Type getType() {
		return Type.UNKNOWN;
	}

	@Override
	public void setMean(Matrix mean) {
		this.mean = mean;
	}

	@Override
	public void setVariance(SymmetricMatrix variance) {
		this.variance = variance;
	}

	@Override
	public boolean isParametric() {
		return false;
	}

	@Override
	public Matrix getRandomRealization() {
		return null;
	}

//	@Override
//	public double getQuantile(double... values) {
//		return -1;
//	}


}
