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

import java.io.Serializable;

import repicea.math.Matrix;
import repicea.math.SymmetricMatrix;
import repicea.stats.CentralMomentsSettable;
import repicea.stats.distributions.UnknownDistribution;

/**
 * The Estimate class is the general class for all estimates.
 * @author Mathieu Fortin - July 2011
 */
public class SimpleEstimate extends Estimate<Matrix, UnknownDistribution> implements CentralMomentsSettable, Serializable {

	private static final long serialVersionUID = 20110912L;
	
	/**
	 * Public constructor 1 for derived classes.
	 */
	public SimpleEstimate() {
		super(new UnknownDistribution());
		estimatorType = EstimatorType.Unknown;
	}

	/**
	 * Public constructor 2 with mean and variance
	 * @param mean a Matrix instance
	 * @param variance a Matrix instance
	 */
	public SimpleEstimate(Matrix mean, SymmetricMatrix variance) {
		this();
		setMean(mean);
		setVariance(variance);
	}

	@Override
	public void setMean(Matrix mean) {
		getDistribution().setMean(mean);
	}
	
	@Override
	public void setVariance(SymmetricMatrix variance) {
		getDistribution().setVariance(variance);
	}


	@Override
	public ConfidenceInterval getConfidenceIntervalBounds(double oneMinusAlpha) {
		return null; // as no specific distribution is assumed
	}
	
	
}
