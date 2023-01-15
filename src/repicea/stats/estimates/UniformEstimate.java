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
import repicea.stats.distributions.UniformDistribution;

/**
 * The UniformEstimate class relies on a uniform distribution.
 * @author Mathieu Fortin - April 2016
 *
 */
@SuppressWarnings("serial")
public class UniformEstimate extends Estimate<UniformDistribution> implements BoundedDistribution {

	public UniformEstimate(Matrix lowerBoundValue, Matrix upperBoundValue) {
		super(new UniformDistribution(lowerBoundValue, upperBoundValue));
	}

	@Override
	public void setLowerBoundValue(Matrix lowerBoundValue) {
		getDistribution().setLowerBoundValue(lowerBoundValue);
	}

	@Override
	public void setUpperBoundValue(Matrix upperBoundValue) {
		getDistribution().setUpperBoundValue(upperBoundValue);
	}

	protected Matrix getQuantileForProbability(double probability) {
		Matrix lowerBound = getDistribution().getLowerBound().getBoundValue();
		Matrix upperBound = getDistribution().getUpperBound().getBoundValue();
		Matrix boundDistance = upperBound.subtract(lowerBound);
		return lowerBound.add(boundDistance.scalarMultiply(probability));
	}

	@Override
	public ConfidenceInterval getConfidenceIntervalBounds(double oneMinusAlpha) {
		Matrix lowerBoundValue = getQuantileForProbability(.5 * (1d - oneMinusAlpha));
		Matrix upperBoundValue = getQuantileForProbability(1d - .5 * (1d - oneMinusAlpha));
		return new ConfidenceInterval(lowerBoundValue, upperBoundValue, oneMinusAlpha);
	}

}
