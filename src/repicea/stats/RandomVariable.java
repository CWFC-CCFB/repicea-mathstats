/*
 * This file is part of the repicea library.
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
package repicea.stats;

import java.io.Serializable;

import repicea.math.AbstractMatrix;
import repicea.stats.estimates.DistributionProvider;

/**
 * The RandomVariable class ensures the instance can provide its two first central moments (i.e. mean and variance) as
 * well as a distribution. The parameter D sets the distribution of the random variables (see the repicea.stats.distributions 
 * package)
 * @author Mathieu Fortin - May 2012
 * @param <D> a Distribution derived instance
 */
public abstract class RandomVariable<M extends AbstractMatrix, V extends AbstractMatrix, D extends Distribution<M,V>> implements MomentGettable<M,V>, 
																	Serializable,
																	DistributionProvider<D> {

	private static final long serialVersionUID = 1L;
	
	private final D distribution;
	
	protected RandomVariable(D distribution) {
		this.distribution = distribution;
	}
	
	@Override
	public D getDistribution() {
		return distribution;
	}
	
	@Override
	public M getMean() {
		return getMeanFromDistribution();
	}

	protected M getMeanFromDistribution() {
		return getDistribution().getMean();
	}
	
	@Override
	public V getVariance() {
		return getVarianceFromDistribution();
	}

	protected V getVarianceFromDistribution() {
		return getDistribution().getVariance();
	}
}
