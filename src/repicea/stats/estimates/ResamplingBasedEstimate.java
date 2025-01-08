package repicea.stats.estimates;

import java.security.InvalidParameterException;
import java.util.List;

import repicea.math.AbstractMatrix;
import repicea.math.Matrix;
import repicea.stats.distributions.AbstractEmpiricalDistribution;

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
@SuppressWarnings({ "serial", "rawtypes" })
abstract class ResamplingBasedEstimate<M extends AbstractMatrix, V extends AbstractMatrix> 
					extends AbstractEstimate<M, V, AbstractEmpiricalDistribution<M,V>> implements NumberOfRealizationsProvider {
	
	/**
	 * Constructor.
	 */
	protected ResamplingBasedEstimate(AbstractEmpiricalDistribution<M, V> dist) {
		super(dist);
		estimatorType = EstimatorType.Resampling;
	}


	/**
	 * This method adds a realization to the empirical distribution. The conformity of the
	 * new realization is checked before adding it.
	 * @param value a Matrix
	 */
	public void addRealization(M value) {
		if (checkConformity(value)) {
			getDistribution().addRealization(value);
		} else {
			throw new InvalidParameterException("The matrix is not conform to previous observations!");
		}
	}

	private boolean checkConformity(M value) {
		List<M> observations = getDistribution().getRealizations();
		if (observations.isEmpty()) {
			return true; 
		} else {
			M firstObservation = observations.get(0);
			return (firstObservation.m_iRows == value.m_iRows && firstObservation.m_iCols == value.m_iCols);
		}
	}


	/**
	 * Provide the quantile associated to a particular probability. 
	 * @param probability the probability level
	 * @return a Matrix instance that contains the quantiles
	 */
	protected abstract Matrix getQuantileForProbability(double probability); 

	/** 
	 * {@inheritDoc}<p>
	 * 
	 * For classes involving complex numbers, the confidence intervals are based on the 
	 * real part of the realizations.
	 */
	@Override
	public ConfidenceInterval getConfidenceIntervalBounds(double oneMinusAlpha) {
		Matrix lowerBoundValue = getQuantileForProbability(.5 * (1d - oneMinusAlpha));
		Matrix upperBoundValue = getQuantileForProbability(1d - .5 * (1d - oneMinusAlpha));
		return new ConfidenceInterval(lowerBoundValue, upperBoundValue, oneMinusAlpha);
	}
	
	/**
	 * This method returns the list of realizations in the empirical distribution.
	 * @return a List of Matrix instance
	 */
	public List<M> getRealizations() {return getDistribution().getRealizations();}

	@Override
	public int getNumberOfRealizations() {return getDistribution().getNumberOfRealizations();}

}
