package repicea.stats.estimates;

import java.security.InvalidParameterException;
import java.util.List;

import repicea.math.Matrix;
import repicea.stats.distributions.EmpiricalDistribution;

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
@SuppressWarnings("serial")
abstract class ResamplingBasedEstimate extends Estimate<EmpiricalDistribution> implements NumberOfRealizationsProvider {
	
	/**
	 * Constructor.
	 */
	protected ResamplingBasedEstimate() {
		super(new EmpiricalDistribution());
		estimatorType = EstimatorType.Resampling;
	}


	/**
	 * This method adds a realization to the empirical distribution. The conformity of the
	 * new realization is checked before adding it.
	 * @param value a Matrix
	 */
	public void addRealization(Matrix value) {
		if (checkConformity(value)) {
			getDistribution().addRealization(value);
		} else {
			throw new InvalidParameterException("The matrix is not conform to previous observations!");
		}
	}

	private boolean checkConformity(Matrix value) {
		List<Matrix> observations = getDistribution().getRealizations();
		if (observations.isEmpty()) {
			return true; 
		} else {
			Matrix firstObservation = observations.get(0);
			if (value == null || firstObservation == null) {
				int u = 0;
			}
			return (firstObservation.m_iRows == value.m_iRows && firstObservation.m_iCols == value.m_iCols);
		}
	}


	/**
	 * This method returns the list of realizations in the empirical distribution.
	 * @return a List of Matrix instance
	 */
	public List<Matrix> getRealizations() {return getDistribution().getRealizations();}

	@Override
	public int getNumberOfRealizations() {return getDistribution().getNumberOfRealizations();}

}
