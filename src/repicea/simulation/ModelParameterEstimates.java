/*
 * This file is part of the repicea-simulation library.
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
package repicea.simulation;

import java.util.ArrayList;
import java.util.List;

import repicea.math.Matrix;
import repicea.math.SymmetricMatrix;
import repicea.stats.Distribution;
import repicea.stats.StatisticalUtility;
import repicea.stats.estimates.GaussianEstimate;

/**
 * The ModelParameterEstimates is a wrapper for parameter estimates. It includes a common variance-
 * covariance matrix for the parameter estimates. 
 *
 * @author Mathieu Fortin - September 2018
 * 
 */
@SuppressWarnings("serial")
public class ModelParameterEstimates extends GaussianEstimate {
	
	protected final List<Integer> estimatedParameterIndices;

	/**
	 * Constructor.
	 * @param mean a vector that corresponds to the mean value
	 * @param variance a symmetric positive definite matrix 
	 */
	public ModelParameterEstimates(Matrix mean, SymmetricMatrix variance) {
		super(mean, variance);
		estimatedParameterIndices = new ArrayList<Integer>();
		setEstimatedParameterIndices();
	}
	
	protected void setEstimatedParameterIndices() {
		for (int i = 0; i < getMean().m_iRows; i++) {
			estimatedParameterIndices.add(i);
		}
	}
	
	/**
	 * This method returns the indices of the parameters that were truly estimated. In the case of 
	 * a SAS implementation, some parameters may actually be fake.. 
	 * @return a List of Integer which is a copy of the original list to avoid modifications.
	 */
	public List<Integer> getTrueParameterIndices() {
		List<Integer> copyList = new ArrayList<Integer>();
		copyList.addAll(estimatedParameterIndices);
		return copyList;
	}

	@Override
	public Matrix getRandomDeviate() {
		Matrix lowerChol = getDistribution().getStandardDeviation();
		Matrix randomVector = StatisticalUtility.drawRandomVector(lowerChol.m_iRows, Distribution.Type.GAUSSIAN);
		Matrix oMat = lowerChol.multiply(randomVector);
		Matrix deviate = getMean().getDeepClone();
		deviate.addElementsAt(estimatedParameterIndices, oMat);
		return deviate;
	}

}