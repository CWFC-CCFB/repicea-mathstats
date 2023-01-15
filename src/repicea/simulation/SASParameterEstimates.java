/*
 * This file is part of the repicea-foresttools library.
 *
 * Copyright (C) 2009-2016 Mathieu Fortin for Rouge-Epicea
 *
 * This library is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public
 * License as published by the Free Software Foundation; either
 * version 2.1 of the License, or (at your option) any later version.
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

import java.security.InvalidParameterException;

import repicea.math.Matrix;
import repicea.math.SymmetricMatrix;

/**
 * The SASParameterEstimates class is customized for SAS outputs. The major difference
 * is related to how the random deviates are calculated. Since SAS produces false estimates,
 * such as 0 for reference class or 1 for offset, the getRandomDeviate() method makes sure
 * the false estimate are not accounted for during the simulation of the random deviates.
 * @author Mathieu Fortin - September 2012
 */
@SuppressWarnings("serial")
public class SASParameterEstimates extends ModelParameterEstimates {


	/**
	 * Constructor.
	 * @param mean a vector that corresponds to the mean value
	 * @param variance a symmetric positive definite matrix 
	 */
	public SASParameterEstimates(Matrix mean, SymmetricMatrix variance) {
		super(mean, variance);
	}

	@Override
	protected void setEstimatedParameterIndices() {
		Matrix mean = getMean();
		for (int i = 0; i < mean.m_iRows; i++) {
			if (mean.getValueAt(i, 0) != 0d && mean.getValueAt(i, 0) != 1d) { 
				estimatedParameterIndices.add(i);
			}
		}
		Matrix variance = getVariance();
		if (variance!= null  && variance.m_iRows != estimatedParameterIndices.size()) {
			throw new InvalidParameterException("SASParameterEstimates: the variance matrix is not compatible with the vector of parameter estimates");
		}
	}
	
	

//	@Override
//	public Matrix getRandomDeviate() {
//		Matrix lowerChol = getDistribution().getStandardDeviation();
//		Matrix randomVector = StatisticalUtility.drawRandomVector(lowerChol.m_iRows, Distribution.Type.GAUSSIAN);
//		Matrix oMat = lowerChol.multiply(randomVector);
//		Matrix deviate = getMean().getDeepClone();
//		deviate.addElementsAt(estimatedParameterIndices, oMat);
//		return deviate;
//	}

//	/**
//	 * This method returns a list of indices for variance extraction. This method is needed
//	 * because the SAS implementation includes some fake parameters.
//	 * @param parameterIndices a List of integer
//	 * @return a List of integer instances
//	 */
//	protected List<Integer> getVarianceIndicesForThoseParameterIndices(List<Integer> parameterIndices) {
//		List<Integer> varianceIndices = new ArrayList<Integer>();
//		for (Integer paramIndex : parameterIndices) {
//			if (!estimatedParameterIndices.contains(paramIndex)) {
//				throw new InvalidParameterException("The list contains some parameter indices that are not valid!");
//			} else {
//				varianceIndices.add(estimatedParameterIndices.indexOf(paramIndex));
//			}
//		}
//		return varianceIndices;
//	}
	

}
