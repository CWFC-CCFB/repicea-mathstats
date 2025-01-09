/*
 * This file is part of the repicea-mathstats library.
 *
 * Copyright (C) 2025 His Majesty the King in right of Canada
 * Author: Mathieu Fortin, Canadian Forest Service
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
import java.util.List;
import java.util.Map;
import java.util.concurrent.ConcurrentHashMap;

import repicea.math.Matrix;
import repicea.math.SymmetricMatrix;
import repicea.math.utility.GaussianUtility;
import repicea.stats.distributions.GaussianDistribution;
import repicea.stats.sampling.PopulationUnit;

/**
 * A class to implement stratified sampling design.<p>
 * The point estimate is that of the total.
 * @author Mathieu Fortin - January 2025
 */
@SuppressWarnings("serial")
public class StratifiedPopulationTotalEstimate extends AbstractEstimate<Matrix, SymmetricMatrix, GaussianDistribution> 
												implements PointEstimate {

	private final Map<String, PopulationTotalEstimate> stratumDesign;
	
	/**
	 * Constructor.
	 * @param stratumNames a List of strings that stand for the stratum names.
	 * @param strataPopulationSizes a List of Double that stand for the stratum sample sizes.
	 */
	public StratifiedPopulationTotalEstimate(List<String> stratumNames, List<Double> strataPopulationSizes) {
		super(new GaussianDistribution(0d, 1d));
		if (stratumNames == null || stratumNames.isEmpty()) {
			throw new InvalidParameterException("The strataNames argument must be non null and not empty!");
		}
		if (strataPopulationSizes == null || strataPopulationSizes.isEmpty()) {
			throw new InvalidParameterException("The strataPopulationSizes argument must be non null and not empty!");
		}
		if (stratumNames.size() != strataPopulationSizes.size()) {
			throw new InvalidParameterException("The strataNames and strataPopulationSizes arguments must have the same size!");
		}
		stratumDesign = new ConcurrentHashMap<String, PopulationTotalEstimate>();
		for (int i = 0; i < stratumNames.size(); i++) {
			String stratumName = stratumNames.get(i);
			double popSize = strataPopulationSizes.get(i);
			if (!stratumDesign.containsKey(stratumName)) {
				stratumDesign.put(stratumName, new PopulationTotalEstimate(popSize));
			} else {
				throw new InvalidParameterException("This stratum name appears twice in the strataNames argument: " + stratumName);
			}
		}
	}
	
	@Override
	public double getPopulationSize() {
		double popSize = 0d;
		for (PopulationTotalEstimate subDomain : stratumDesign.values()) {
			popSize += subDomain.getPopulationSize();
		}
		return popSize;
	}

	public void addObservation(String stratumName, PopulationUnit pu) {
		if (!stratumDesign.containsKey(stratumName)) {
			throw new InvalidParameterException("This stratum name has not been specified in the constructor: " + stratumName);
		}
		stratumDesign.get(stratumName).addObservation(pu);
	}

	@Override
	protected Matrix getMeanFromDistribution() {
		Matrix total = null;
		for (PopulationTotalEstimate estimate : stratumDesign.values()) {
			total = total == null ? 
					estimate.getMean() : 
						total.add(estimate.getMean());
		}
		return total;
	}

	@Override
	protected SymmetricMatrix getVarianceFromDistribution() {
		SymmetricMatrix variance = null;
		for (PopulationTotalEstimate estimate : stratumDesign.values()) {
			variance = variance == null ? 
					estimate.getVariance() : 
						(SymmetricMatrix) variance.add(estimate.getVariance());
		}
		return variance;
	}

	protected final Matrix getQuantileForProbability(double probability) {
		Matrix stdDev = getVariance().diagonalVector().elementWisePower(.5); 
		double quantile = GaussianUtility.getQuantile(probability);
		return getMean().add(stdDev.scalarMultiply(quantile));
	}
	
	@Override
	public ConfidenceInterval getConfidenceIntervalBounds(double oneMinusAlpha) {
		Matrix lowerBoundValue = getQuantileForProbability(.5 * (1d - oneMinusAlpha));
		Matrix upperBoundValue = getQuantileForProbability(1d - .5 * (1d - oneMinusAlpha));
		return new ConfidenceInterval(lowerBoundValue, upperBoundValue, oneMinusAlpha);
	}

}
