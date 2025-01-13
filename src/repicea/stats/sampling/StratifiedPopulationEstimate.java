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
package repicea.stats.sampling;

import java.security.InvalidParameterException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.concurrent.ConcurrentHashMap;

import repicea.math.Matrix;
import repicea.math.SymmetricMatrix;
import repicea.stats.estimates.Estimate;

/**
 * A class to implement stratified sampling design.<p>
 * The point estimate is that of the total.
 * @author Mathieu Fortin - January 2025
 */
@SuppressWarnings("serial")
public class StratifiedPopulationEstimate extends AbstractPointEstimate 
												implements FinitePopulationPointEstimate {

	private final Map<String, FinitePopulationEstimate> stratumDesign;
	
	/**
	 * Constructor.
	 * @param stratumNames a List of strings that stand for the stratum names.
	 * @param strataPopulationSizes a List of Double that stand for the stratum sample sizes.
	 */
	public StratifiedPopulationEstimate(List<String> stratumNames, List<Double> strataPopulationSizes) {
		stratumDesign = new ConcurrentHashMap<String, FinitePopulationEstimate>();
		if (stratumNames == null || stratumNames.isEmpty()) {
			throw new InvalidParameterException("The strataNames argument must be non null and not empty!");
		}
		if (strataPopulationSizes == null || strataPopulationSizes.isEmpty()) {
			throw new InvalidParameterException("The strataPopulationSizes argument must be non null and not empty!");
		}
		if (stratumNames.size() != strataPopulationSizes.size()) {
			throw new InvalidParameterException("The strataNames and strataPopulationSizes arguments must have the same size!");
		}
		for (int i = 0; i < stratumNames.size(); i++) {
			String stratumName = stratumNames.get(i);
			double popSize = strataPopulationSizes.get(i);
			if (!stratumDesign.containsKey(stratumName)) {
				stratumDesign.put(stratumName, new FinitePopulationEstimate(popSize));
			} else {
				throw new InvalidParameterException("This stratum name appears twice in the strataNames argument: " + stratumName);
			}
		}
	}
	
	@Override
	public double getPopulationSize() {
		double popSize = 0d;
		for (FinitePopulationEstimate subDomain : stratumDesign.values()) {
			popSize += subDomain.getPopulationSize();
		}
		return popSize;
	}

	/**
	 * Validate the population unit before adding it to the
	 * observations. <p>
	 * 
	 * For this class, the method checks if the stratum name is found in 
	 * the stratumDesign member.
	 * 
	 * @param obs a Matrix instance instance
	 * @param obsId the observation id
	 * @param stratumName the stratum name
	 */
	@Override
	protected final void validateUnit(Matrix obs, String obsId, String stratumName) {
		if (!stratumDesign.containsKey(stratumName)) {
			throw new InvalidParameterException("This stratum has not been specified in the constructor: " + stratumName);
		}
	}
	
	@Override
	public final void addObservation(Matrix obs, String obsId, String stratumName) {
		super.addObservation(obs, obsId, stratumName);
		stratumDesign.get(stratumName).addObservation(obs, obsId);
	}

	
	@Override
	public final Matrix getTotal() {
		Matrix total = null;
		for (FinitePopulationEstimate estimate : stratumDesign.values()) {
			total = total == null ? 
					estimate.getTotal() : 
						total.add(estimate.getTotal());
		}
		return total;
	}

	@Override
	public final SymmetricMatrix getVarianceOfTotal() {
		SymmetricMatrix totalVariance = null;
		for (FinitePopulationEstimate estimate : stratumDesign.values()) {
			totalVariance = totalVariance == null ? 
					estimate.getVarianceOfTotal() : 
						(SymmetricMatrix) totalVariance.add(estimate.getVarianceOfTotal());
		}
		return totalVariance;
	}

	
	
	@Override
	protected final Matrix getMeanFromDistribution() {
		return getTotal().scalarMultiply(1d / getPopulationSize());
	}

	@Override
	protected final SymmetricMatrix getVarianceFromDistribution() {
		double squaredPopulationSize = getPopulationSize() * getPopulationSize();
		return getVarianceOfTotal().scalarMultiply(1d / squaredPopulationSize);
	}

	@Override
	protected boolean isMergeableEstimate(Estimate<?, ?, ?> estimate) {
		if (getClass().equals(estimate.getClass())) {
			StratifiedPopulationEstimate pe = (StratifiedPopulationEstimate) estimate;
			if (nCols == pe.nCols) {
				if (nRows == pe.nRows) {
					if (stratumDesign.size() == pe.stratumDesign.size()) {
						for (String stratumName : stratumDesign.keySet()) {
							if (pe.stratumDesign.get(stratumName) == null) {
								return false;
							} else {
								FinitePopulationEstimate thisPopTotEstimate = stratumDesign.get(stratumName);
								FinitePopulationEstimate thatPopTotEstimate = pe.stratumDesign.get(stratumName);
								if (!thisPopTotEstimate.isMergeableEstimate(thatPopTotEstimate)) {
									return false;
								}
							}
						}
						return true;
					}
				}
			}
		}
		return false;
	}

	@Override
	protected StratifiedPopulationEstimate add(PointEstimate pointEstimate) {
		if (isMergeableEstimate(pointEstimate)) {
			StratifiedPopulationEstimate newEstimate = getEmptyEstimate();
			StratifiedPopulationEstimate thatEstimate = (StratifiedPopulationEstimate) pointEstimate;
			
			for (String stratumName : stratumDesign.keySet()) {
				FinitePopulationEstimate thisPopTotEstimate = stratumDesign.get(stratumName);
				FinitePopulationEstimate thatPopTotEstimate = thatEstimate.stratumDesign.get(stratumName);
				newEstimate.stratumDesign.put(stratumName, (FinitePopulationEstimate) thisPopTotEstimate.add(thatPopTotEstimate));
			}
			return newEstimate;
		} else {
			throw new InvalidParameterException("Incompatible point estimates!");
		}
	}

	@Override
	protected StratifiedPopulationEstimate subtract(PointEstimate pointEstimate) {
		if (isMergeableEstimate(pointEstimate)) {
			StratifiedPopulationEstimate newEstimate = getEmptyEstimate();
			StratifiedPopulationEstimate thatEstimate = (StratifiedPopulationEstimate) pointEstimate;
			
			for (String stratumName : stratumDesign.keySet()) {
				FinitePopulationEstimate thisPopTotEstimate = stratumDesign.get(stratumName);
				FinitePopulationEstimate thatPopTotEstimate = thatEstimate.stratumDesign.get(stratumName);
				newEstimate.stratumDesign.put(stratumName, (FinitePopulationEstimate) thisPopTotEstimate.subtract(thatPopTotEstimate));
			}
			return newEstimate;
		} else {
			throw new InvalidParameterException("Incompatible point estimates!");
		}
	}

	@Override
	protected StratifiedPopulationEstimate multiply(double scalar) {
		StratifiedPopulationEstimate newEstimate = getEmptyEstimate();
	
		for (String stratumName : stratumDesign.keySet()) {
			FinitePopulationEstimate thisPopTotEstimate = stratumDesign.get(stratumName);
			newEstimate.stratumDesign.put(stratumName, (FinitePopulationEstimate) thisPopTotEstimate.multiply(scalar));
		}
		return newEstimate;
	}

	@Override
	public int getSampleSize() {
		int sampleSize = 0;
		for (FinitePopulationEstimate pte : stratumDesign.values()) {
			sampleSize += pte.getSampleSize();
		}
		return sampleSize;
	}
	
	@Override
	protected final Map<String, Matrix> getObservations() {
		Map<String, Matrix> outputMap = new HashMap<String, Matrix>();
		for (FinitePopulationEstimate pte : stratumDesign.values()) {
			outputMap.putAll(pte.getObservations());
		}
		return outputMap;
	}

	@Override
	protected final List<String> getPopulationUnitIds() {
		List<String> puIds = new ArrayList<String>();
		for (FinitePopulationEstimate pte : stratumDesign.values()) {
			puIds.addAll(pte.getPopulationUnitIds());
		}
		return puIds;
	}
	
	final Map<String, String> getPopulationUnitToStrataMap() {
		Map<String, String> outputMap = new HashMap<String, String>();
		for (String stratumName : stratumDesign.keySet()) {
			FinitePopulationEstimate pte = stratumDesign.get(stratumName);
			for (String puId : pte.getPopulationUnitIds()) {
				outputMap.put(puId, stratumName);
			}
		}
		return outputMap;
	}

	@Override
	StratifiedPopulationEstimate getEmptyEstimate() {
		List<String> stratumNames = new ArrayList<String>();
		List<Double> stratumSizes = new ArrayList<Double>();
		for (String n : stratumDesign.keySet()) {
			stratumNames.add(n);
			stratumSizes.add(stratumDesign.get(n).getPopulationSize());
		}
		return new StratifiedPopulationEstimate(stratumNames, stratumSizes);
	}
	
}
