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
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.concurrent.ConcurrentHashMap;

import repicea.math.Matrix;
import repicea.math.SymmetricMatrix;
import repicea.stats.sampling.PopulationUnit;

/**
 * A class to implement stratified sampling design.<p>
 * The point estimate is that of the total.
 * @author Mathieu Fortin - January 2025
 */
@SuppressWarnings("serial")
public class StratifiedPopulationTotalEstimate extends AbstractPointEstimate {

	private final Map<String, PopulationTotalEstimate> stratumDesign;
	
	
	private StratifiedPopulationTotalEstimate() {
		super();
		stratumDesign = new ConcurrentHashMap<String, PopulationTotalEstimate>();
	}
	
	/**
	 * Constructor.
	 * @param stratumNames a List of strings that stand for the stratum names.
	 * @param strataPopulationSizes a List of Double that stand for the stratum sample sizes.
	 */
	public StratifiedPopulationTotalEstimate(List<String> stratumNames, List<Double> strataPopulationSizes) {
		this();
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
		if (pu == null) {
			throw new InvalidParameterException("The pu argument must be non null!");
		}
		if (nCols == 0) {
			nCols = pu.getData().m_iCols;
		}
		if (nRows == 0) {
			nRows = pu.getData().m_iRows;
		}
		if (pu.getData().m_iCols != nCols || pu.getData().m_iRows != nRows) {
			throw new InvalidParameterException("The observation is incompatible with what was already observed!");
		} else {
			stratumDesign.get(stratumName).addObservation(pu);
		}
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

	@Override
	protected boolean isMergeableEstimate(Estimate<?, ?, ?> estimate) {
		if (getClass().equals(estimate.getClass())) {
			StratifiedPopulationTotalEstimate pe = (StratifiedPopulationTotalEstimate) estimate;
			if (nCols == pe.nCols) {
				if (nRows == pe.nRows) {
					if (stratumDesign.size() == pe.stratumDesign.size()) {
						for (String stratumName : stratumDesign.keySet()) {
							if (pe.stratumDesign.get(stratumName) == null) {
								return false;
							} else {
								PopulationTotalEstimate thisPopTotEstimate = stratumDesign.get(stratumName);
								PopulationTotalEstimate thatPopTotEstimate = pe.stratumDesign.get(stratumName);
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
	protected StratifiedPopulationTotalEstimate add(PointEstimate pointEstimate) {
		if (isMergeableEstimate(pointEstimate)) {
			StratifiedPopulationTotalEstimate newEstimate = new StratifiedPopulationTotalEstimate();
			StratifiedPopulationTotalEstimate thatEstimate = (StratifiedPopulationTotalEstimate) pointEstimate;
			
			for (String stratumName : stratumDesign.keySet()) {
				PopulationTotalEstimate thisPopTotEstimate = stratumDesign.get(stratumName);
				PopulationTotalEstimate thatPopTotEstimate = thatEstimate.stratumDesign.get(stratumName);
				newEstimate.stratumDesign.put(stratumName, thisPopTotEstimate.add(thatPopTotEstimate));
			}
			return newEstimate;
		} else {
			throw new InvalidParameterException("Incompatible point estimates!");
		}
	}

	@Override
	protected StratifiedPopulationTotalEstimate subtract(PointEstimate pointEstimate) {
		if (isMergeableEstimate(pointEstimate)) {
			StratifiedPopulationTotalEstimate newEstimate = new StratifiedPopulationTotalEstimate();
			StratifiedPopulationTotalEstimate thatEstimate = (StratifiedPopulationTotalEstimate) pointEstimate;
			
			for (String stratumName : stratumDesign.keySet()) {
				PopulationTotalEstimate thisPopTotEstimate = stratumDesign.get(stratumName);
				PopulationTotalEstimate thatPopTotEstimate = thatEstimate.stratumDesign.get(stratumName);
				newEstimate.stratumDesign.put(stratumName, thisPopTotEstimate.subtract(thatPopTotEstimate));
			}
			return newEstimate;
		} else {
			throw new InvalidParameterException("Incompatible point estimates!");
		}
	}

	@Override
	protected StratifiedPopulationTotalEstimate multiply(double scalar) {
		StratifiedPopulationTotalEstimate newEstimate = new StratifiedPopulationTotalEstimate();
	
		for (String stratumName : stratumDesign.keySet()) {
			PopulationTotalEstimate thisPopTotEstimate = stratumDesign.get(stratumName);
			newEstimate.stratumDesign.put(stratumName, thisPopTotEstimate.multiply(scalar));
		}
		return newEstimate;
	}

	@Override
	public int getSampleSize() {
		int sampleSize = 0;
		for (PopulationTotalEstimate pte : stratumDesign.values()) {
			sampleSize += pte.getSampleSize();
		}
		return sampleSize;
	}
	
	@Override
	protected final Map<String, PopulationUnit> getObservations() {
		Map<String, PopulationUnit> outputMap = new HashMap<String, PopulationUnit>();
		for (PopulationTotalEstimate pte : stratumDesign.values()) {
			outputMap.putAll(pte.getObservations());
		}
		return outputMap;
	}

	@Override
	protected final List<String> getPopulationUnitIds() {
		List<String> puIds = new ArrayList<String>();
		for (PopulationTotalEstimate pte : stratumDesign.values()) {
			puIds.addAll(pte.getPopulationUnitIds());
		}
		return puIds;
	}
	
	final Map<String, String> getPopulationUnitToStrataMap() {
		Map<String, String> outputMap = new HashMap<String, String>();
		for (String stratumName : stratumDesign.keySet()) {
			PopulationTotalEstimate pte = stratumDesign.get(stratumName);
			for (String puId : pte.getPopulationUnitIds()) {
				outputMap.put(puId, stratumName);
			}
		}
		return outputMap;
	}
	
}
