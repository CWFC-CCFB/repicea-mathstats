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
import java.util.Collections;
import java.util.List;
import java.util.Map;
import java.util.concurrent.ConcurrentHashMap;

import repicea.math.Matrix;
import repicea.stats.estimates.Estimate;

@SuppressWarnings("serial")
public abstract class AbstractSimplePointEstimate extends AbstractPointEstimate {

	private final Map<String, Matrix> observations;
	protected boolean needsToBeRecalculated;
	
	/**
	 * Basic constructor without population size.
	 */
	protected AbstractSimplePointEstimate() {
		super();
		observations = new ConcurrentHashMap<String, Matrix>();
		needsToBeRecalculated = true;
	}

	/**
	 * Create a List with the ordered sample plot ids.
	 * @return a List of strings
	 */
	@Override
	protected final List<String> getPopulationUnitIds() {
		List<String> puIds = new ArrayList<String>();
		puIds.addAll(observations.keySet());
		Collections.sort(puIds);
		return puIds;
	}
	
	/**
	 * Validate the population unit before adding it to the
	 * observations.
	 * @param obs a Matrix instance instance
	 * @param obsId the observation id
	 * @param stratumName useless argument for this class. Can be set to null.
	 */
	@Override
	protected final void validateUnit(Matrix obs, String obsId, String stratumName) {
		if (observations.containsKey(obsId)) {
			throw new InvalidParameterException("The sample id " + obsId + " is already contained in the observation map!");
		}
	}

	/**
	 * Add an observation to the sample.
	 * 
	 * @param obs a Matrix instance instance
	 * @param obsId the observation id
	 * @param stratumName useless argument for this class. Can be set to null.
	 */
	@Override
	public void addObservation(Matrix obs, String obsId, String stratumName) {
		super.addObservation(obs, obsId, stratumName);
		observations.put(obsId, obs);
		needsToBeRecalculated = true;
	}

	/**
	 * Add an observation to the sample.<p>
	 * 
	 * This method call the {@link AbstractSimplePointEstimate#addObservation(Matrix, String, String)}
	 * method with the argument stratumName set to null.
	 * 
	 * @param obs a Matrix instance instance
	 * @param obsId the observation id
	 */
	public final void addObservation(Matrix obs, String obsId) {
		addObservation(obs, obsId, null);
	}

	@Override
	protected final Map<String, Matrix> getObservations() {return observations;}

	@Override
	protected boolean isMergeableEstimate(Estimate<?,?,?> estimate) {
		if (getClass().equals(estimate.getClass())) {
			AbstractSimplePointEstimate pe = (AbstractSimplePointEstimate) estimate;
			if (getPopulationUnitIds().equals(pe.getPopulationUnitIds())) {	// make sure we have the same sample ids
				if (nRows == pe.nRows) {
					if (nCols == pe.nCols) {
						return true;
					}
				}
			}
		}
		return false;
	}

	@Override
	protected abstract AbstractSimplePointEstimate add(PointEstimate pointEstimate);

	@Override
	protected abstract AbstractSimplePointEstimate subtract(PointEstimate pointEstimate);

	@Override
	protected abstract AbstractSimplePointEstimate multiply(double scalar);

	@Override
	public final int getSampleSize() {return observations.size();}
}
