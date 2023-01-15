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
package repicea.simulation.disturbances;

import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import repicea.simulation.MonteCarloSimulationCompliantObject;

public class DisturbanceOccurrences {
	
	protected final List<Integer> datesYr;
	protected final DisturbanceAffectedProvider parms;
	protected final Map<String, Boolean> outcomeMap;
	
	/**
	 * Constructor for multiple occurrences
	 * @param parms
	 * @param datesYr a List of dates (integers)
	 */
	public DisturbanceOccurrences(DisturbanceAffectedProvider parms, List<Integer> datesYr) {
		this(parms);
		for (Integer dateYr : datesYr) {
			addOccurrenceDateYr(dateYr);
		}
	}

	
	/**
	 * Basic constructor
	 * @param parms
	 */
	public DisturbanceOccurrences(DisturbanceAffectedProvider parms) {
		this.parms = parms;
		this.datesYr = new ArrayList<Integer>();
		this.outcomeMap = new HashMap<String, Boolean>();
	}
	
	/**
	 * Constructor for single occurrences
	 * @param parms
	 * @param dateYr
	 */
	public DisturbanceOccurrences(DisturbanceAffectedProvider parms, int dateYr) {
		this(parms);
		addOccurrenceDateYr(dateYr);
	}
	
	
	/**
	 * This method checks whether a particular plot is affected by the disturbance. It assumes
	 * that the method isThereADisturbance of the DisturbanceParameters class has returned true. 
	 * @param plot a MonteCarloSimulationCompliantObject instance
	 * @return a boolean
	 */
	public synchronized boolean isThisPlotAffected(MonteCarloSimulationCompliantObject plot) {
		if (!outcomeMap.containsKey(plot.getSubjectId())) {
			boolean outcome;
			if (parms == null) { // former implementation
				outcome = true;
			} else {
				outcome = parms.isThisPlotAffected(plot);
			}
			outcomeMap.put(plot.getSubjectId(), outcome);
		}
		return outcomeMap.get(plot.getSubjectId());
	}

	
	public void addOccurrenceDateYr(int dateYr) {
		datesYr.add(dateYr);
		Collections.sort(datesYr); // to have them is ascending order.
	}
	
	/**
	 * Return the date of the latest occurrence prior or concurrent to the currentDateYr.
	 * @param currentDateYr the reference date (yr)
	 * @return an integer which is the date or -1 if there is no occurrence
	 */
	public int getLastOccurrenceDateYrToDate(int currentDateYr) {
		int latestOccurrenceDate = -1;
		for (Integer o : datesYr) {
			if (o <= currentDateYr) {
				if (latestOccurrenceDate == -1 || o > latestOccurrenceDate) {
					latestOccurrenceDate = o;
				}
			} else {
				break;	// means we are already beyond the currentDateYr parameter and there is no chance of having a date prior or concurrent to currentDateYr
			}
		}
		return latestOccurrenceDate;
	}

	
	/**
	 * Return the date of the latest occurrence.
	 * @return an integer which is the date or -1 if there is no occurrence
	 */
	public int getLastOccurrenceDateYr() {
		if (!datesYr.isEmpty()) {
			return datesYr.get(datesYr.size() - 1);
		} else {
			return -1;
		}
	}
	
	
	/**
	 * Returns true if there is at least one occurrence. 
	 * @return a boolean
	 */
	public boolean isThereAnyOccurrence() {
		return !datesYr.isEmpty();
	}

	/**
	 * Return the number of occurrences.
	 * @return an integer
	 */
	public int getNumberOfOccurrences() {
		return datesYr.size();
	}
}
