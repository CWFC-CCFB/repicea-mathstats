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
package repicea.simulation.climate;

import java.security.InvalidParameterException;
import java.util.ArrayList;

import repicea.simulation.climate.REpiceaClimateVariableMap.ClimateVariable;

public class REpiceaClimateChangeTrend extends ArrayList<REpiceaClimateChangeTrendSegment> {

	public void addSegment(int startDateYr, int endDateYr, REpiceaClimateVariableChangeMap changeMap) {
		if (startDateYr >= endDateYr) {
			throw new InvalidParameterException("The start date must be prior to the end date!");
		}
		if (isEmpty() || getEndDateYr() ==  startDateYr) {
			add(new REpiceaClimateChangeTrendSegment(startDateYr, endDateYr, changeMap));
		} else {
			throw new InvalidParameterException("There is a mismatch between the start date and the end date of the previous segment!");
		}
	}

	/**
	 * Returns the last end date of the trend.
	 * @return an integer
	 */
	public int getEndDateYr() {
		if (isEmpty()) {
			return -1;
		} else {
			return get(size()-1).endDateYr;
		}
	}
	
	/**
	 * Returns the start date of the trend.
	 * @return an integer
	 */
	public int getStartDateYr() {
		if (isEmpty()) {
			return -1;
		} else {
			return get(0).startDateYr;
		}
	}
	
	/**
	 * Provides the change in the variable between the start and the end dates.
	 * @param startDateYr
	 * @param endDateYr
	 * @param variable a ClimateVariable enum
	 * @return the change
	 */
	public double getChangeFromTo(int startDateYr, int endDateYr, ClimateVariable variable) {
		double increase = 0d;
		for (REpiceaClimateChangeTrendSegment segment : this) {
			int nbYears = 0;
			if (startDateYr <= segment.startDateYr && endDateYr > segment.startDateYr) { // then this segment in part or in whole should be considered
				nbYears = endDateYr > segment.endDateYr ? segment.endDateYr - segment.startDateYr : endDateYr - segment.startDateYr;
			} else if (startDateYr <= segment.endDateYr && endDateYr > segment.startDateYr) {
				nbYears = endDateYr > segment.endDateYr ? segment.endDateYr - startDateYr : endDateYr - startDateYr;
			}
			if (segment.changeMap.containsKey(variable)) {
				increase += nbYears * segment.changeMap.get(variable);
			}
		}
		return increase;
	}
	
	/**
	 * Return the average change over a particular time interval with respect to a reference date.
	 * @param referenceDateYr
	 * @param startDateYr
	 * @param endDateYr
	 * @param variable a ClimateVariable enum
	 * @return the average change
	 */
	public double getAverageChangeOverThisPeriod(int referenceDateYr, int startDateYr, int endDateYr, ClimateVariable variable) {
		if (startDateYr >= endDateYr) {
			throw new InvalidParameterException("The start date must be prior to the end date!");
		}
		double area = 0d;
		for (REpiceaClimateChangeTrendSegment segment : this) {
			int nbYears = 0;
			double h1,h2;
			if (startDateYr <= segment.startDateYr && endDateYr > segment.startDateYr) { // then this segment in part or in whole should be considered
				int endYear = endDateYr > segment.endDateYr ? segment.endDateYr : endDateYr;
				nbYears = endYear - segment.startDateYr;
				h1 = getChangeFromTo(referenceDateYr, segment.startDateYr, variable);
				h2 = h1 + getChangeFromTo(segment.startDateYr, endYear, variable);
				area += (h1 + h2) * .5 * nbYears;
			} else if (startDateYr <= segment.endDateYr && endDateYr > segment.startDateYr) {
				int endYear = endDateYr > segment.endDateYr ? segment.endDateYr : endDateYr;
				nbYears = endYear - startDateYr;
				h1 = getChangeFromTo(referenceDateYr, startDateYr, variable);
				h2 = h1 + getChangeFromTo(startDateYr, endYear, variable);
				area += (h1 + h2) * .5 * nbYears;
			} 
		}
		return area / (endDateYr - startDateYr);
	}
	
	
}
