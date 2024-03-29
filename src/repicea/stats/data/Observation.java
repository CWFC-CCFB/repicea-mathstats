/*
 * This file is part of the repicea library.
 *
 * Copyright (C) 2009-2012 Mathieu Fortin for Rouge-Epicea
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
package repicea.stats.data;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.LinkedHashMap;
import java.util.List;

import com.fasterxml.jackson.annotation.JsonProperty;

import repicea.data.Record;

@SuppressWarnings("rawtypes")
public class Observation implements Record, Comparable {

	static List<Integer> ComparableFields = new ArrayList<Integer>();
	
	List<Object> values;
	
	protected Observation(Object[] obj) {
		values = new ArrayList<Object>();
		values.addAll(Arrays.asList(obj));
	}
		
	@SuppressWarnings("unchecked")
	@Override
	public int compareTo(Object o) {
		for (Integer index : ComparableFields) {
			Comparable thisValue = (Comparable) values.get(index);
			Comparable thatValue = (Comparable) ((Observation) o).values.get(index);
			int comparisonResult = thisValue.compareTo(thatValue);
			if (comparisonResult < 0) {
				return -1;
			} else if (comparisonResult > 0) {
				return 1;
			}
		}
		return 0;
	}

	/**
	 * Converts this observation to an array of Object instances
	 * @return an Array of Object instances
	 */
	public Object[] toArray() {return values.toArray();}


	/**
	 * Checks if two observations have the same values.
	 * @param obs an Observation instance
	 * @return a boolean
	 */
	public boolean isEqualToThisObservation(Observation obs) {
		if (obs == null) {
			return false;
		} else {
			if (values.size() != obs.values.size()) {
				return false;
			} 
			for (int i = 0; i < values.size(); i++) {
				Object thisValue = values.get(i);
				Object thatValue = obs.values.get(i);
				Class thisClass = thisValue.getClass();
				if (!thisClass.equals(thatValue.getClass())) {
					return false;
 				} else {
 					if (thisClass.equals(Double.class)) {
 						if (Math.abs((Double) thisValue - (Double) thatValue) > 1E-8) {
 							return false;
 						}
 					} else if (!thisValue.equals(thatValue)) {
 						return false;
 					}
 				}
			}
			return true;
		}
	}
	
	/**
	 * Return the value at a particular index within this observation instance.
	 * @param index the index an integer.
	 * @return an Object instance
	 */
	public Object getValueAt(int index) {
		return values.get(index);
	}

	@JsonProperty("values")
	@Override
	public Object[] getContentArray() {
		return toArray();
	}
	
}
