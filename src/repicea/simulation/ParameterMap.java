/*
 * This file is part of the repicea-foresttools library.
 *
 * Copyright (C) 2009-2014 Mathieu Fortin for Rouge-Epicea
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

import java.io.Serializable;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.TreeMap;

import repicea.math.Matrix;

@SuppressWarnings("serial")
public class ParameterMap implements Serializable {

	public static int DEFAULT_INDEX = 0;

	private final Map<Integer, Map<Integer, List<Double>>> parametersMap;
	
	/**
	 * Constructor.
	 */
	public ParameterMap() {
		parametersMap = new HashMap<Integer, Map<Integer, List<Double>>>();
	} 

	/**
	 * This method returns a Matrix instance from the list of double that were stored in the map according to the 
	 * two indices.
	 * @param index1 first index
	 * @param index2 second index
	 * @return a Matrix instance
	 */
	public Matrix get(int index1, int index2) {
		Map<Integer, List<Double>> subMap = parametersMap.get(index1);
		if (subMap == null) {
			return null;
		} else {
			List<Double> oTempVector = subMap.get(index2);
			if (oTempVector == null) {
				return null;
			} else {
				Matrix mat = new Matrix(oTempVector);
				return mat;
			}
		}
	}

	/**
	 * This method returns the same than get(index1, index2) in case of a single index implementation.
	 * @param index2 the index
	 * @return a Matrix instance
	 */
	public Matrix get(int index2) {
		return get(DEFAULT_INDEX, index2);
	}

	/**
	 * This method returns the same than get(index1, index2) in case of no index implementation.
	 * @return a Matrix instance
	 */
	public Matrix get() {
		return get(DEFAULT_INDEX, DEFAULT_INDEX);
	}

	/**
	 * This method stores a double in the map at the specified indices.
	 * @param index1 the first index
	 * @param index2 the second index
	 * @param value a double 
	 */
	protected void addParameter(int index1, int index2, double value) {
		Map<Integer, List<Double>> oMap = parametersMap.get(index1);
		if (oMap == null) {
			oMap = new TreeMap<Integer, List<Double>>();
			parametersMap.put(index1, oMap);
		} 
		List<Double> vector = oMap.get(index2);
		if (vector == null) {
			vector = new ArrayList<Double>();
			oMap.put(index2, vector);
		}
		vector.add(value);
	}


	/**
	 * This method does the same than addParameter(index1, index2, vec) except it is for a one-index implementation.
	 * @param index2 the index 
	 * @param value a double 
	 */
	protected void addParameter(int index2, double value) {
		addParameter(DEFAULT_INDEX, index2, value);
	}

	/**
	 * This method does the same than addParameter(index1, index2, vec) except it is for a no-index implementation.
	 * @param value a double 
	 */
	protected void addParameter(double value) {
		addParameter(DEFAULT_INDEX, DEFAULT_INDEX, value);
	}


	
}
