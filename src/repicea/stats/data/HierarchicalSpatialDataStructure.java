/*
 * This file is part of the repicea-statistics library.
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

import java.util.List;
import java.util.Map;

import repicea.math.Matrix;


public interface HierarchicalSpatialDataStructure extends HierarchicalStatisticalDataStructure {
	
	/**
	 * Return the distance map between the observations. <br>
	 * <br>
	 * The elements of the list correspond to distance types. The first key represents 
	 * the observation index. The second key represent the distant observation index, 
	 * whereas the value is the the distance to this distant observation.
	 * @return a Map instance
	 */
	public List<Map<Integer, Map<Integer, Double>>> getDistancesBetweenObservations();
	
//	/**
//	 * This method returns the angle map between the observations. The keys represent the observation index, whereas the double 
//	 * is the angle in radians. 
//	 * @return a Map instance
//	 */
//	public Map<Integer, Map<Integer, Matrix>> getAngleBetweenObservations();
	
	/**
	 * This method sets the distance field names. <br>
	 * <br> 
	 * The number of elements in the first list sets the different types of distance. For instance, we might want to
	 * account for two types of distances (e.g. spatial distance and temporal distance). The inner list contains the fields to 
	 * compute the distance.
	 * @param fields a List of List that contains the names of the fields that serve to compute the Euclidian distance.
	 */
	public void setDistanceFields(List<List<String>> fields);
	
}
