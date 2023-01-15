/*
 * This file is part of the repicea-statistics library.
 *
 * Copyright (C) 2009-2018 Mathieu Fortin for Rouge-Epicea
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

import repicea.math.Matrix;

/**
 * This class represents a single sample unit that would
 * have been observed in a finite population.
 *
 * @author Mathieu Fortin - May 2018
 */
public abstract class PopulationUnit {

	private final Matrix data;
	private final String sampleId;
	
	/**
	 * Constructor.
	 * @param sampleId a string that stands for the sample id
	 * @param obs a Matrix instance. Must be a column vector
	 */
	protected PopulationUnit(String sampleId, Matrix obs) {
		if (sampleId == null) {
			throw new InvalidParameterException("The sample argument must be non null!");
		}
		if (!obs.isColumnVector()) {
			throw new InvalidParameterException("The obs parameter should be a column vector!");
		}
		this.sampleId = sampleId;
		this.data = obs;
	}

	/**
	 * Return the response that was observed in 
	 * this unit of the population. The response can be multivariate.
	 * @return a Matrix instance
	 */
	public Matrix getData() {return data;}

	/**
	 * Provide the sample id which is a unique string for this sample unit.
	 * @return a String
	 */
	public String getSampleId() {return sampleId;}
}
