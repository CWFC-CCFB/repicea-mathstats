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
 * This class implements the sampling with uneven inclusion probabilities across the units.
 * 
 * @author Mathieu Fortin - May 2018
 */
public class PopulationUnitWithUnequalInclusionProbability extends PopulationUnit {

	private final double inclusionProbability;
	
	/**
	 * Constructor. IMPORTANT: the inclusion probability is the individual inclusion probability 
	 * of this unit. For instance, that would be the area of the plot divided by the total area, 
	 * and not the plot area multiplied by the sample size. This product is actually handled internally.
	 * @param obs the response that was observed in the unit (a column vector)
	 * @param sampleId a string that stands for the sample id
	 * @param inclusionProbability the probability that this unit is part of the sample.
	 */
	public PopulationUnitWithUnequalInclusionProbability(String sampleId, Matrix obs, double inclusionProbability) {
		super(sampleId, obs);
		if (inclusionProbability <= 0 || inclusionProbability >= 1d) {
			throw new InvalidParameterException("The inclusion probability must be larger than 0 and smaller than 1");
		}
		this.inclusionProbability = inclusionProbability;
	}

	/**
	 * This method returns the sampling probability of this population unit.
	 * @return a double
	 */
	public double getInclusionProbability() {return inclusionProbability;}
}
