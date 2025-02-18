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

import repicea.math.Matrix;

/**
 * This class is the most basic population unit that can be sampled. It assumes that
 * the inclusion probability is even across the units.
 * 
 * @author Mathieu Fortin - May 2018
 * @deprecated Use the PopulationUnit class instead of this class.
 */
@Deprecated
class PopulationUnitWithEqualInclusionProbability extends PopulationUnit {

	/**
	 * Constructor.
	 * @param sampleId a string that stands for the sample id
	 * @param obs a Matrix instance. Must be a column vector
	 */
	private PopulationUnitWithEqualInclusionProbability(String sampleId, Matrix obs) {
		super(sampleId, obs);
	}

}
