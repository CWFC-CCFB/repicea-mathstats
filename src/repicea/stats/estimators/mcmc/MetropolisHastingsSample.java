/*
 * This file is part of the repicea library.
 *
 * Copyright (C) 2009-2021 Mathieu Fortin for Rouge Epicea.
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
package repicea.stats.estimators.mcmc;

import repicea.math.Matrix;

class MetropolisHastingsSample implements Comparable<MetropolisHastingsSample> {

	final Matrix parms;
	final double llk;
	
	MetropolisHastingsSample(Matrix parms, double lk) {
		this.parms = parms;
		this.llk = lk;
	}

	@Override
	public int compareTo(MetropolisHastingsSample arg0) {
		if (llk > arg0.llk) {
			return 1;
		} else if (llk < arg0.llk) {
			return -1;
		} else {
			return 0;
		}
	}

	@Override
	public String toString() {
		return "LLK=" + llk + ", " + parms.toString();
	}
	
}
