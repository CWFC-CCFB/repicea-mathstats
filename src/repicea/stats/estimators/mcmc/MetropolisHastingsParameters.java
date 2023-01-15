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

public class MetropolisHastingsParameters {

	public int nbBurnIn = 10000;
	public int nbRealizations = 500000 + nbBurnIn;
	public int nbInternalIter = 100000;
	public int oneEach = 50;
	public int nbInitialGrid = 10000;	

	public MetropolisHastingsParameters() {}

	@Override
	public MetropolisHastingsParameters clone() {
		try {
			return (MetropolisHastingsParameters) super.clone();
		} catch (CloneNotSupportedException e) {
			e.printStackTrace();
			return null;
		}
	}

}
