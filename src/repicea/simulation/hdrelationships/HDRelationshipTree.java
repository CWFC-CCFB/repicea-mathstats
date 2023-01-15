/*
 * This file is part of the repicea-simulation library.
 *
 * Copyright (C) 2009-2015 Mathieu Fortin for Rouge-Epicea
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
package repicea.simulation.hdrelationships;

import repicea.simulation.HierarchicalLevel;
import repicea.simulation.MonteCarloSimulationCompliantObject;
import repicea.simulation.covariateproviders.treelevel.HeightMProvider;
import repicea.stats.distributions.GaussianErrorTermList.IndexableErrorTerm;

public abstract interface HDRelationshipTree extends MonteCarloSimulationCompliantObject, 
													IndexableErrorTerm, 
													HeightMProvider {

	@Override
	default public HierarchicalLevel getHierarchicalLevel() {return HierarchicalLevel.TREE;}

	/**
	 * This method returns the error group in case of different error correlation structure. For instance, if coniferous species 
	 * have a correlation structure that differs from that of broadleaved species.
	 * @return an Enum that defines the group
	 */
	public Enum<?> getHDRelationshipTreeErrorGroup();
	
}
