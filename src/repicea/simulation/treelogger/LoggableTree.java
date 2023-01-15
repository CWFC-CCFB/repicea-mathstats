/*
 * This file is part of the repicea-simulation library.
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
package repicea.simulation.treelogger;

import repicea.simulation.covariateproviders.treelevel.BarkProportionProvider;
import repicea.simulation.covariateproviders.treelevel.CommercialVolumeM3Provider;
import repicea.simulation.covariateproviders.treelevel.ExpansionFactorProvider;
import repicea.simulation.covariateproviders.treelevel.SpeciesNameProvider;
import repicea.simulation.covariateproviders.treelevel.TreeWeightProvider;

/**
 * This interface applies at the tree level and only serves to determine the
 * object that are compatible with the tree logger in the repicea-foresttools
 * library.
 * 
 * @author Mathieu Fortin and Jean-Francois Lavoie - January 2012
 */
public interface LoggableTree extends CommercialVolumeM3Provider,
										SpeciesNameProvider, 
										ExpansionFactorProvider,
										TreeWeightProvider,
										BarkProportionProvider {

}
