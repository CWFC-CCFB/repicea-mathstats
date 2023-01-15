/*
 * This file is part of the repicea-simulation library.
 *
 * Copyright (C) 2009-2014 Mathieu Fortin for Rouge-Epicea
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
package repicea.simulation.allometrycalculator;

import repicea.simulation.covariateproviders.treelevel.BasalAreaM2Provider;
import repicea.simulation.covariateproviders.treelevel.DbhCmProvider;
import repicea.simulation.covariateproviders.treelevel.ExpansionFactorProvider;
import repicea.simulation.covariateproviders.treelevel.SquaredDbhCmProvider;

/**
 * The LightAllometryCalculableTree interface implements the basic methods for calculating the number of stems per hectare,
 * the basal area per hectare and the mean quadratic diameter.
 * @author Mathieu Fortin - July 2014
 */
public interface LightAllometryCalculableTree extends DbhCmProvider,
														SquaredDbhCmProvider,
														ExpansionFactorProvider,
														BasalAreaM2Provider {

}
