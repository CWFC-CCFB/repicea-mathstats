/*
 * This file is part of the repicea library.
 *
 * Copyright (C) 2009-2013 Mathieu Fortin for Rouge-Epicea
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
package repicea.simulation.covariateproviders.plotlevel;

import java.util.Collection;

import repicea.simulation.covariateproviders.treelevel.TreeStatusProvider.StatusClass;

/**
 * This interface ensures the stand can provide collections of trees according
 * to a particular tree status.
 * @author Mathieu Fortin - September 2013
 */
public interface TreeStatusCollectionsProvider {

	/**
	 * This method returns the Collection of Tree instances with a particular 
	 * status class.
	 * @param statusClass a StatusClass enum
	 * @return a Collection of trees which can be empty
	 */
	@SuppressWarnings("rawtypes")
	public Collection getTrees(StatusClass statusClass);

}
