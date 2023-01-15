/*
 * This file is part of the repicea library.
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
package repicea.simulation.covariateproviders.treelevel;

/**
 * This interface ensures that the tree instance can provide the standard deviation of its dbh. 
 * This interface is useful if the tree instance represents many trees.
 * @author Mathieu Fortin - October 2018
 */
public interface DbhCmStandardDeviationProvider {

	/**
	 * This method returns the standard deviation of the DBH (cm)
	 * @return a double
	 */
	public double getDbhCmStandardDeviation();
}
