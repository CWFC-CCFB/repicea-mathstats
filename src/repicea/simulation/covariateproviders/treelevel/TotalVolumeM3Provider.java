/*
 * This file is part of the repicea library.
 *
 * Copyright (C) 2009-2020 Mathieu Fortin for Rouge-Epicea
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
 * This interface ensures the tree instance can provide its
 * total volume. This volume can be either over or under bark. 
 * The default method assumes it is under bark.
 * @author Mathieu Fortin - August 2020
 */
public interface TotalVolumeM3Provider {

	/**
	 * This method calculates the total volume of the tree.
	 * @return a double - the volume (m3)
	 */
	public double getTotalVolumeM3();

	public boolean isTotalVolumeOverbark();

}
