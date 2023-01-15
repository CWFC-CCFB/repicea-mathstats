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
 * This interface ensures the tree instance can provide its commercial volume. 
 * This volume can be either under or over bark. 
 * @author Mathieu Fortin - August 2020
 *
 */
public interface CommercialVolumeM3Provider {

	/**
	 * This method calculates the commercial volume, that is the volume from
	 * the stump up to a small-end diameter that defines the commercial limit.
	 * The volume is calculated WITHOUT any expansion factor, that is the volume of
	 * a single tree, regardless of the expansion factor.
	 * @return a double -- the volume (m3)
	 */
	public double getCommercialVolumeM3();

	public boolean isCommercialVolumeOverbark();

}
