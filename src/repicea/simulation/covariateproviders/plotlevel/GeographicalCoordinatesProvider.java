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
package repicea.simulation.covariateproviders.plotlevel;

/**
 * This interface ensures the plot instance can provide its geographical location.
 * @author Mathieu Fortin - June 2014
 */
public interface GeographicalCoordinatesProvider extends ElevationMProvider {

	/**
	 * This method returns the latitude of the plot in degrees.
	 * @return a double
	 */
	public double getLatitudeDeg();
	
	/**
	 * This method returns the longitude of the plot in degrees.
	 * @return a double
	 */
	public double getLongitudeDeg();
		
	
}
