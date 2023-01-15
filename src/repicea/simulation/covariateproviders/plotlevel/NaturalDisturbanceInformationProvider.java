/*
 * This file is part of the repicea library.
 *
 * Copyright (C) 2009-2019 Mathieu Fortin for Rouge-Epicea
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

import repicea.simulation.disturbances.DisturbanceTypeProvider.DisturbanceType;

/**
 * This interface ensures the instance can provide either the time since the last
 * disturbance or at least the time since the first known date.
 * @author Mathieu Fortin - March 2019
 */
public interface NaturalDisturbanceInformationProvider {

	/**
	 * This method returns the time since the last disturbance or null if this time is unknown.
	 * @param currentDateYrs the current date
	 * @return an Integer instance or null
	 */
	public Integer getTimeSinceLastDisturbanceYrs(DisturbanceType type, int currentDateYrs);
	
	/**
	 * This method returns the time since the first known date. This is a work-around if the 
	 * getTimeSinceLastDisturbanceYrs method returns null.
	 * @param currentDateYrs the current date
	 * @return a positive integer
	 */
	public int getTimeSinceFirstKnownDateYrs(int currentDateYrs);
	
}
