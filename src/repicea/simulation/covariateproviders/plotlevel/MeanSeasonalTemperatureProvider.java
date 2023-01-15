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

/**
 * This interface ensures the plot or the stand instance can provide its seasonal temperature over the interval.
 * @author Mathieu Fortin - July 2019
 *
 */
public interface MeanSeasonalTemperatureProvider {

	
	/**
	 * Returns the monthly mean temperature above 6 degrees Celsius.
	 * @return a double
	 */
	public double getMeanAnnualTempAbove6C();

}
