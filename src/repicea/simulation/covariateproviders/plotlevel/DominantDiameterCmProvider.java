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
 * This interface ensures the stand instance can provide its own
 * dominant diameter, ie. usually the average diameter of the 100 
 * thickest trees per hectare.
 * @author Mathieu Fortin - November 2012
 */
public interface DominantDiameterCmProvider {

	/**
	 * This method returns the stand dominant diameter (cm)
	 * @return a double
	 */
	public double getDominantDiameterCm();
}
