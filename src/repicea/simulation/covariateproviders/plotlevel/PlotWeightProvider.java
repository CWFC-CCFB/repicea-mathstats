/*
 * This file is part of the repicea library.
 *
 * Copyright (C) 2009-2018 Mathieu Fortin for Rouge-Epicea
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
 * This interface ensures that the plot can provide its sampling weight 
 * @author Mathieu Fortin - December 2018
 */
public interface PlotWeightProvider {

	/**
	 * This method returns the weight of the plot instance. By default, it returns 1.
	 * @return a double
	 */
	default public double getWeight() {
		return 1d;
	}
	
	/**
	 * This method sets the weight of the plot instance.
	 * @param weight a double
	 */
	public void setWeight(double weight);
	
}
