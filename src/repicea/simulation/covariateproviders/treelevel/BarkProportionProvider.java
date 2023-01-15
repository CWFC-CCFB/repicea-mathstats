/*
 * This file is part of the repicea library.
 *
 * Copyright (C) 2020 Mathieu Fortin  
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
 * This interface ensures the instance can provide the bark proportion with respect to
 * the wood volume.
 * @author Mathieu Fortin - August 2020
 */
public interface BarkProportionProvider {

	/**
	 * Returns the bark proportion with respect to the wood volume.
	 * @return the bark proportion 
	 */
	public double getBarkProportionOfWoodVolume();

	
}
