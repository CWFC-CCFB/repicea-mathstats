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
package repicea.simulation.hdrelationships;

/**
 * This interface ensures the REpiceaPredictor instance can predict tree heights
 * @author Mathieu Fortin - June 2019
 */
public interface HeightPredictor<Stand extends HDRelationshipStand, Tree extends HDRelationshipTree> {


	/**
	 * Predicts the height for individual trees and also implements the Monte Carlo simulation automatically. In case of 
	 * exception, it also returns -1. If the predicted height is lower than 1.3, this method returns 1.3.
	 * @param stand a HDRelationshipStand-derived instance
	 * @param tree a HDRelationshipTree-derived instance
	 * @return tree height (m)
	 */
	public double predictHeightM(Stand stand, Tree tree);

}
