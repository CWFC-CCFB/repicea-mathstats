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
 * This interface ensures the plot instance can tell whether or not it is going
 * to be sprayed against defoliator in the upcoming interval.
 * @author Mathieu Fortin - July 2019
 */
public interface SprayedAgainstDefoliatorProvider {

	/**
	 * This method returns a boolean if the stand is going to be sprayed in
	 * the upcoming interval.
	 * @return a boolean
	 */
	public boolean isSprayed();

}
