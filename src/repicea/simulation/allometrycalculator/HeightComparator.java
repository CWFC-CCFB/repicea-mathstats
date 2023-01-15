/*
 * This file is part of the repicea library.
 *
 * Copyright (C) 2009-2015 Mathieu Fortin for Rouge-Epicea
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
package repicea.simulation.allometrycalculator;

import java.util.Comparator;

import repicea.simulation.covariateproviders.treelevel.HeightMProvider;

/**
 * The DbhComparator class implements the Comparator interface for HeightMProvider instances.
 * @author Mathieu Fortin - June 2015
 */
public class HeightComparator implements Comparator<HeightMProvider> {

	private final boolean ascending;
	
	/**
	 * General constructor.
	 * @param ascending true for ascending sorting or false for descending sorting
	 */
	public HeightComparator(boolean ascending) {
		this.ascending = ascending;
	}
	
	/**
	 * Default constructor for ascending comparison.
	 */
	public HeightComparator() {
		this(true);
	}
	
	@Override
	public int compare(HeightMProvider arg0, HeightMProvider arg1) {
		if (ascending) {
			if (arg0.getHeightM() < arg1.getHeightM()) {
				return -1;
			} else if (arg0.getHeightM() == arg1.getHeightM()) {
				return 0;
			} else {
				return 1;
			}
		} else {
			if (arg0.getHeightM() < arg1.getHeightM()) {
				return 1;
			} else if (arg0.getHeightM() == arg1.getHeightM()) {
				return 0;
			} else {
				return -1;
			}
		}
	}

}
