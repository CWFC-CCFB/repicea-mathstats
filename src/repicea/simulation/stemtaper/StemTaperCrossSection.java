/*
 * This file is part of the repicea-foresttools library.
 *
 * Copyright (C) 2009-2014 Mathieu Fortin for Rouge-Epicea
 *
 * This library is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public
 * License as published by the Free Software Foundation; either
 * version 2.1 of the License, or (at your option) any later version.
 *
 * This library is distributed with the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied
 * warranty of MERCHANTABILITY or FITNESS FOR A
 * PARTICULAR PURPOSE. See the GNU Lesser General Public
 * License for more details.
 *
 * Please see the license at http://www.gnu.org/copyleft/lesser.html.
 */
package repicea.simulation.stemtaper;

/**
 * This class simply describe a cross section in tree. It includes a height and a diameter. It also implements the Comparable
 * interface in order to sort the section within a given tree.
 * @author Mathieu Fortin - July 2011
 */
@SuppressWarnings("rawtypes")
public class StemTaperCrossSection implements Comparable {

	private double height;
	private double diameter;

	/**
	 * Constructor.
	 * @param height the height of the section along the bole (m)
	 * @param diameter the diameter of the section (mm)
	 */
	public StemTaperCrossSection(double height, double diameter) {
		this.height = height;
		this.diameter = diameter;
	}

	@Override
	public int compareTo(Object arg0) {
		StemTaperCrossSection otherSection = (StemTaperCrossSection) arg0;
		if (height < otherSection.height) {
			return -1;
		} else if (height == otherSection.height) {
			return 0;
		} else {
			return 1;
		}
	}

	/**
	 * This method returns the height of the section (m).
	 * @return a double
	 */
	public double getSectionHeight() {return height;}
	
	/**
	 * This method returns the diameter of the section (mm).
	 * @return a double
	 */
	public double getSectionDiameter() {return diameter;}
	
}
