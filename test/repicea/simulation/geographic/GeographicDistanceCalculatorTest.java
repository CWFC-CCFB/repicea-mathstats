/*
 * This file is part of the repicea library.
 *
 * Copyright (C) 2009-2022 Mathieu Fortin for Rouge-Epicea
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
package repicea.simulation.geographic;

import org.junit.Assert;
import org.junit.Test;

import repicea.math.Matrix;

public class GeographicDistanceCalculatorTest {
	
	@Test
	public void testDistanceLatitudeDifference() {
		double observed = GeographicDistanceCalculator.getCircumferenceKmAtThisLatitudeDeg(60);
		double expected = GeographicDistanceCalculator.EarthCircumferenceKmEquator * .5;
		Assert.assertEquals("Comparing distance in latitudes", expected, observed, 1E-8);
	}
	
	
	@Test
	public void testRatioDegreesKm() {
		double observed = GeographicDistanceCalculator.getRatioLongitudeDegKmAtThisLatitude(60);
		double expected = 360d / GeographicDistanceCalculator.getCircumferenceKmAtThisLatitudeDeg(60);
		Assert.assertEquals("Comparing ratios", expected, observed, 1E-8);
	}
	
	@Test
	public void testDistanceBetweenTheseTwoCoordinates() {
		Matrix latitudes = new Matrix(3,1,49,1);
		Matrix longitudes = new Matrix(3,1,-71,2);
		Matrix distances = GeographicDistanceCalculator.getDistanceBetweenTheseCoordinates(latitudes, longitudes);
		Assert.assertTrue("Comparing diagonal", !distances.diagonalVector().anyElementDifferentFrom(0d));
		Assert.assertEquals("Comparing value at 0, 1", 182.3664179743962, distances.getValueAt(0, 1), 1E-8);
	}
	

	
}
