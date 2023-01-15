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

import java.security.InvalidParameterException;

import repicea.math.Matrix;
import repicea.math.SymmetricMatrix;

/**
 * A class of static methods for the calculation of approximate distances between geographical coordinates.
 * @author Mathieu Fortin - August 2022
 */
public class GeographicDistanceCalculator {

	static final double EarthCircumferenceKmEquator = 40075.0167;
	static final double EarthCircumferenceKmPole = 40007.863;

	/**
	 * Calculate an approximate Earth circumference at a given latitude.
	 * @param latitudeDeg the latitude in degrees
	 * @return the circumference (km)
	 */
	public static double getCircumferenceKmAtThisLatitudeDeg(double latitudeDeg) {
		checkLatitudeDeg(latitudeDeg);
		double latitudeRad = latitudeDeg / 360 * 2 * Math.PI;
		return(Math.cos(latitudeRad) * EarthCircumferenceKmEquator);
	} 

	private static void checkLatitudeDeg(double latitudeDeg) {
		if (latitudeDeg > 90 || latitudeDeg < -90) {
			throw new InvalidParameterException("The latitudeDeg argument must range between -90 and +90!");
		}
	}
	
	private static void checkLongitudeDeg(double longitudeDeg) {
		if (longitudeDeg > 180 || longitudeDeg < -180) {
			throw new InvalidParameterException("The longitudeDeg argument must range between -180 and +180!");
		}
	}

	/**
	 * Return the ratio degrees:km for the the longitude at a particular latitude.
	 * @param latitudeDeg the latitude in degrees
	 * @return the ratio
	 */
	public static double getRatioLongitudeDegKmAtThisLatitude(double latitudeDeg) {
		checkLatitudeDeg(latitudeDeg);
		return 360d / getCircumferenceKmAtThisLatitudeDeg(latitudeDeg);
	}

	/**
	 * Return the distance from the difference of two latitudes.
	 * @param latitudeDeg1 latitude (degrees) of the first point
	 * @param latitudeDeg2 latitude (degrees) of the second point
	 * @return the distance (km)
	 */
	public static double getKmFromLatitudeDifferences(double latitudeDeg1, double latitudeDeg2) {
		checkLatitudeDeg(latitudeDeg1);
		checkLatitudeDeg(latitudeDeg2);
		return Math.abs((latitudeDeg1 - latitudeDeg2) * EarthCircumferenceKmPole / 360);
	}
	
	
	/**
	 * Return the approximate distance between two geographical coordinates. <br>
	 * <br>
	 * The distance is the Euclidean distance. So the approximation is good only for small 
	 * differences in latitude and longitude.
	 * 
	 * @param latitudeDeg1 the latitude (degrees) of the first coordinate
	 * @param longitudeDeg1 the longitude (degrees) of the first coordinate
	 * @param latitudeDeg2 the latitude (degrees) of the second coordinate
	 * @param longitudeDeg2 the longitude (degrees) of the second coordinate
	 * @return the distance (km)
	 */
	public static double getGeographicalDistance(double latitudeDeg1, double longitudeDeg1, double latitudeDeg2, double longitudeDeg2) {
		
		checkLongitudeDeg(longitudeDeg1);
		checkLongitudeDeg(longitudeDeg2);
		
		double diffY = getKmFromLatitudeDifferences(latitudeDeg1, latitudeDeg2);
		double meanLat = (latitudeDeg2 + latitudeDeg1) * .5;
		double longitudeRatio = getRatioLongitudeDegKmAtThisLatitude(meanLat);
		double diffX = (longitudeDeg2 - longitudeDeg1) / longitudeRatio;
		return Math.sqrt(diffY * diffY + diffX * diffX);
	}

	/**
	 * Calculate the approximate distances between a set of coordinates.
	 * @param latitudeDeg the latitudes of the coordinates (a Matrix instance that is a column vector) 
	 * @param longitudeDeg the longitudes of the coordinates (a Matrix instance that is a column vector)
	 * @return a SymmetricMatrix instance
	 */
	public static SymmetricMatrix getDistanceBetweenTheseCoordinates(Matrix latitudeDeg, Matrix longitudeDeg) {
		if (!latitudeDeg.isColumnVector() || !latitudeDeg.isTheSameDimension(longitudeDeg)) {
			throw new InvalidParameterException("The latitudeDeg and longitudeDeg Matrix instances must be column vectors of the same size!");
		}
		SymmetricMatrix outputMatrix = new SymmetricMatrix(latitudeDeg.m_iRows);
		for (int i = 0; i < latitudeDeg.m_iRows; i++) {
			for (int j = i; j < latitudeDeg.m_iRows; j++) {
				if (i == j) {
					outputMatrix.setValueAt(i, j, 0d);
				} else {
					double distanceKm = getGeographicalDistance(latitudeDeg.getValueAt(i, 0),
																longitudeDeg.getValueAt(i, 0),
																latitudeDeg.getValueAt(j, 0),
																longitudeDeg.getValueAt(j, 0));
					outputMatrix.setValueAt(i, j, distanceKm);
				}
			}
		}
		return outputMatrix;
	}
	
}
