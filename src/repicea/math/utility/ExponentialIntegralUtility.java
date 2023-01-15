/*
 * This file is part of the repicea library.
 *
 * Copyright (C) 2022 Mathieu Fortin for Rouge-Epicea
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
package repicea.math.utility;

import java.security.InvalidParameterException;

/**
 * A class that implements the Exponential Integral function.
 * @author Mathieu Fortin - June 2022
 */
public class ExponentialIntegralUtility {

	private final static double EULER_CONSTANT = 0.57721566490153286060651209008240243104215933593992;
	private static int K_MAX = 30;
	
	/**
	 * Returns the result of the E1(z) function using the approximation of Abramowitz and Stegun, p. 229, 5.1.11. 
	 * This approximation is good if -2.5 < z < 2.5. 
	 * @param z a double different from 0
	 * @return a double
	 */
	public static double getE1(double z, Double logTransformed) {
		checkValue(z);
		if (z < 0) {
			return -getAbramowitzAndStegunApproximationForEi(-z, logTransformed);
		}
		if (z == 0d && (logTransformed == null || Double.isInfinite(logTransformed))) {
			return Double.POSITIVE_INFINITY;
		}
		return getAbramowitzAndStegunApproximationForE1(z, logTransformed);
	}
	
	/**
	 * Returns the result of the E1(z) function using the approximation of Abramowitz and Stegun, p. 229, 5.1.11. 
	 * This approximation is good if -2.5 < z < 2.5. 
	 * @param z a double different from 0
	 * @return a double
	 */
	public static double getE1(double z) {
		return getE1(z, null);
	}


	private static void checkValue(double z) {
		if (z < -2.5 || z > 2.5) {
			System.out.println("Warning: the Abramowitz And Stegun Approximation for the exponential integral is valid only for -2.5 < z 2.5 !");
		}
	}
	
	/**
	 * Returns the result of the Ei(z) function using the approximation of Abramowitz and Stegun, p. 229, 5.1.11. 
	 * @param z a double different from 0
	 * @return a double
	 */
	public static double getEi(double z) {
		return getEi(z, null);
	}


	/**
	 * Returns the result of the Ei(z) function using the approximation of Abramowitz and Stegun, p. 229, 5.1.11. 
	 * @param z a double different from 0
	 * @return a double
	 */
	public static double getEi(double z, Double logTransformed) {
		checkValue(z);
		if (z < 0) {
			return -getAbramowitzAndStegunApproximationForE1(-z, logTransformed);
		} 
		if (z == 0d && (logTransformed == null || Double.isInfinite(logTransformed))) {
			return Double.NEGATIVE_INFINITY;
		}
		return getAbramowitzAndStegunApproximationForEi(z, logTransformed);
	}

	private static double getAbramowitzAndStegunApproximationForE1(double z, Double logTransformed) {
		if (z < 0) {
			throw new InvalidParameterException("The argument z must be equal to or greater than zero!");
		}
		if (z == 0d && logTransformed != null) {
			return -EULER_CONSTANT - logTransformed;
		} else {
			double sum = 0;
			double powMinusZ = 1;
			long factK = 1;
			for (int k = 1; k < K_MAX; k++) {
				powMinusZ *= -z;
				factK *= k;
				sum += powMinusZ / (k * factK);
			}
			return -EULER_CONSTANT - Math.log(z) - sum;
		}
	}

	
	
	private static double getAbramowitzAndStegunApproximationForEi(double z, Double logTransformed) {
		if (z < 0) {
			throw new InvalidParameterException("The argument z must be equal to or greater than zero!");
		}
		if (z == 0d && logTransformed != null) {
			return EULER_CONSTANT + logTransformed;
		} else {
			double sum = 0;
			double powMinusZ = 1;
			long factK = 1;
			for (int k = 1; k < K_MAX; k++) {
				powMinusZ *= z;
				factK *= k;
				sum += powMinusZ / (k * factK);
			}
			return EULER_CONSTANT + Math.log(z) + sum;
		}
	}
	
	

//	public static void main(String[] args) {
//		long startTime = System.currentTimeMillis();
////		for (int i = 0; i < 10000; i++)
//		System.out.println("E1(z) = " + ExponentialIntegralFunction.getE1(2));
//		System.out.println("Ei(z) = " + ExponentialIntegralFunction.getEi(2));
//		System.out.println("Time = " + (System.currentTimeMillis() - startTime) + " ms.");
//	}

}
