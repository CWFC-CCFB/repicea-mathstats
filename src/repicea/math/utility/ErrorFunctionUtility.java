/*
 * This file is part of the repicea-mathstats library.
 *
 * Copyright (C) 2023 His Majesty the King in Right of Canada
 * Author: Mathieu Fortin, Canadian Wood Fibre Centre, CFS
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

/**
 * This class implements static methods to approximate the value of the error function (erf).
 * @author Mathieu Fortin - August 2023
 */
public class ErrorFunctionUtility {

	
	private static final double p = 0.3275911;
	private static final double a1 = 0.254829592;
	private static final double a2 = -0.284496736;
	private static final double a3 = 1.421413741;
	private static final double a4 = -1.453152027;
	private static final double a5 = 1.061405429;
	

	/**
	 * Calculate an approximation of the error function.<p>
	 * 
	 * This approximation relies on one of Abramowitz and Stegun's approximations.
	 * 
	 * @param x a double
	 * @return a double bounded to [-1,1]
	 */
	public static double erf(double x) {
		boolean neg = MathUtility.isNegative(x);
		double y = neg ? -x : x;
		if (y > 5) {
			return neg ? -1 : 1;
		} else {
			double t = 1 /(1 + p*y);
			double t2 = t * t;
			double approx = 1 - 
					(a1 * t + 
					 a2 * t2 + 
					 a3 * t2 * t +
					 a4 * t2 * t2 +
					 a5 * t2 * t2 * t) * Math.exp(-y*y);
			return neg ? - approx : approx;
		}
	}
	
}
