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
package repicea.math.utility;

import java.security.InvalidParameterException;

import repicea.math.Matrix;

/**
 * This class implements simple static methods in mathematics.
 * @author Mathieu Fortin - February 2019
 */
public class MathUtility {


	/**
	 * Compute the Euclidean distance between two points. <br>
	 * <br>
	 * This method assumes that the checks have been performed on the coordinates argument. Basically,
	 * these matrices should be column vectors of the same size. Each one of them represents a dimensions.
	 * 
	 * @param i the index of the first point
	 * @param j the index of the second point
	 * @param coordinates A series of column matrices that stand for the coordinates. 
	 * @return
	 */
	public final static double getEuclideanDistance(int i, int j, Matrix... coordinates) {
		double squareDiffSum = 0d;
		for (int k = 0; k < coordinates.length; k++) {
			Matrix c = coordinates[k];
			double diff = c.getValueAt(i, 0) - c.getValueAt(j, 0); 
			squareDiffSum += diff * diff;
		}
		return Math.sqrt(squareDiffSum);
	}

	
	/**
	 * This method returns the factorial of parameter i.
	 * @param i an integer
	 * @return the result as a double
	 */
	public static double Factorial(int i) {
		if (i < 0) {
			throw new InvalidParameterException("Parameter i must be equal to or greater than 0!");
		} else  if (i==0) {
			return 1;
		} else {
			double result = 1;
			for (int j = 1; j <= i; j++) {
				result *= j;
			}
			return result;
		}
	}
	
	
	/**
	 * This method returns the ratio of two factorial factorial of parameter i.
	 * @param i an integer
	 * @return the result as a double
	 */
	public static double FactorialRatio(int i, int j) {
		if (i < 0 || j < 0) {
			throw new InvalidParameterException("Parameters i and j must be equal to or greater than 0!");
		} else if (j > i) {
			throw new InvalidParameterException("Parameter j must be smaller than parameter i!");
		} else {
			if (j == 0) {
				j = 1;
			}
			double result = 1;
			for (int k = i; k > j; k--) {
				result *= k;
			}
			return result;
		}
	}
	
	/**
	 * Makes it possible to determine if a double is negative or not. Useful to identify -0.0 for example.
	 * After Peter Lawrey on <a href=https://stackoverflow.com/questions/10399801/how-to-check-if-double-value-is-negative-or-not> StackOverFlow </a>
	 * @param d
	 * @return
	 */
	public static boolean isNegative(double d) {
	     return Double.doubleToRawLongBits(d) < 0;
	}

}
