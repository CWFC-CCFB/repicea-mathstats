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
package repicea.math.derivative;

import java.security.InvalidParameterException;

import repicea.math.MathematicalFunction;

/**
 * This class implements several numerical differentiation techniques.
 * @author Mathieu Fortin - August 2023
 */
public class NumericalDerivativeCalculator {

	private static double EPS = Math.ulp(1.0);
	private static double SQRT_EPS = Math.sqrt(EPS);
	
	public static enum Technique {
		/**
		 * The Newton's difference quotient.<p>
		 * The derivative is approximated as (f(x+h) - f(x))/h.  
		 */
		NewtonDifferenceQuotient,
		/**
		 * The Symmetric difference quotient.<p>
		 * The derivative is approximated as (f(x+h) - f(x-h))/2h.  
		 */
		SymmetricDifferenceQuotient;
	}
		
	private static void setValueInFunction(MathematicalFunction f, int index, boolean isParameter, double value) {
		if (isParameter) f.setParameterValue(index, value);
		else f.setVariableValue(index, value);
	}
	
	/**
	 * Compute the numerical first derivative of a mathematical function.
	 * @param t a Technique enum 
	 * @param f a MathematicalFunction instance
	 * @param index the index of the parameter or variable to be differentiated
	 * @param isParameter a boolean true if the derivative is calculated with respect to a parameter or false if it is with respect to a variable
	 * @return a double
	 */
	public static double computeNumericalFirstDerivative(Technique t, MathematicalFunction f, int index, boolean isParameter) {
		if (t == null) {
			throw new InvalidParameterException("Parameter t cannot be null!");
		}
		if (f == null) {
			throw new InvalidParameterException("Parameter f cannot be null!");
		}
		double originalValue = isParameter ? f.getParameterValue(index) : f.getVariableValue(index);
		double slope, h, xph, dx;
		switch(t) {
		case NewtonDifferenceQuotient:
			h = originalValue == 0 ? SQRT_EPS : SQRT_EPS * originalValue;
			xph = originalValue + h;
			dx = xph - originalValue;
			double F_x = f.getValue();
			setValueInFunction(f, index, isParameter, xph);
			slope = (f.getValue() - F_x) / dx;
			break;
		case SymmetricDifferenceQuotient:
			h = originalValue == 0 ? SQRT_EPS : SQRT_EPS * originalValue;
//			h = SQRT_EPS / 2;
			xph = originalValue + h;
			double xmh = originalValue - h;
			dx = xph - xmh;
			setValueInFunction(f, index, isParameter, xmh);
			double F_xmh = f.getValue();
			setValueInFunction(f, index, isParameter, xph);
			slope = (f.getValue() - F_xmh) / dx;
			break;
		default:
			throw new InvalidParameterException("This technique is not supported yet: " + t.name());
		}
		setValueInFunction(f, index, isParameter, originalValue);
		return(slope);
	}

//	TODO the square of h at the denominator is too small for the machine precision MF20230816
//	/**
//	 * Compute the numerical second derivative of a mathematical function.<p>
//	 * This differentiation is based on Newton's difference quotient.
//	 * @param f a MathematicalFunction instance
//	 * @param index the index of the parameter or variable to be differentiated
//	 * @param isParameter a boolean true if the derivative is calculated with respect to a parameter or false if it is with respect to a variable
//	 * @return a double
//	 */
//	public static double computeNumericalSecondDerivative(MathematicalFunction f, int index, boolean isParameter) {
//		if (f == null) {
//			throw new InvalidParameterException("Parameter f cannot be null!");
//		}
//		double originalValue = isParameter ? f.getParameterValue(index) : f.getVariableValue(index);
//		double F1_x = f.getValue();
//		double h = originalValue == 0 ? SQRT_EPS : SQRT_EPS * originalValue;
//		double xph = originalValue + h;
//		setValueInFunction(f, index, isParameter, xph);
//		double F1_xph = f.getValue();
//		double xp2h = originalValue + 2 * h;
//		setValueInFunction(f, index, isParameter, xp2h);
//		double F1_xp2h = f.getValue();
//		
//		double dx = (xp2h - originalValue) * .5 ;
//		
////		double F1_x = computeNumericalFirstDerivative(Technique.NewtonDifferenceQuotient, f, index, isParameter);
////		setValueInFunction(f, index, isParameter, xph);
////		double F1_xph = computeNumericalFirstDerivative(Technique.NewtonDifferenceQuotient, f, index, isParameter);
//		double numerator = F1_xp2h - 2 * F1_xph + F1_x;
//		double slope = numerator / (dx * dx);
//		setValueInFunction(f, index, isParameter, originalValue);
//		return(slope);
//	}

	
	
	
	
		
}
