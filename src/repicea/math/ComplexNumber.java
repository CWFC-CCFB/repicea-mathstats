/*
 * This file is part of the repicea-mathstats library.
 *
 * Copyright (C) 2023 His Majesty the King in Right of Canada
 * Author: Mathieu Fortin, Canadian Wood Fibre Centre, Canadian Forest Service
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
package repicea.math;

import java.security.InvalidParameterException;

/**
 * The ComplexNumber class implements a constructor and the methods 
 * to create complex numbers and perform operations on them.
 * @author Mathieu Fortin - Sept 2023
 */
@SuppressWarnings("serial")
public final class ComplexNumber extends Number {

	final double realPart;
	final double imaginaryPart;
	
	/**
	 * Constructor.
	 * @param realPart a double that stands for the real part of the complex number.
	 * @param imaginaryPart a double that stands for the imaginary part of the complex number.
	 */
	public ComplexNumber(double realPart, double imaginaryPart) {
		if (Double.isNaN(imaginaryPart)) {throw new InvalidParameterException("The imaginaryPart argument must be non null double!");}
		if (imaginaryPart == 0) {throw new InvalidParameterException("The imaginaryPart argument must be non null double!");}
		this.realPart = realPart;
		this.imaginaryPart = imaginaryPart;
	}
	
	/**
	 * Add a Number instance to this complex number.
	 * @param d a Number instance
	 * @return a ComplexNumber instance if the imaginary part is different from 0. Otherwise, a double.
	 */
	public Number add(Number d) {
		double iPart = d instanceof ComplexNumber ? 
				imaginaryPart + ((ComplexNumber) d).imaginaryPart :
					imaginaryPart;
		return iPart != 0 ? 
				new ComplexNumber(realPart + d.doubleValue(), iPart) : 
					realPart + d.doubleValue();
	}

	/**
	 * Subtract a Number instance to this complex number.
	 * @param d a Number instance
	 * @return a ComplexNumber instance if the imaginary part is different from 0. Otherwise, a double.
	 */
	public Number subtract(Number d) {
		double iPart = d instanceof ComplexNumber ? 
				imaginaryPart - ((ComplexNumber) d).imaginaryPart :
					imaginaryPart;
		return iPart != 0 ? 
				new ComplexNumber(realPart - d.doubleValue(), iPart) : 
					realPart - d.doubleValue();
	}

	/**
	 * Compute the product of this random number by a Number instance. <p>
	 * @param d a Number instance which can be a ComplexNumber instance.
	 * @return a ComplexNumber instance
	 */
	public ComplexNumber multiply(Number d) {
		if (d instanceof ComplexNumber) {
			double rPart = this.realPart * ((ComplexNumber) d).realPart - this.imaginaryPart * ((ComplexNumber) d).imaginaryPart;
			double iPart = this.realPart * ((ComplexNumber) d).imaginaryPart + this.imaginaryPart * ((ComplexNumber) d).realPart;
			return new ComplexNumber(rPart, iPart);
		} else {
			return new ComplexNumber(this.realPart * d.doubleValue(), this.imaginaryPart * d.doubleValue());
		}
	}
	
	/**
	 * Compute the square of this complex number.
	 * @return a ComplexNumber instance
	 */
	public ComplexNumber square() {
		return (ComplexNumber) multiply(this);
	}
	
	/**
	 * Compute the exponential of this complex number.<p>
	 * The method is based on Euler's formula.
	 * @return a ComplexNumber instance
	 */
	public ComplexNumber exp() {
		double multiplier = Math.exp(realPart);
		return new ComplexNumber(multiplier * Math.cos(imaginaryPart), multiplier * Math.sin(imaginaryPart));
	}
	

	/**
	 * Compute the natural logarithm of this complex number.<p>
	 * The method relies on the polar form to make the calculation.
	 * @return a ComplexNumber instance
	 */
	public ComplexNumber ln() {
		double r = Math.sqrt(realPart * realPart + imaginaryPart * imaginaryPart);
		double phi = Math.atan(imaginaryPart / realPart);
		return new ComplexNumber(Math.log(r), phi);
	}
	

	/**
	 * Compute the square root of this complex number.<p>
	 * It may occasionally return a double if the imaginary part is equal to 0.
	 * @return a Number instance
	 */
	public Number sqrt() {
		double modulus = Math.sqrt(realPart * realPart + imaginaryPart * imaginaryPart);
		double gamma = Math.sqrt((realPart + modulus) * .5);
		double delta = imaginaryPart > 0 ?
				Math.sqrt((-realPart + modulus) * .5) :
					-Math.sqrt((-realPart + modulus) * .5);
		return delta != 0 ? 
				new ComplexNumber(gamma, delta) :
					gamma;
	}
	
	
	/**
	 * Returns the real part of this complex number as a double.
	 * @return a double
	 */
	@Override
	public double doubleValue() {
		return realPart;
	}

	/**
	 * Returns the real part of this complex number as a float.<p>
	 * This method implies a truncation.
	 * @return a float
	 */
	@Override
	public float floatValue() {
		return (float) realPart;
	}

	/**
	 * Returns the real part of this complex number as an integer.<p>
	 * This method implies a truncation.
	 * @return an integer
	 */
	@Override
	public int intValue() {
		return (int) realPart;
	}

	/**
	 * Returns the real part of this complex number as a long.<p>
	 * This method implies a truncation.
	 * @return a long
	 */
	@Override
	public long longValue() {
		return (long) realPart;
	}

	@Override
	public String toString() {
		return "" + realPart + " + i" + imaginaryPart; 
	}
	
	@Override
	public boolean equals(Object c) {
		if (c instanceof ComplexNumber) {
			ComplexNumber cn = (ComplexNumber) c;
			return realPart == cn.realPart && imaginaryPart == cn.imaginaryPart; 
		}
		return false;
	}
	
}
