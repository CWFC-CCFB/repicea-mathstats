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
import java.util.List;

/**
 * The ComplexNumber class implements a constructor and the methods 
 * to create complex numbers and perform operations on them.
 * @author Mathieu Fortin - Sept 2023
 */
@SuppressWarnings("serial")
public final class ComplexNumber extends Number {

	final double realPart;
	final double imaginaryPart;
	final double absoluteValue;
	
	/**
	 * Constructor.
	 * @param realPart a double that stands for the real part of the complex number.
	 * @param imaginaryPart a double that stands for the imaginary part of the complex number.
	 */
	public ComplexNumber(double realPart, double imaginaryPart) {
		if (Double.isNaN(imaginaryPart)) {throw new InvalidParameterException("The imaginaryPart argument must be a double!");}
		this.realPart = realPart;
		this.imaginaryPart = imaginaryPart;
		this.absoluteValue = Math.sqrt(realPart * realPart + imaginaryPart * imaginaryPart); 
	}
	
	/**
	 * Add a Number instance to this complex number.
	 * @param d a Number instance
	 * @return a ComplexNumber instance if the imaginary part is different from 0. Otherwise, a double.
	 */
	public ComplexNumber add(Number d) {
		double iPart = d instanceof ComplexNumber ? 
				imaginaryPart + ((ComplexNumber) d).imaginaryPart :
					imaginaryPart;
		return new ComplexNumber(realPart + d.doubleValue(), iPart);
	}

	/**
	 * Subtract a Number instance to this complex number.
	 * @param d a Number instance
	 * @return a ComplexNumber instance
	 */
	public ComplexNumber subtract(Number d) {
		double iPart = d instanceof ComplexNumber ? 
				imaginaryPart - ((ComplexNumber) d).imaginaryPart :
					imaginaryPart;
		return new ComplexNumber(realPart - d.doubleValue(), iPart);
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
	public ComplexNumber log() {
		double r = absoluteValue;
		double phi = Math.atan(imaginaryPart / realPart);
		return new ComplexNumber(Math.log(r), phi);
	}
	
	
	/**
	 * Return the absolute value of this complex number.<p>
	 * That is the the square root of Re(z)^2 + Im(z)^2. The absolute value
	 * is also known as the modulus or magnitude.
	 * @return a double
	 */
	public double getAbsoluteValue() {
		return absoluteValue;
	}

	/**
	 * Return the reciprocal of this complex number, that is 1/z.
	 * @return a ComplexNumber instance
	 */
	public ComplexNumber getReciprocal() {
		double r2 = absoluteValue * absoluteValue;
		return new ComplexNumber(realPart/r2, -imaginaryPart/r2);
	}
	
	
	/**
	 * Return the division of this complex number by n.
	 * @param n a Number instance (can be a ComplexNumber instance)
	 * @return a ComplexNumber instance
	 */
	public ComplexNumber divide(Number n) {
		if (n instanceof ComplexNumber) {
			ComplexNumber cN = (ComplexNumber) n;
			double r2 = cN.absoluteValue * cN.absoluteValue;
			double real = realPart * cN.realPart + imaginaryPart * cN.imaginaryPart;
			double imag = imaginaryPart * cN.realPart - realPart * cN.imaginaryPart;
			return new ComplexNumber(real / r2, imag / r2);
		} else {
			return new ComplexNumber(realPart / n.doubleValue(), imaginaryPart / n.doubleValue());
		}
	}
	

	/**
	 * Return the complex conjugate of this complex number. <p>
	 * The complex conjugate of z is Re(z) - Im(z).
	 * @return a ComplexNumber instance
	 */
	public ComplexNumber getComplexConjugate() {
		return new ComplexNumber(realPart, -imaginaryPart);
	}
	
	/**
	 * Compute the square root of this complex number.<p>
	 * It may occasionally return a double if the imaginary part is equal to 0.
	 * @return a ComplexNumber instance
	 */
	public ComplexNumber sqrt() {
		double modulus = Math.sqrt(realPart * realPart + imaginaryPart * imaginaryPart);
		double gamma = Math.sqrt((realPart + modulus) * .5);
		double delta = imaginaryPart > 0 ?
				Math.sqrt((-realPart + modulus) * .5) :
					-Math.sqrt((-realPart + modulus) * .5);
		return new ComplexNumber(gamma, delta);
	}
	
	
	/**
	 * Return the real part of this complex number as a double.
	 * @return a double
	 */
	@Override
	public double doubleValue() {
		return realPart;
	}

	/**
	 * Return the real part of this complex number as a float.<p>
	 * This method implies a truncation.
	 * @return a float
	 */
	@Override
	public float floatValue() {
		return (float) realPart;
	}

	/**
	 * Return the real part of this complex number as an integer.<p>
	 * This method implies a truncation.
	 * @return an integer
	 */
	@Override
	public int intValue() {
		return (int) realPart;
	}

	/**
	 * Return the real part of this complex number as a long.<p>
	 * This method implies a truncation.
	 * @return a long
	 */
	@Override
	public long longValue() {
		return (long) realPart;
	}

	@Override
	public String toString() {
		return "" + realPart + " + " + imaginaryPart + "i"; 
	}
	
	@Override
	public boolean equals(Object c) {
		if (c instanceof ComplexNumber) {
			ComplexNumber cn = (ComplexNumber) c;
			return realPart == cn.realPart && imaginaryPart == cn.imaginaryPart; 
		}
		return false;
	}

	/**
	 * Estimate the expectation of a sample of realizations of a complex random variable.
	 * @param sample a List of ComplexNumber instances
	 * @return a ComplexNumber instance
	 */
	public static ComplexNumber getExpectation(List<ComplexNumber> sample) {
		if (sample == null || sample.isEmpty()) {
			throw new InvalidParameterException("The sample argument must be a non empty List instance!");
		}
		double real = 0;
		double imag = 0;
		for (ComplexNumber c : sample) {
			real += c.realPart;
			imag += c.imaginaryPart;
		}
		return new ComplexNumber(real / sample.size(), imag / sample.size());
	}

	/**
	 * Estimate the variance of a sample of realizations of a complex random variable.<p>
	 * The variance is actually the sum of the variance of the real part and the variance
	 * of the imaginary part.
	 * @param sample a List of ComplexNumber instances
	 * @return a double
	 */
	public static double getVariance(List<ComplexNumber> sample) {
		if (sample == null || sample.size() <= 1) {
			throw new InvalidParameterException("The sample argument must be a List instance of at least two elements!");
		}
		ComplexNumber mean = getExpectation(sample);
		double diff2Real = 0;
		double diff2Imag = 0;
		for (ComplexNumber c : sample) {
			double dReal = c.realPart - mean.realPart;
			diff2Real += dReal * dReal;
			double dImag = c.imaginaryPart - mean.imaginaryPart;
			diff2Imag += dImag * dImag;
		}
		double df = sample.size() - 1;
		return diff2Real / df + diff2Imag / df;
	}
	
}
