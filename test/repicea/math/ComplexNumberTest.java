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
import java.util.ArrayList;
import java.util.List;

import org.junit.Assert;
import org.junit.Test;

public class ComplexNumberTest {

	@Test
	public void performSquareAndThenSquareRoot() {
		ComplexNumber c1 = new ComplexNumber(4,1.4);
		ComplexNumber sqr = c1.square();
		ComplexNumber c2 = (ComplexNumber) sqr.sqrt();
		Assert.assertEquals("Comparing real parts", c1.realPart, c2.realPart, 1E-15);
		Assert.assertEquals("Comparing imaginary parts", c1.imaginaryPart, c2.imaginaryPart, 1E-15);
		Assert.assertTrue("Testing equals method", c1.equals(c2));
	}

	@Test
	public void performExpAndLn() {
		ComplexNumber c1 = new ComplexNumber(4,1.4);
		ComplexNumber exp = c1.exp();
		ComplexNumber c2 = exp.log();
		Assert.assertEquals("Comparing real parts", c1.realPart, c2.realPart, 1E-15);
		Assert.assertEquals("Comparing imaginary parts", c1.imaginaryPart, c2.imaginaryPart, 1E-15);
		Assert.assertTrue("Testing equals method", c1.equals(c2));
	}
	
	@Test
	public void comparePow2AndSquare() {
		ComplexNumber c1 = new ComplexNumber(4,1.4);
		ComplexNumber pow2 = c1.pow(2d);
		ComplexNumber sqr = c1.square();
		Assert.assertEquals("Test real part", pow2.realPart, sqr.realPart, 1E-10);
		Assert.assertEquals("Test imag part", pow2.imaginaryPart, sqr.imaginaryPart, 1E-10);
	}
	
	@Test
	public void testNumberInterfaceMethods() {
		ComplexNumber c1 = new ComplexNumber(2.5, 1);
		Assert.assertEquals("Test double", c1.realPart, c1.doubleValue(), 1E-15);
		Assert.assertEquals("Test float", c1.realPart, c1.floatValue(), 1E-8);
		Assert.assertEquals("Test integer", 2, c1.intValue());
		Assert.assertEquals("Test long", 2, c1.longValue());
	}
	
	@Test
	public void constructorWithInvalidParameters1() {
		try {
			new ComplexNumber(4, Double.NaN);
			Assert.fail("Should have thrown an InvalidParameterException instance.");
		} catch (InvalidParameterException e) {}
	}
	
	@Test
	public void performTwiceReciprocal() {
		ComplexNumber c1 = new ComplexNumber(4,1.4);
		ComplexNumber c2 = c1.getReciprocal();
		Assert.assertTrue("Checking c2 is different from c1", !c1.equals(c2));
		ComplexNumber c3 = c2.getReciprocal();
		Assert.assertEquals("Comparing real parts", c1.realPart, c3.realPart, 1E-15);
		Assert.assertEquals("Comparing imaginary parts", c1.imaginaryPart, c3.imaginaryPart, 1E-15);
	}

	@Test
	public void performReciprocalAndThenMultiply() {
		ComplexNumber c1 = new ComplexNumber(4,1.4);
		ComplexNumber c2 = c1.getReciprocal();
		Assert.assertTrue("Checking c2 is different from c1", !c1.equals(c2));
		ComplexNumber c3 = c2.multiply(c1);
		Assert.assertEquals("Comparing real parts", 1d, c3.realPart, 1E-15);
		Assert.assertEquals("Comparing imaginary parts", 0d, c3.imaginaryPart, 1E-15);
	}
	
	@Test
	public void testAbsoluteValue() {
		ComplexNumber c1 = new ComplexNumber(4,1.4);
		ComplexNumber c2 = c1.multiply(c1.getComplexConjugate());
		Assert.assertEquals("Comparing real parts", c1.absoluteValue * c1.absoluteValue, c2.realPart, 1E-15);
		Assert.assertEquals("Comparing imaginary parts", 0d, c2.imaginaryPart, 1E-15);
	}
	
	@Test
	public void testDivideVsReciprocal() {
		ComplexNumber c1 = new ComplexNumber(4,1.4);
		ComplexNumber c2 = new ComplexNumber(5,0.8);
		ComplexNumber c3 = c1.divide(c2);
		ComplexNumber c4 = c1.multiply(c2.getReciprocal());
		Assert.assertEquals("Comparing real parts", c3.realPart, c4.realPart, 1E-15);
		Assert.assertEquals("Comparing imaginary parts", c3.imaginaryPart, c4.imaginaryPart, 1E-15);
	}
	
	@Test
	public void performDivideAndThenMultiply() {
		ComplexNumber c1 = new ComplexNumber(4,1.4);
		ComplexNumber c2 = new ComplexNumber(5,0.8);
		ComplexNumber c3 = c1.divide(c2);
		ComplexNumber c4 = c3.multiply(c2);
		Assert.assertEquals("Comparing real parts", c1.realPart, c4.realPart, 1E-15);
		Assert.assertEquals("Comparing imaginary parts", c1.imaginaryPart, c4.imaginaryPart, 1E-15);
	}

	@Test
	public void testExpectation() {
		List<ComplexNumber> sample = new ArrayList<ComplexNumber>();
		sample.add(new ComplexNumber(2,1));
		sample.add(new ComplexNumber(-1,1.5));
		sample.add(new ComplexNumber(3,7));
		ComplexNumber c = ComplexNumber.getExpectation(sample);
		Assert.assertEquals("Comparing real part", 4d/3, c.realPart, 1E-15);
		Assert.assertEquals("Comparing imaginary part", 9.5/3, c.imaginaryPart, 1E-15);
	}
	
	@Test
	public void testVariance() {
		List<ComplexNumber> sample = new ArrayList<ComplexNumber>();
		sample.add(new ComplexNumber(2,1));
		sample.add(new ComplexNumber(-1,1.5));
		sample.add(new ComplexNumber(3,7));
		double c = ComplexNumber.getVariance(sample);
		Assert.assertEquals("Testing variance", 4.3333333333333333333333 + 11.083333333333333333, c, 1E-15);
	}
	
	
}
