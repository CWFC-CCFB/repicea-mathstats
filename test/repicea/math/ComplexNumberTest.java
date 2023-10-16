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
	
}
