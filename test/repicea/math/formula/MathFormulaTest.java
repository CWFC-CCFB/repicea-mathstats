/*
 * This file is part of the repicea-statistics library.
 *
 * Copyright (C) 2009-2012 Mathieu Fortin for Rouge-Epicea
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
package repicea.math.formula;

import java.security.InvalidParameterException;
import java.util.LinkedHashMap;

import org.junit.Assert;
import org.junit.Test;


public class MathFormulaTest {

	@Test
	public void simpleTest1() {
		MathFormula mathFormula = new MathFormula("2^(2 + 3)", null, null);
		double result = mathFormula.calculate();
		Assert.assertEquals("Testing result", Math.pow(2, 2 + 3), result, 1E-15);
	}

	@Test
	public void simpleTest2() {
		MathFormula mathFormula = new MathFormula("2 - (2 + 3)", null, null);
		double result = mathFormula.calculate();
		Assert.assertEquals("Testing result", -3, result, 1E-15);
	}

	
	@Test
	public void simpleTest3() {
		MathFormula mathFormula = new MathFormula("2^(2 + (3*.5))", null, null);
		double result = mathFormula.calculate();
		Assert.assertEquals("Testing result", Math.pow(2, 2 + (3*.5)), result, 1E-15);
	}

	
	@Test
	public void testWithExponential1() {
		MathFormula mathFormula = new MathFormula("exp(2 + 3)", 
				null, 
				null);
		System.out.println(mathFormula.toString());
		double result = mathFormula.calculate();
		System.out.println("Result is " + result);
		Assert.assertEquals("Testing result", Math.exp(5), result, 1E-15);
	}

	@Test
	public void testWithExponential2() {
		MathFormula mathFormula = new MathFormula("exp(2 + (3*.5))", 
				null, 
				null);
		System.out.println(mathFormula.toString());
		double result = mathFormula.calculate();
		System.out.println("Result is " + result);
		Assert.assertEquals("Testing result", Math.exp(2 + 3*.5), result, 1E-15);
	}

	
	@Test
	public void testWithExponential3() {
		MathFormula mathFormula = new MathFormula("2 + exp(2 + 3)", 
				null, 
				null);
		System.out.println(mathFormula.toString());
		double result = mathFormula.calculate();
		System.out.println("Result is " + result);
		Assert.assertEquals("Testing result", 2 + Math.exp(5), result, 1E-15);
	}

	@Test
	public void testWithLogarithm1() {
		MathFormula mathFormula = new MathFormula("2 + log(2 + 3)", 
				null, 
				null);
		System.out.println(mathFormula.toString());
		double result = mathFormula.calculate();
		System.out.println("Result is " + result);
		Assert.assertEquals("Testing result", 2 + Math.log(5), result, 1E-15);
	}

	@Test
	public void testWithParameterizedLogarithm1() {
		LinkedHashMap<String, Double> parameters = new LinkedHashMap<String, Double>();
		parameters.put("b1", 2d);
		MathFormula mathFormula = new MathFormula("2 + log(b1 + 3)", 
				parameters, 
				null);
		System.out.println(mathFormula.toString());
		double result = mathFormula.calculate();
		System.out.println("Result is " + result);
		Assert.assertEquals("Testing result", 2 + Math.log(5), result, 1E-15);
		mathFormula.setParameter("b1", 1d);
		result = mathFormula.calculate();
		System.out.println("Result is " + result);
		Assert.assertEquals("Testing result", 2 + Math.log(4), result, 1E-15);
	}

	@Test
	public void testWithParameterizedExponential() {
		LinkedHashMap<String, Double> parameters = new LinkedHashMap<String, Double>();
		parameters.put("b1", 2d);
		MathFormula mathFormula = new MathFormula("exp(b1)", 
				parameters, 
				null);
		System.out.println(mathFormula.toString());
		double result = mathFormula.calculate();
		System.out.println("Result is " + result);
		Assert.assertEquals("Testing result", Math.exp(2), result, 1E-15);
		mathFormula.setParameter("b1", 1d);
		result = mathFormula.calculate();
		System.out.println("Result is " + result);
		Assert.assertEquals("Testing result", Math.exp(1), result, 1E-15);
	}

	@Test
	public void testWithNestedLongNamedOperators() {
		MathFormula mathFormula = new MathFormula("2 + log(exp(2 + 3))", 
				null, 
				null);
		System.out.println(mathFormula.toString());
		double result = mathFormula.calculate();
		System.out.println("Result is " + result);
		Assert.assertEquals("Testing result", 2 + 5, result, 1E-15);
	}

	@Test
	public void testWithUnknownLongNamedOperator() {
		try {
			new MathFormula("2 + offset(2 + 3)", null, null);
			Assert.fail("The MathFormula instance should have thrown an Exception!");
		} catch(InvalidParameterException e) {}
	}

}
