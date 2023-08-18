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

import org.junit.Assert;
import org.junit.Test;

public class ErrorFunctionUtilityTest {
	
	@Test
	public void errorFunctionTest1() {
		double actual = ErrorFunctionUtility.erf(0);
		Assert.assertEquals("Testing approximation for x = 0", 0, actual, 1.5E-7);
	}
		
	@Test
	public void errorFunctionTest2() {
		double actual = ErrorFunctionUtility.erf(0.5);
		Assert.assertEquals("Testing approximation for x = 0.5", 0.5204999, actual, 1.5E-7);
	}
	
	@Test
	public void errorFunctionTest3() {
		double actual = ErrorFunctionUtility.erf(-0.5);
		Assert.assertEquals("Testing approximation for x = -0.5", -0.5204999, actual, 1.5E-7);
	}

	@Test
	public void errorFunctionTest4() {
		double actual = ErrorFunctionUtility.erf(2.5);
		Assert.assertEquals("Testing approximation for x = 2.5", 0.999593, actual, 1.5E-7);
	}
	
	@Test
	public void errorFunctionTest5() {
		double actual = ErrorFunctionUtility.erf(-2.5);
		Assert.assertEquals("Testing approximation for x = -2.5", -0.999593, actual, 1.5E-7);
	}

	@Test
	public void errorFunctionTest6() {
		double actual = ErrorFunctionUtility.erf(5);
		Assert.assertEquals("Testing approximation for x = 5", 1, actual, 1.5E-7);
	}

	
	@Test
	public void errorFunctionTest7() {
		double actual = ErrorFunctionUtility.erf(-5);
		Assert.assertEquals("Testing approximation for x = -5", -1, actual, 1.5E-7);
	}

	@Test
	public void errorFunctionTest8() {
		double actual = ErrorFunctionUtility.erf(10);
		Assert.assertEquals("Testing approximation for x = 10", 1, actual, 1.5E-7);
	}
	
	@Test
	public void errorFunctionTest9() {
		double actual = ErrorFunctionUtility.erf(-10);
		Assert.assertEquals("Testing approximation for x = -10", -1, actual, 1.5E-7);
	}

}
