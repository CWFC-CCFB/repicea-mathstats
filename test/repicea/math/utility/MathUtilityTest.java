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

import org.junit.Assert;
import org.junit.Test;


public class MathUtilityTest {

	@Test
	public void testSimpleFactorials() {
		double actual = MathUtility.Factorial(5);
		Assert.assertEquals("Testing factorial 5", 120, actual, 1E-8);

		actual = MathUtility.Factorial(0);
		Assert.assertEquals("Testing factorial 0", 1, actual, 1E-8); 
		
		actual = MathUtility.Factorial(10);
		Assert.assertEquals("Testing factorial 10", 3628800, actual, 1E-8);
	}
	
	
	@Test
	public void testFactorialRatios() {
		double actual = MathUtility.FactorialRatio(5, 3);
		Assert.assertEquals("Testing 5!/3!", 20, actual, 1E-8);
		
		actual = MathUtility.FactorialRatio(10, 7);
		Assert.assertEquals("Testing 10!/7!", 720, actual, 1E-8);

		actual = MathUtility.FactorialRatio(1, 0);
		Assert.assertEquals("Testing 1!/0!", 1, actual, 1E-8);
	}

}
