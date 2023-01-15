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

import repicea.math.utility.ExponentialIntegralUtility;

public class ExponentialIntegralUtilityTest {

	@Test
	public void testE1() {
		double observed = ExponentialIntegralUtility.getE1(2);
		Assert.assertEquals("Testing the E1(z) function", 0.04890051070808066, observed, 1E-8);
	}
	
	@Test
	public void testEi() {
		double observed = ExponentialIntegralUtility.getEi(2);
		Assert.assertEquals("Testing the Ei(z) function", 4.954234356001867, observed, 1E-8);
	}

	@Test
	public void testE1_negative() {
		double observed = ExponentialIntegralUtility.getE1(-2);
		Assert.assertEquals("Testing the E1(z) function", -4.954234356001867, observed, 1E-8);
	}
	
	@Test
	public void testEi_negative() {
		double observed = ExponentialIntegralUtility.getEi(-2);
		Assert.assertEquals("Testing the Ei(z) function", -0.04890051070808066, observed, 1E-8);
	}

}
