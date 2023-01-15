/*
 * This file is part of the repicea library.
 *
 * Copyright (C) 2009-2022 Mathieu Fortin for Rouge-Epicea
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
package repicea.stats.distributions.utility;

import org.junit.Assert;
import org.junit.Test;

public class WeibullUtilityTest {

	@Test
	public void simpleWeibullDensityTest() {
		double d = WeibullUtility.getProbabilityDensity(1, 1, 1);
		Assert.assertEquals("Testing Weibull density", 0.36787944117144233, d, 1E-8);
	}

	@Test
	public void simpleWeibullDensity2Test() {
		double d = WeibullUtility.getProbabilityDensity(1, 2, 1);
		Assert.assertEquals("Testing Weibull density", 0.7357588823428847, d, 1E-8);
	}

	
	@Test
	public void simpleWeibullCDFTest() {
		double d = WeibullUtility.getCumulativeProbability(1, 1, 1);
		Assert.assertEquals("Testing Weibull density", 0.6321205588285577, d, 1E-8);
	}

	@Test
	public void simpleWeibullCDF2Test() {
		double d = WeibullUtility.getCumulativeProbability(1, 1, 2);
		Assert.assertEquals("Testing Weibull density", 0.3934693402873666, d, 1E-8);
	}

}
