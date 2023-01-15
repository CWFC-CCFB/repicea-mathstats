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
package repicea.stats.data;

import org.junit.Assert;
import org.junit.Test;

public class ObservationTest {

	
	@Test
	public void testEquality1() {
		Object[] values1 = new Object[] {"Allo", 2.345332, 1};
		Object[] values2 = new Object[] {"Allo", 2.345332, 1};
		Observation obs1 = new Observation(values1);
		Observation obs2 = new Observation(values2);
		Assert.assertTrue("Testing that the two observations are equal",
				obs1.isEqualToThisObservation(obs2));
	}

	@Test
	public void testEquality2() {
		Object[] values1 = new Object[] {"Allo", 2.345332, 1};
		Observation obs1 = new Observation(values1);
		Observation obs2 = null;
		Assert.assertTrue("Testing that one observation compared to null returns false",
				!obs1.isEqualToThisObservation(obs2));
	}

	@Test
	public void testEquality3() {
		Object[] values1 = new Object[] {"Allo", 2.345332, 1};
		Object[] values2 = new Object[] {"Allo", 2.345332, 1, "allo"};
		Observation obs1 = new Observation(values1);
		Observation obs2 = new Observation(values2);
		Assert.assertTrue("Testing that the two observations are unequal",
				!obs1.isEqualToThisObservation(obs2));
	}


	@Test
	public void testEqualityForDoubles() {
		Object[] values1 = new Object[] {"Allo", 2.3453320000000000000000, 1};
		Object[] values2 = new Object[] {"Allo", 2.3453320000000000000001, 1};
		Observation obs1 = new Observation(values1);
		Observation obs2 = new Observation(values2);
		Assert.assertTrue("Testing that the two doubles are equal",
				obs1.isEqualToThisObservation(obs2));
	}

	@Test
	public void testInequalityForDoubles() {
		Object[] values1 = new Object[] {"Allo", 2.34533200, 1};
		Object[] values2 = new Object[] {"Allo", 2.34533202, 1};
		Observation obs1 = new Observation(values1);
		Observation obs2 = new Observation(values2);
		Assert.assertTrue("Testing that the two doubles are unequal",
				!obs1.isEqualToThisObservation(obs2));
	}

}
