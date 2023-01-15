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
package repicea.stats.model.glm;

import java.security.InvalidParameterException;

import org.junit.Assert;
import org.junit.Test;

import repicea.stats.model.glm.Family.GLMDistribution;

public class FamilyTest {

	@Test
	public void testFamilyException() {
		try {
			Family.createFamily(GLMDistribution.Bernoulli, LinkFunction.Type.Log, null);
			Assert.fail("Should have thrown an exception because the type is incompatible with the family!");
		} catch(InvalidParameterException e) {}
		try {
			Family.createFamily(GLMDistribution.NegativeBinomial, LinkFunction.Type.Logit, null);
			Assert.fail("Should have thrown an exception because the type is incompatible with the family!");
		} catch(InvalidParameterException e) {}
	}
}
