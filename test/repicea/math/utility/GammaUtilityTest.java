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

import java.io.IOException;

import org.junit.Assert;
import org.junit.Test;

import repicea.io.javacsv.CSVReader;
import repicea.util.ObjectUtility;

public class GammaUtilityTest {

	@Test
	public void testValues() throws IOException {
		String filename = ObjectUtility.getPackagePath(getClass()) + "gammaTest.csv";
		CSVReader reader = new CSVReader(filename);
		Object[] record;
		int i = 0;
		while ((record = reader.nextRecord()) != null) {
			double x = Double.parseDouble(record[1].toString());
			double expectedValue = Double.parseDouble(record[2].toString());
			double actualValue = GammaUtility.gamma(x);
			Assert.assertEquals("Testing observation" + i, expectedValue, actualValue, 1E-8);
			i++;
		}
		reader.close();
		System.out.println("GammaUtility successfully tested on " + i + " observations");
	}
	
	
	@Test
	public void testInverseGammaFunction() {
		for (double d = 2; d < 15; d+=0.5) {
			double gammaValue = GammaUtility.gamma(d);
			double actual = GammaUtility.inverseGamma(gammaValue);
			System.out.println("Expected = " + d + "; Gamma value = " + gammaValue + "; Actual = " + actual);
			double tolerance = 1E-2;
			if (d == 2) {
				tolerance = 2.5E-2;
			}
			Assert.assertEquals("Testing value d = " + d, d, actual, tolerance);
		}
	}
	
	@Test
	public void testDigammaImplementation() {
		double observed = GammaUtility.digamma(0.5);
		double expected = -1.96351002602142;
		Assert.assertEquals("Testing value of digamma implementation", expected, observed, 1E-12);
		observed = GammaUtility.digamma(5);
		expected = 1.5061176684318;
		Assert.assertEquals("Testing value of digamma implementation", expected, observed, 1E-12);
	}

	@Test
	public void testTrigammaImplementation() {
		double observed = GammaUtility.trigamma(0.5);
		double expected = 4.93480220054468;
		Assert.assertEquals("Testing value of trigamma implementation", expected, observed, 1E-12);
		observed = GammaUtility.trigamma(5);
		expected = 0.221322955737115;
		Assert.assertEquals("Testing value of trigamma implementation", expected, observed, 1E-12);
	}

}
