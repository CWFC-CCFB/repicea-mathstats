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

import java.util.ArrayList;
import java.util.List;

import org.junit.Assert;
import org.junit.Test;

/**
 * Tests for the behavior of the ComplexMatrix class.
 * @author Mathieu Fortin - November 2023
 */
public class ComplexMatrixTest {

	@Test
	public void complexMatricesProductTest() {
		List<ComplexNumber> cNumbers = new ArrayList<ComplexNumber>();
		cNumbers.add(new ComplexNumber(1,2));
		cNumbers.add(new ComplexNumber(1,4));
		ComplexMatrix cMatrix = new ComplexMatrix(cNumbers); 
		ComplexMatrix product = cMatrix.multiply(cMatrix.transpose());
		Assert.assertEquals("Testing value at 1,1", new ComplexNumber(-3,4), product.getValueAt(0, 0));
		Assert.assertEquals("Testing value at 1,2", new ComplexNumber(-7,6), product.getValueAt(0, 1));
		Assert.assertEquals("Testing value at 2,1", new ComplexNumber(-7,6), product.getValueAt(1, 0));
		Assert.assertEquals("Testing value at 2,2", new ComplexNumber(-15,8), product.getValueAt(1, 1));
	}
	
	
}
