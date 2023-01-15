/*
 * This file is part of the repicea-statistics library.
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
package repicea.math;

import org.junit.Assert;
import org.junit.BeforeClass;
import org.junit.Test;

public class DiagonalMatrixTest {

	static DiagonalMatrix DiagMat;
	
	@BeforeClass
	public static void initialize() {
		DiagonalMatrix sm = new DiagonalMatrix(3);
		sm.setValueAt(0, 0, 1);
		sm.setValueAt(1, 1, 4);
		sm.setValueAt(2, 2, 6);
		DiagMat = sm;
	}
	
	
	@Test
	public void testInternalArrayConstruction() {
		DiagonalMatrix sm = new DiagonalMatrix(4);
		Assert.assertEquals("Testing the array implementation", 1, sm.m_afData.length);
		Assert.assertEquals("Testing the array implementation", 4, sm.m_afData[0].length);
	}

	@Test
	public void testValueAttribution() {
		SymmetricMatrix sm = DiagMat.getDeepClone();
		Assert.assertEquals("Testing getValueAt", 4, sm.getValueAt(1,1), 1E-12);
		Assert.assertEquals("Testing getValueAt", 0, sm.getValueAt(2,1), 1E-12);
		
		sm.setValueAt(2, 2, 15);
		Assert.assertEquals("Testing getValueAt", 15, sm.getValueAt(2,2), 1E-12);
		try {
			sm.setValueAt(0, 1, 0);
			Assert.fail("Should have thrown an exception!");
		} catch(Exception e) {}
//		System.out.println(sm.toString());
	}

	@Test
	public void testMatrixMultiplicationByItself() {
		DiagonalMatrix sm = DiagMat.getDeepClone();
		
		Matrix smPow2 = sm.multiply(sm);
		
		Assert.assertTrue("Testing if the resulting matrix smPow2 is a SymmetricMatrix instance", smPow2 instanceof DiagonalMatrix);
		Assert.assertEquals("Testing a value of the resulting matrix", sm.getValueAt(1, 1) * sm.getValueAt(1, 1), smPow2.getValueAt(1, 1), 1E-8);
	}
	
	
	@Test 
	public void testInverseMatrix() {
		DiagonalMatrix sm = DiagMat.getDeepClone();
		DiagonalMatrix invSM = sm.getInverseMatrix();
		Matrix smTimesInvSM = sm.multiply(invSM);
		
		Assert.assertTrue("Checking if the in-class getInverseMatrix method produces a DiagonalMatrix instance", invSM instanceof DiagonalMatrix);

		Assert.assertTrue("Checking if the product of the inverse by the original matrix returns the identity matrix", smTimesInvSM.equals(Matrix.getIdentityMatrix(sm.m_iRows)));
	}

	
}
