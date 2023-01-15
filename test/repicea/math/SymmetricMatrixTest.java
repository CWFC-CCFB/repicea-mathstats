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
package repicea.math;

import org.junit.Assert;
import org.junit.BeforeClass;
import org.junit.Test;


public class SymmetricMatrixTest {

	static SymmetricMatrix SymMat;
	
	@BeforeClass
	public static void initialize() {
		SymmetricMatrix sm = new SymmetricMatrix(3);
		sm.setValueAt(0, 0, 1);
		sm.setValueAt(0, 1, 2);
		sm.setValueAt(0, 2, 3);
		sm.setValueAt(1, 1, 4);
		sm.setValueAt(1, 2, 5);
		sm.setValueAt(2, 2, 6);
		SymMat = sm;
	}
	
	
	@Test
	public void testInternalArrayConstruction() {
		SymmetricMatrix sm = new SymmetricMatrix(4);
		Assert.assertEquals("Testing the array implementation", 4, sm.m_afData[0].length);
		Assert.assertEquals("Testing the array implementation", 1, sm.m_afData[3].length);
	}

	@Test
	public void testValueAttribution() {
		SymmetricMatrix sm = SymMat.getDeepClone();
		Assert.assertEquals("Testing getValueAt", 2, sm.getValueAt(1,0), 1E-12);
		Assert.assertEquals("Testing getValueAt", 5, sm.getValueAt(2,1), 1E-12);
		
		sm.setValueAt(2, 1, 15);
		Assert.assertEquals("Testing getValueAt", 15, sm.getValueAt(1,2), 1E-12);
//		System.out.println(sm.toString());
	}

	@Test
	public void testMatrixMultiplicationByItself() {
		SymmetricMatrix sm = SymMat.getDeepClone();
		
//		System.out.println("Matrix sm = " + sm.toString());
		Matrix smPow2 = sm.multiply(sm);
		
		Matrix m = new Matrix(3,3);
		m.setValueAt(0, 0, 1);
		m.setValueAt(0, 1, 2);
		m.setValueAt(0, 2, 3);
		m.setValueAt(1, 0, 2);
		m.setValueAt(1, 1, 4);
		m.setValueAt(1, 2, 5);
		m.setValueAt(2, 0, 3);
		m.setValueAt(2, 1, 5);
		m.setValueAt(2, 2, 6);
//		System.out.println("Matrix m = " + m.toString());
	
		Matrix mPow2 = m.multiply(m);
		Assert.assertTrue("Testing if the resulting matrix smPow2 is a SymmetricMatrix instance", smPow2 instanceof SymmetricMatrix);
		Assert.assertTrue("Testing if the resulting matrix mPow2 is a not SymmetricMatrix instance", !(mPow2 instanceof SymmetricMatrix));
		Assert.assertTrue("Testing that the two matrices are equal", smPow2.equals(mPow2));
	}
	
	
	@Test 
	public void testKroneckerProductOfSymmetricMatrices() {
		SymmetricMatrix sm = new SymmetricMatrix(2);
		sm.setValueAt(0, 0, 10);
		sm.setValueAt(0, 1, 20);
		sm.setValueAt(1, 1, 40);

		Matrix m = new Matrix(2,2);
		m.setValueAt(0, 0, 10);
		m.setValueAt(0, 1, 20);
		m.setValueAt(1, 0, 20);
		m.setValueAt(1, 1, 40);

		Matrix kProd = SymMat.getKroneckerProduct(m);
		Assert.assertTrue("Checking if the regular Kronecker product is symmetric", kProd.isSymmetric());
		
		Matrix kProd2 = SymMat.getKroneckerProduct(sm);
		Assert.assertTrue("Checking if the in-class Kronecker product implementation produces a SymmetricMatrix instance", kProd2 instanceof SymmetricMatrix);

		Assert.assertTrue("Checking if the two resulting matrices are equal", kProd2.equals(kProd));
	}
	
	
	@Test 
	public void testInverseMatrix() {
		Matrix sm = SymMat.getDeepClone();
		Matrix invSM = sm.getInverseMatrix();
		Matrix smTimesInvSM = sm.multiply(invSM);
		
		Assert.assertTrue("Checking if the in-class getInverseMatrix method produces a SymmetricMatrix instance", invSM instanceof SymmetricMatrix);

		Assert.assertTrue("Checking if the product of the inverse by the original matrix returns the identity matrix", smTimesInvSM.equals(Matrix.getIdentityMatrix(sm.m_iRows)));
	}

	@Test 
	public void testInverseMatrixWithMultiplyBlocks() {
		Matrix sm = SymMat.matrixDiagBlock(SymMat);
		Assert.assertTrue("Checking if the sm object is a SymmetricMatrix instance", sm instanceof SymmetricMatrix);
		
		Matrix invSM = sm.getInverseMatrix();
		Matrix smTimesInvSM = sm.multiply(invSM);
		
		Assert.assertTrue("Checking if the in-class getInverseMatrix method produces a SymmetricMatrix instance", invSM instanceof SymmetricMatrix);

		Assert.assertTrue("Checking if the product of the inverse by the original matrix returns the identity matrix", smTimesInvSM.equals(Matrix.getIdentityMatrix(sm.m_iRows)));
	}

	@Test
	public void testIsserlisMatrix() {
		Matrix m = new Matrix(3,3);
		m.setValueAt(0, 0, 1);
		m.setValueAt(0, 1, 2);
		m.setValueAt(0, 2, 3);
		m.setValueAt(1, 0, 2);
		m.setValueAt(1, 1, 4);
		m.setValueAt(1, 2, 5);
		m.setValueAt(2, 0, 3);
		m.setValueAt(2, 1, 5);
		m.setValueAt(2, 2, 6);
		
		Assert.assertTrue("Checking if the m instance is symmetric", m.isSymmetric());
		
		Matrix isserlis = m.getIsserlisMatrixOnlyForTestPurpose();
//		System.out.println(isserlis);
		Assert.assertTrue("Checking if the m.getIsserlisMatrix method produces a Matrix instance that is is symmetric", isserlis.isSymmetric());
		
		Matrix isserlis2 = SymMat.getIsserlisMatrix();
		Assert.assertTrue("Checking if the SymMat.getIsserlisMatrix method produces a SymmetricMatrix instance", isserlis2 instanceof SymmetricMatrix);

		Assert.assertTrue("Checking if the two matrices are equal", isserlis.equals(isserlis2));

	}
	
	@Test
	public void symSquare_SquareSymTest() {
		Matrix symSquareVector = SymMat.symSquare();
		Assert.assertEquals("Testing the vector length.", 6, symSquareVector.getNumberOfElements());
		SymmetricMatrix originalMatrix = symSquareVector.squareSym();
		Assert.assertTrue("Testing if the recomposed matrix is equal to the original", SymMat.equals(originalMatrix));
	}
	
}
