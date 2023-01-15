/*
 * This file is part of the repicea-statistics library.
 *
 * Copyright (C) 2009-2012 Mathieu Fortin for Rouge-Epicea
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

import java.io.File;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

import org.junit.Assert;
import org.junit.Test;

import repicea.io.FormatField;
import repicea.io.javacsv.CSVField;
import repicea.io.javacsv.CSVWriter;
import repicea.lang.REpiceaSystem;
import repicea.serial.xml.XmlDeserializer;
import repicea.serial.xml.XmlMarshallException;
import repicea.serial.xml.XmlSerializer;
import repicea.stats.StatisticalUtility;
import repicea.stats.StatisticalUtility.TypeMatrixR;
import repicea.util.ObjectUtility;

/**
 * This test class performs some tests on matrix calculation.
 * @author Mathieu Fortin - September 2013
 */
public class MatrixTest {

	/**
	 * This test is performed on the calculation of the inverse of a blocked matrix.
	 */
	@Test
	public void blockedInversedMatrixTest() {
		
		Matrix mat = new Matrix(9,9);
		mat.setValueAt(0, 0, 5.49);
		mat.setValueAt(0, 4, 1.85);
		mat.setValueAt(1, 1, 3.90);
		mat.setValueAt(2, 2, 2.90);
		mat.setValueAt(2, 3, 1.02);
		mat.setValueAt(2, 5, 0.70);
		mat.setValueAt(2, 6, 0.76);
		mat.setValueAt(2, 7, 0.77);
		mat.setValueAt(2, 8, 0.80);
		mat.setValueAt(3, 3, 3.20);
		mat.setValueAt(3, 5, 0.89);
		mat.setValueAt(3, 6, 0.87);
		mat.setValueAt(3, 7, 0.89);
		mat.setValueAt(3, 8, 0.93);
		mat.setValueAt(4, 4, 4.55);
		mat.setValueAt(5, 5, 2.70);
		mat.setValueAt(5, 6, 0.66);
		mat.setValueAt(5, 7, 0.67);
		mat.setValueAt(5, 8, 0.70);
		mat.setValueAt(6, 6, 2.69);
		mat.setValueAt(6, 7, 0.66);
		mat.setValueAt(6, 8, 0.69);
		mat.setValueAt(7, 7, 2.70);
		mat.setValueAt(7, 8, 0.70);
		mat.setValueAt(8, 8, 2.76);
		
		for (int i = 0; i < mat.m_iRows; i++) {
			for (int j = i; j < mat.m_iCols; j++) {
				if (i != j) {
					mat.setValueAt(j, i, mat.getValueAt(i, j));
				}
			}
		}
		
		Matrix invMat = mat.getInverseMatrix();
		Matrix ident = mat.multiply(invMat);
		Matrix diff = ident.subtract(Matrix.getIdentityMatrix(ident.m_iCols)).getAbsoluteValue();
		boolean equalToIdentity = !diff.anyElementLargerThan(1E-15);
		Assert.assertEquals(true, equalToIdentity);
	}

	
	/**
	 * This test is performed on the calculation of the inverse of a blocked matrix.
	 */
	@Test
	public void toArrayTest() {
		
		Matrix mat = new Matrix(9,9);
		mat.setValueAt(0, 0, 5.49);
		mat.setValueAt(0, 4, 1.85);
		mat.setValueAt(1, 1, 3.90);
		mat.setValueAt(2, 2, 2.90);
		mat.setValueAt(2, 3, 1.02);
		mat.setValueAt(2, 5, 0.70);
		mat.setValueAt(2, 6, 0.76);
		mat.setValueAt(2, 7, 0.77);
		mat.setValueAt(2, 8, 0.80);
		mat.setValueAt(3, 3, 3.20);
		mat.setValueAt(3, 5, 0.89);
		mat.setValueAt(3, 6, 0.87);
		mat.setValueAt(3, 7, 0.89);
		mat.setValueAt(3, 8, 0.93);
		mat.setValueAt(4, 4, 4.55);
		mat.setValueAt(5, 5, 2.70);
		mat.setValueAt(5, 6, 0.66);
		mat.setValueAt(5, 7, 0.67);
		mat.setValueAt(5, 8, 0.70);
		mat.setValueAt(6, 6, 2.69);
		mat.setValueAt(6, 7, 0.66);
		mat.setValueAt(6, 8, 0.69);
		mat.setValueAt(7, 7, 2.70);
		mat.setValueAt(7, 8, 0.70);
		mat.setValueAt(8, 8, 2.76);
		
		double[][] array = mat.toArray();
		Matrix mat2 = new Matrix(array);
		boolean isThereAnyDiff = mat.subtract(mat2).getAbsoluteValue().anyElementLargerThan(1E-8);
		Assert.assertTrue("Testing that there is no difference between the two matrices", !isThereAnyDiff);
	}
	
	/**
	 * This test is performed on the calculation of the inverse of a large symmetric matrix with many zero cells.
	 */
	@Test
	public void inversionWithZeroCellsTest() {
		
		Matrix coordinates = new Matrix(20,1,0,1);
		
//		Matrix rMatrix = StatisticalUtility.constructRMatrix(coordinates, 2, 0.2, TypeMatrixR.LINEAR);
		Matrix rMatrix = StatisticalUtility.constructRMatrix(Arrays.asList(new Double[] {2d, 0.2}), TypeMatrixR.LINEAR, coordinates);
		Matrix invMat = rMatrix.getInverseMatrix();
	
		Matrix ident = rMatrix.multiply(invMat);

		Matrix diff = ident.subtract(Matrix.getIdentityMatrix(ident.m_iCols)).getAbsoluteValue();

		boolean equalToIdentity = !diff.anyElementLargerThan(1E-10);
		
		Assert.assertEquals(true, equalToIdentity);
	}

	@Test
	public void testMemoryManagement() {
		List<Matrix> myArrayList = new ArrayList<Matrix>();
		for (int i = 0; i < 10000; i++) {
			myArrayList.add(new Matrix(700,1, false)); // old implementation
		}
		double currentMemoryLoad = REpiceaSystem.getCurrentMemoryLoadMb();
		System.out.println("Current memory load with old implementation = " + currentMemoryLoad + " Mb");
		myArrayList.clear();
		for (int i = 0; i < 10000; i++) {
			myArrayList.add(new Matrix(700,1)); // new implementation
		}
		double newMemoryLoad = REpiceaSystem.getCurrentMemoryLoadMb();
		System.out.println("Current memory load with new implementation = " + newMemoryLoad + " Mb");
		Assert.assertTrue("Testing that former implementation takes at least twice the memory space", currentMemoryLoad > newMemoryLoad * 2);
	}
	
		@Test
	public void deserializationOfAMatrix() throws FileNotFoundException, XmlMarshallException {
		String pathname = ObjectUtility.getPackagePath(getClass()) + "serializedMatrix.xml";
//		Matrix myMatrix = new Matrix(2,4);
//		XmlSerializer serializer = new XmlSerializer(pathname);
//		serializer.writeObject(myMatrix);
		
		XmlDeserializer deserializer = new XmlDeserializer(pathname);
		Matrix myDeserializedMatrix = (Matrix) deserializer.readObject();

		Assert.assertEquals("Testing deserialized matrix nb rows", 2, myDeserializedMatrix.m_iRows);
		Assert.assertEquals("Testing deserialized matrix nb cols", 4, myDeserializedMatrix.m_iCols);
	}
	
	@Test
	public void serializationWithAndWithoutCompression() throws XmlMarshallException {
		Matrix mat = new Matrix(100,100);
		String filename1 = ObjectUtility.getPackagePath(getClass()) + "serializedWithCompression.xml";
		XmlSerializer ser1 = new XmlSerializer(filename1);
		ser1.writeObject(mat);
		File file1 = new File(filename1);
		long file1size = file1.length();
		String filename2 = ObjectUtility.getPackagePath(getClass()) + "serializedWithoutCompression.xml";
		XmlSerializer ser2 = new XmlSerializer(filename2, false);
		ser2.writeObject(mat);
		File file2 = new File(filename2);
		long file2size = file2.length();
		double ratio = (double) file2size / file1size; 
		Assert.assertTrue("Testing compression ratio", ratio > 75);
	}

	public void speedTestInversionMatrix(int iMax) throws IOException {
		String filename = ObjectUtility.getPackagePath(getClass()) + "inversionTimes.csv";
		CSVWriter writer = new CSVWriter(new File(filename), false);
		List<FormatField> fields = new ArrayList<FormatField>();
		fields.add(new CSVField("dimension"));
		fields.add(new CSVField("time"));
		writer.setFields(fields);
		
		long startingTime;
		Object[] record = new Object[2];
		int size;
		for (int i = 1; i < iMax; i++) {
			size = i * 10;
			record[0] = size;
			startingTime = System.currentTimeMillis();
			Matrix coordinates = new Matrix(size,1,0,1);
//			Matrix rMatrix = StatisticalUtility.constructRMatrix(coordinates, 2, 0.2, TypeMatrixR.LINEAR);
			Matrix rMatrix = StatisticalUtility.constructRMatrix(Arrays.asList(new Double[] {2d, 0.2}), TypeMatrixR.LINEAR, coordinates);
			rMatrix.getInverseMatrix();
			record[1] = (System.currentTimeMillis() - startingTime) * .001;
			writer.addRecord(record);
		}
		writer.close();
	}
	
	public void speedTestMatrixMultiplication() {
		int i = 1000;
		Matrix oMat = Matrix.getIdentityMatrix(i);
		long startingTime;
		startingTime = System.currentTimeMillis();
		oMat.multiply(oMat);
		System.out.println("Elapsed time = " + ((System.currentTimeMillis() - startingTime) * .001));
	}
	
//	public static void main(String[] args) throws IOException {
//		MatrixTests test = new MatrixTests();
////		test.speedTestInversionMatrix(100);
//		test.speedTestMatrixMultiplication();
//	}
	
}
