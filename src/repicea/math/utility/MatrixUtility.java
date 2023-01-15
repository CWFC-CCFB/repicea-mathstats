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
package repicea.math.utility;

import java.security.InvalidParameterException;
import java.util.ArrayList;
import java.util.List;

import repicea.math.Matrix;

/**
 * This class provides some special matrix operations. 
 * 
 * @author Mathieu Fortin - May 2012
 */
public class MatrixUtility {
	
//	/**
//	 * Add a matrix to a reference matrix. The result is returned in the reference matrix. NOTE:
//	 * No new matrix instance is created.
//	 * @param a a Matrix instance (the reference)
//	 * @param b a Matrix instance
//	 */
//	public static void add(Matrix a, Matrix b) {
//		if (a.m_iCols != b.m_iCols || a.m_iRows != b.m_iRows) {
//			throw new InvalidParameterException("Matrices a and b are not the same size!");
//		}
//		for (int i = 0 ; i < a.m_iRows; i++) {
//			for (int j = 0; j < a.m_iCols; j++) {
//				a.setValueAt(i, j, a.getValueAt(i, j) + b.getValueAt(i, j));
//			}
//		}
//	}
	
//	/**
//	 * Subtract a matrix from a reference matrix. The result is returned in the reference matrix. NOTE:
//	 * No new matrix instance is created.
//	 * @param a a Matrix instance (the reference)
//	 * @param b a Matrix instance
//	 */
//	public static void subtract(Matrix a, Matrix b) {
//		if (a.m_iCols != b.m_iCols || a.m_iRows != b.m_iRows) {
//			throw new InvalidParameterException("Matrices a and b are not the same size!");
//		}
//		for (int i = 0 ; i < a.m_iRows; i++) {
//			for (int j = 0; j < a.m_iCols; j++) {
//				a.setValueAt(i, j, a.getValueAt(i, j) - b.getValueAt(i, j));
//			}
//		}
//		
//	}
	
//	/**
//	 * Return the element wise product of two matrices. The result is returned in the first matrix. NOTE:
//	 * No new matrix instance is created.
//	 * @param a a Matrix instance (the reference)
//	 * @param b a Matrix instance
//	 */
//	public static void elementWiseMultiply(Matrix a, Matrix b) {
//		if (a.m_iCols != b.m_iCols || a.m_iRows != b.m_iRows) {
//			throw new InvalidParameterException("Matrices a and b are not the same size!");
//		}
//		for (int i = 0 ; i < a.m_iRows; i++) {
//			for (int j = 0; j < a.m_iCols; j++) {
//				a.setValueAt(i, j, a.getValueAt(i, j) * b.getValueAt(i, j));
//			}
//		}
//		
//	}

//	/**
//	 * Multiply all the element of a matrix by a double. The result is returned in the same matrix. NOTE:
//	 * No new matrix instance is created.
//	 * @param a a Matrix instance 
//	 * @param b a double instance
//	 */
//	public static void scalarMultiply(Matrix a, double b) {
//		for (int i = 0 ; i < a.m_iRows; i++) {
//			for (int j = 0; j < a.m_iCols; j++) {
//				a.setValueAt(i, j, a.getValueAt(i, j) * b);
//			}
//		}
//	}

	/**
	 * Perform a special addition in which only the elements different from 0 and 1 
	 * are involved. NOTE: this method is used with SAS output. 
	 * @param originalMatrix the matrix of parameters
	 * @param matrixToAdd the matrix of parameter deviates
	 * @return the new parameters in a new Matrix instance
	 */
	public static Matrix performSpecialAdd(Matrix originalMatrix, Matrix matrixToAdd) {
		Matrix oMat = originalMatrix.getDeepClone();
		List<Integer> oVector = new ArrayList<Integer>();
		oVector.clear();
		
		for (int i = 0; i < originalMatrix.m_iRows; i++) {
			if (oMat.getValueAt(i, 0) != 0d && oMat.getValueAt(i, 0) != 1d) { 
				oVector.add(i);
			}
		}

		if (oVector.size() != matrixToAdd.m_iRows) {
			throw new InvalidParameterException("The number of rows do not match!");
		} else {
			for (int j = 0; j < oVector.size(); j++) {
				double newValue = oMat.getValueAt(oVector.get(j), 0) + matrixToAdd.getValueAt(j, 0);
				oMat.setValueAt(oVector.get(j), 0, newValue);
			}
		}
		return oMat;
	}
	
	/**
	 * Combine two row vectors of dummy variables. Useful for regressions.
	 * @param mat1 the first row vector
	 * @param mat2 the second row vector
	 * @return the resulting matrix
	 */
	public static Matrix combineMatrices(Matrix mat1, Matrix mat2) {
		if (mat1.m_iRows == mat2.m_iRows) {
			int nbCols = mat1.m_iCols * mat2.m_iCols;
			Matrix oMat = new Matrix(mat1.m_iRows, nbCols);
			for (int i = 0; i < mat1.m_iRows; i++) {
				for (int j = 0; j < mat1.m_iCols; j++) {
					for (int j_prime = 0; j_prime < mat2.m_iCols; j_prime++) {
						oMat.setValueAt(i, j*mat2.m_iCols+j_prime, mat1.getValueAt(i, j) * mat2.getValueAt(i, j_prime));
					}
				}
			}
			return oMat;
		} else {
			throw new UnsupportedOperationException("The two matrices do not have the same number of rows!");
		}
	}

}
