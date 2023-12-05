/*
 * This file is part of the repicea library.
 *
 * Copyright (C) 2023 His Majesty the King in Right of Canada
 * Author: Mathieu Fortin (Canadian Forest Service)
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

import java.security.InvalidParameterException;
import java.util.List;

import repicea.util.DeepCloneable;

/**
 * The AbstractMatrix class is an abstract class that supports the 
 * Matrix and ComplexMatrix classes.
 * @author Mathieu Fortin - November 2023
 *
 * @param <P> an AbstractMatrix class
 */
@SuppressWarnings("rawtypes")
public abstract class AbstractMatrix<P extends AbstractMatrix> implements DeepCloneable {

	public final int m_iRows;
	public final int m_iCols;

	protected AbstractMatrix(int iRows, int iCols) {
		if (iRows <= 0 || iCols <= 0) {
			throw new InvalidParameterException("The number of rows or columns must be equal to or greater than 1!");
		}
		m_iRows = iRows;
		m_iCols = iCols;
	}


	/**
	 * Add matrix m to the current matrix.
	 * @param m the matrix to be added
	 * @return the result in a new Matrix instance
	 */
	public abstract P add(P m);

	
	/**
	 * Compute the exponential of the elements of this matrix.
	 * @return the results in a Matrix instance
	 */
	public abstract P expMatrix();

	
	/**
	 * This method returns a submatrix of this matrix. 
	 * @param startRow the index of the first row (included)
	 * @param endRow the index of the last row (included)
	 * @param startColumn the index of the first column (included)
	 * @param endColumn the index of the last column (included)
	 * @return the submatrix in a Matrix instance
	 */
	public abstract P getSubMatrix(int startRow, int endRow, int startColumn, int endColumn);

	/**
	 * This method returns a sub matrix whose elements correspond to the indices listed in 
	 * the row index list and the column index list.
	 * 
	 * @param rowIndex a List of integers (if null all the rows are selected)
	 * @param columnIndex a List of integers (if null all the columns are selected)
	 * @param sortIndices a boolean true to enable the sorting of the indices
	 * @return a Matrix instance
	 */
	public abstract P getSubMatrix(List<Integer> rowIndex, List<Integer> columnIndex, boolean sortIndices);

	
	/**
	 * This method returns a sub matrix whose elements correspond to the indices listed in 
	 * the row index list and the column index list. <p>
	 *
	 * This method sorts the indices before constructing the sub matrices. So if rowIndex = {1,3,2},
	 * the rows of resulting submatrix will correspond to rows 1, 2, 3 in this order. It is a proxy for 
	 * getSubMatrix(rowIndex, columnIndex, true). 
	 *  
	 * @see Matrix#getSubMatrix(List, List, boolean)
	 * @param rowIndex a List of integers (if null all the rows are selected)
	 * @param columnIndex a List of integers (if null all the columns are selected)
	 * @return a Matrix instance
	 */
	public abstract P getSubMatrix(List<Integer> rowIndex, List<Integer> columnIndex);

	/**
	 * Compute the elements of the matrix at a given power.
	 * @param power a double
	 * @return a Matrix instance
	 */
	public abstract P elementWisePower(double power);
		
	
	/**
	 * This method compute the elementwise product of this x m
	 * @param m the matrix that contains the elements to be multiplied with.
	 * @return a Matrix instance
	 */
	public abstract P elementWiseMultiply(P m);
	
	/**
	 * This method checks if this is a column vector
	 * @return a boolean that is true if this is a column vector
	 */
	public final boolean isColumnVector() {return m_iCols == 1;}
	
	/**
	 * This method checks if this is a row vector
	 * @return a boolean that is true if this is a row vector
	 */
	public final boolean isRowVector() {return m_iRows == 1;}
	
	/**
	 * This method checks if this is a square matrix
	 * @return true if the matrix is square or false otherwise
	 */
	public final boolean isSquare() {return m_iRows == m_iCols;}
	
	/**
	 * This method checks whether or not this and m have the same dimensions
	 * @param m an AbstractMatrix-derived instance
	 * @return boolean
	 */
	public final boolean isTheSameDimension(P m) {
		boolean output = false;
		if (m_iCols == m.m_iCols) {
			if (m_iRows == m.m_iRows) {
				output = true;
			}
		}
		return output;
	}

	/**
	 * Add the scalar d to all the elements of the current matrix.
	 * @param d the scalar to be added
	 * @return the result in a new Matrix instance
	 */
	public abstract P scalarAdd(double d);
	
	/**
	 * Multiply the elements of the current matrix by the scalar d.
	 * @param d the multiplier
	 * @return the result in a new Matrix instance
	 */
	public abstract P scalarMultiply(double d);


	/**
	 * Replace some elements of the matrix by those that are contained in matrix m.
	 * @param m a Matrix instance 
	 * @param i the row index of the first element to be changed
	 * @param j the column index of the first element to be changed
	 */
	public abstract void setSubMatrix(P m, int i, int j);


	/**
	 * Compute the logarithm of the elements of this matrix 
	 * @return the results in a Matrix instance
	 * @throws UnsupportedOperationException if one element of the matrix is smaller than or equal to 0
	 */
    public abstract P logMatrix();
    
	/**
	 * Subtract matrix m from this matrix.
	 * @param m the matrix to be subtracted
	 * @return the result in a new Matrix instance
	 */
	public abstract P subtract(P m);
	
	
	/**
	 * Create a transposed matrix.
	 * @return the transposed matrix in a new Matrix instance
	 */
	public abstract P transpose();


	@Override
	public abstract P getDeepClone();


	/**
	 * Returns a representation of the matrix content.
	 */
	@Override
	public final String toString() {
		StringBuilder outputString = new StringBuilder();
		outputString.append("{");
		for (int i = 0; i < m_iRows; i ++) {
			outputString.append(convertArrayToString(i));
			if (i == m_iRows - 1) {
				outputString.append("}");
			} else {
				if (isColumnVector()) {
					outputString.append(", ");
				} else {
					outputString.append(", \n");
				}
			}
			if (outputString.length() > 5000) {
				outputString.append("...");
				break;
			}
		}
		return outputString.toString();
	}

	/**
	 * Convert a particular row of the matrix into a string.
	 * @param rowIndex the index of the row to be converted
	 * @return a String instance
	 */
	protected abstract String convertArrayToString(int rowIndex);

	
	
}
