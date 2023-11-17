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
import java.util.ArrayList;
import java.util.Collections;
import java.util.List;

/**
 * The ComplexMatrix class implements a special class of 
 * matrices for complex numbers.
 * @author Mathieu Fortin - November 2023
 *
 */
public class ComplexMatrix extends AbstractMatrix<ComplexMatrix> {
	
	final ComplexNumber[][] m_afData;
	
	/**
	 * Constructor 1. Creates a matrix from a two-dimension array.
	 * @param data a two-dimension array of double
	 */
	public ComplexMatrix(ComplexNumber[][] data) {
		this(data.length, data[0].length);
		for (int i = 0; i < m_iRows; i++)
			for (int j = 0; j < m_iCols; j++)
				setValueAt(i, j, data[i][j].clone());
	}
	
	/**
	 * Constructor 2. Creates a column vector from an array of double
	 * @param data an array of double instances.
	 */
	public ComplexMatrix(ComplexNumber[] data) {
		this(data.length, 1);
		for (int i = 0; i < m_iRows; i++)
			setValueAt(i, 0, data[i].clone());
	}
	
	/**
	 * Constructor 3. Creates a column vector with all the values found in the List instance.
	 * @param list a List of Number-derived instances
	 */
	public ComplexMatrix(List<ComplexNumber> list) {
		this(list.size(), 1);
		for (int i = 0; i < m_iRows; i++) {
			setValueAt(i, 0, list.get(i).clone());
		}
	}

	/**
	 * Basic constructor. Creates a matrix with all elements set to 0.
	 * @param iRows number of rows
	 * @param iCols number of columns
	 */
	private ComplexMatrix(int iRows, int iCols) {
		super(iRows, iCols);
		if (iCols == 1) {
			m_afData = contructInternalArray(iCols, iRows);		// the array is stored as a row vector for better memory management
		} else {
			m_afData = contructInternalArray(iRows, iCols);
		}
	}

	protected ComplexNumber[][] contructInternalArray(int iRows, int iCols) {
		return new ComplexNumber[iRows][iCols];
	}
	
	private boolean isNewImplementationForColumnVector() {
		return isColumnVector() && !isRowVector() && m_afData.length == 1;	// second condition is to avoid side effect when dealing with a 1x1 matrix
	}
	
	/**
	 * Set the value at row i and column j.
	 * @param i the row index 
	 * @param j the column index
	 * @param value the value to be set in the cell
	 */
	public void setValueAt(int i, int j, ComplexNumber value) {
		if (isNewImplementationForColumnVector()) {	// the vector is actually transposed for a better memory management
			m_afData[j][i] = value.clone();
		} else {
			m_afData[i][j] = value.clone();
		}
	}

	/**
	 * Return the value at row i and column j.
	 * @param i the row index
	 * @param j the column index
	 * @return the ComplexNumber entry
	 */
	public ComplexNumber getValueAt(int i, int j) {
		return isNewImplementationForColumnVector() ? m_afData[j][i] : m_afData[i][j];
	}
	
	@Override
	public ComplexMatrix add(ComplexMatrix m) {
		if (!isTheSameDimension(m)) {
			throw new UnsupportedOperationException("This instance and the Matrix m are not of the same dimension!");
		}
		ComplexMatrix mat = new ComplexMatrix(m_iRows, m_iCols);
		for (int i = 0; i < m_iRows; i++) {
			for (int j = 0; j < m_iCols; j++) {
				mat.setValueAt(i, j, getValueAt(i, j).add(m.getValueAt(i, j)));
			}
		}
		return mat;
	}

	@Override
	public ComplexMatrix expMatrix() {
		ComplexMatrix matrix = new ComplexMatrix(m_iRows, m_iCols);
		for (int i = 0; i < matrix.m_iRows; i++) {
			for (int j = 0; j < matrix.m_iCols; j++) {
				matrix.setValueAt(i, j, getValueAt(i, j).exp());
			}
		}
		return matrix;
	}
	
	/**
	 * This method returns a submatrix of this matrix. 
	 * @param startRow the index of the first row (included)
	 * @param endRow the index of the last row (included)
	 * @param startColumn the index of the first column (included)
	 * @param endColumn the index of the last column (included)
	 * @return the submatrix in a Matrix instance
	 */
	public final ComplexMatrix getSubMatrix(int startRow, int endRow, int startColumn, int endColumn) {
		int iRows = endRow - startRow + 1;
		int iCols = endColumn - startColumn + 1;
		ComplexMatrix mat = new ComplexMatrix(iRows, iCols);
		for (int i = 0; i < iRows; i++) {
			for (int j = 0; j < iCols; j++) {
				mat.setValueAt(i, j, getValueAt(startRow + i, startColumn + j).clone());
			}
		}
		return mat;		
	}

	/**
	 * This method returns a sub matrix whose elements correspond to the indices listed in 
	 * the row index list and the column index list.
	 * 
	 * @param rowIndex a List of integers (if null all the rows are selected)
	 * @param columnIndex a List of integers (if null all the columns are selected)
	 * @param sortIndices a boolean true to enable the sorting of the indices
	 * @return a Matrix instance
	 */
	public final ComplexMatrix getSubMatrix(List<Integer> rowIndex, List<Integer> columnIndex, boolean sortIndices) { 
		if (rowIndex != null && !rowIndex.isEmpty()) {
			if (sortIndices)
				Collections.sort(rowIndex);
		} else {
			rowIndex = new ArrayList<Integer>();
			for (int i = 0; i < m_iRows; i++) {
				rowIndex.add(i);
			}
		}
		
		if (columnIndex != null && !columnIndex.isEmpty()) {
			if (sortIndices)
				Collections.sort(columnIndex);
		} else {
			columnIndex = new ArrayList<Integer>();
			for (int j = 0; j < m_iCols; j++) {
				columnIndex.add(j);
			}
		}
		
		ComplexMatrix outputMatrix = new ComplexMatrix(rowIndex.size(), columnIndex.size());
		for (int i = 0; i < rowIndex.size(); i++) {
			for (int j = 0; j < columnIndex.size(); j++) {
				outputMatrix.setValueAt(i, j, getValueAt(rowIndex.get(i), columnIndex.get(j)).clone());
			}
		}
	
		return outputMatrix;
	}

	
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
	 * @return a ComplexMatrix instance
	 */
	public final ComplexMatrix getSubMatrix(List<Integer> rowIndex, List<Integer> columnIndex) { 
		if (rowIndex != null && !rowIndex.isEmpty()) {
			Collections.sort(rowIndex);
		} else {
			rowIndex = new ArrayList<Integer>();
			for (int i = 0; i < m_iRows; i++) {
				rowIndex.add(i);
			}
		}
		
		if (columnIndex != null && !columnIndex.isEmpty()) {
			Collections.sort(columnIndex);
		} else {
			columnIndex = new ArrayList<Integer>();
			for (int j = 0; j < m_iCols; j++) {
				columnIndex.add(j);
			}
		}
		
		ComplexMatrix outputMatrix = new ComplexMatrix(rowIndex.size(), columnIndex.size());
		for (int i = 0; i < rowIndex.size(); i++) {
			for (int j = 0; j < columnIndex.size(); j++) {
				outputMatrix.setValueAt(i, j, getValueAt(rowIndex.get(i), columnIndex.get(j)).clone());
			}
		}
	
		return outputMatrix;
	}

//	/**
//	 * This method checks if this is a column vector
//	 * @return a boolean that is true if this is a column vector
//	 */
//	public final boolean isColumnVector() {return m_iCols == 1;}
//	
//	/**
//	 * This method checks if this is a row vector
//	 * @return a boolean that is true if this is a row vector
//	 */
//	public final boolean isRowVector() {return m_iRows == 1;}
//	
//	/**
//	 * This method checks if this is a square matrix
//	 * @return true if the matrix is square or false otherwise
//	 */
//	public final boolean isSquare() {return m_iRows == m_iCols;}
	
//	/**
//	 * This method checks whether or not this and m have the same dimensions
//	 * @param m a ComplexMatrix instance
//	 * @return boolean
//	 */
//	public final boolean isTheSameDimension(ComplexMatrix m) {
//		boolean output = false;
//		if (m_iCols == m.m_iCols) {
//			if (m_iRows == m.m_iRows) {
//				output = true;
//			}
//		}
//		return output;
//	}

	/**
	 * Compute the logarithm of the elements of this matrix 
	 * @return the results in a Matrix instance
	 * @throws UnsupportedOperationException if one element of the matrix is smaller than or equal to 0
	 */
	@Override
    public ComplexMatrix logMatrix() {
    	ComplexMatrix matrix = new ComplexMatrix(m_iRows, m_iCols);
    	for (int i = 0; i < matrix.m_iRows; i++) {
    		for (int j = 0; j < matrix.m_iCols; j++) {
    			matrix.setValueAt(i, j, getValueAt(i, j).log());
    		}
    	}
   		return matrix;
    }

	
    @Override
	public ComplexMatrix scalarAdd(double d) {
		ComplexMatrix mat = new ComplexMatrix(m_iRows, m_iCols);
		for (int i = 0; i < m_iRows; i++) {
			for (int j = 0; j < m_iCols; j++) {
				mat.setValueAt(i, j, getValueAt(i, j).add(d));
			}
		}
		return mat;
	}
	
	/**
	 * Multiply the elements of the current matrix by the scalar d.
	 * @param d the multiplier
	 * @return the result in a new Matrix instance
	 */
    @Override
	public ComplexMatrix scalarMultiply(double d) {
		ComplexMatrix mat = new ComplexMatrix(m_iRows, m_iCols);
		for (int i = 0; i < m_iRows; i++) {
			for (int j = 0; j < m_iCols; j++) {
				mat.setValueAt(i, j, getValueAt(i, j).multiply(d));
			}
		}
		return mat;
	}

    @Override
	public final void setSubMatrix(ComplexMatrix m, int i, int j) {
		for (int ii = 0; ii < m.m_iRows; ii++) {
			for (int jj = 0; jj < m.m_iCols; jj++) {
				setValueAt(i + ii, j + jj, m.getValueAt(ii, jj));
			}
		}
	}

    @Override
	public ComplexMatrix subtract(ComplexMatrix m) {
		if (!isTheSameDimension(m)) {
			throw new UnsupportedOperationException("This instance and the Matrix m are not of the same dimension!");
		}
		ComplexMatrix mat = new ComplexMatrix(m_iRows, m_iCols);
		for (int i = 0; i < m_iRows; i++) {
			for (int j = 0; j < m_iCols; j++) {
				mat.setValueAt(i, j, getValueAt(i, j).subtract(m.getValueAt(i, j)));
			}
		}
		return mat;
	}

    @Override
	public ComplexMatrix transpose() {
		ComplexMatrix matrix = new ComplexMatrix(m_iCols, m_iRows);
		for (int i = 0; i < m_iRows; i++) {
			for (int j = 0; j < m_iCols; j++) {
				matrix.setValueAt(j, i, getValueAt(i, j).clone());
			}
		}
		return matrix;
	}
	

	/**
	 * Compute the sum of all the elements in the Matrix instance.
	 * @return a double
	 */
	public final ComplexNumber getSumOfElements() {
		ComplexNumber sum = new ComplexNumber(0d,0d);
		for (int i = 0; i < m_iRows; i++) {
			for (int j = 0; j < m_iCols; j++) {
				sum = sum.add(getValueAt(i, j));
			}
		}
		return sum;
	}

	/**
	 * Compute the sum of the elements of a submatrix. The submatrix bounds are determined
	 * through the parameters.
	 * @param startRow the index of the starting row
	 * @param endRow the index of the ending row 
	 * @param startColumn the index of the starting column
	 * @param endColumn the index of the ending column
	 * @return the sum (double)
	 */
	public final ComplexNumber getSumOfElements(int startRow, int endRow, int startColumn, int endColumn) {
		if (endRow >= this.m_iRows || endColumn >= this.m_iCols) {
			throw new InvalidParameterException("The specified end row or end column exceeds the capacity of the matrix!");
		} else if (startRow < 0 || startRow > endRow) {
			throw new InvalidParameterException("The specified start row is either negative or larger than the end row!");
		} else if (startColumn < 0 || startColumn > endColumn) {
			throw new InvalidParameterException("The specified start column is either negative or larger than the end column!");
		}
		ComplexNumber sum = new ComplexNumber(0d,0d);
		for (int i = startRow; i <= endRow; i++) {
			for (int j = startColumn; j <= endColumn; j++) {
				sum = sum.add(getValueAt(i, j));
			}
		}
		return sum;
	}
	
	/**
	 * Return the number of elements in a Matrix object.
	 * @return the number of elements (integer)
	 */
	public final int getNumberOfElements() {
		return m_iRows * m_iCols;
	}
	
	@Override
	public ComplexMatrix getDeepClone() {
		ComplexMatrix oMat = new ComplexMatrix(m_iRows, m_iCols);
		for (int i = 0; i < m_iRows; i++) {
			for (int j = 0; j < m_iCols; j++) {
				oMat.setValueAt(i, j, getValueAt(i, j).clone());
			}
		}
		return oMat;
	}
	
//	/**
//	 * Returns a representation of the matrix content.
//	 */
//	@Override
//	public final String toString() {
//		String outputString = "{";
//		for (int i = 0; i < m_iRows; i ++) {
//			outputString += convertArrayToString(i);
//			if (i == m_iRows - 1) {
//				outputString += "}";
//			} else {
//				if (isColumnVector()) {
//					outputString += ", ";
//				} else {
//					outputString += ", \n";
//				}
//			}
//			if (outputString.length() > 5000) {
//				outputString += "...";
//				break;
//			}
//		}
//		return outputString;
//	}
//
//	private String convertArrayToString(int i) {
//		String outputString = "";
//		for (int j = 0; j < m_iCols; j++) {
//			if (j > 0) {
//				outputString = outputString.concat(" ");
//			}
//			double absValue = Math.abs(getValueAt(i, j));
//			if (absValue > 0.1 && absValue < 1E3) {
//				outputString = outputString.concat("[" + SimpleDecimalFormatter.format(getValueAt(i, j)) + "]");
//			} else {
//				outputString = outputString.concat("[" + ScientificFormatter.format(getValueAt(i, j)) + "]");
//			}
//		}
//		return outputString;
//	}
	
	@Override
	public final boolean equals(Object obj) {
		if (obj instanceof ComplexMatrix) {
			ComplexMatrix mat = (ComplexMatrix) obj;
			if (mat.m_iCols != m_iCols || mat.m_iRows != m_iRows) {
				return false;
			} else {
				double jLength = -1;
				for (int i = 0; i < m_afData.length; i++) {
					if (jLength == -1) {
						jLength = m_afData[i].length;
					}
					for (int j = 0; j < jLength; j++) {
						if (!getValueAt(i, j).equals(mat.getValueAt(i, j))) {
							return false;
						}
					}
				}
				return true;
			}
		} else {
			return false;
		}
	}
	
}
