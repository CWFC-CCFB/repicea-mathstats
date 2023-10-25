/*
 * This file is part of the repicea library.
 *
 * Copyright (C) 2009-2019 Mathieu Fortin for Rouge-Epicea
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

import java.io.Serializable;
import java.security.InvalidParameterException;
import java.text.DecimalFormat;
import java.text.NumberFormat;
import java.util.ArrayList;
import java.util.Collections;
import java.util.List;

import repicea.util.DeepCloneable;


/**
 * This class implement most of the basic function in linear algebra
 * Authors: Jean-Francois Lavoie and Mathieu Fortin (June 2009)
 */
public class Matrix implements Serializable, DeepCloneable {

	private static final long serialVersionUID = 20100804L;
		
	protected static int SizeBeforeSwitchingToLUDecompositionInDeterminantCalculation = 6;
	
	public static int NB_ROWS_BEYOND_WHICH_MATRIX_INVERSION_TAKES_TOO_MUCH_TIME = 600;
	
	private static final NumberFormat ScientificFormatter = new DecimalFormat("0.##E0");
	
	private static final NumberFormat SimpleDecimalFormatter = NumberFormat.getNumberInstance();
	static {SimpleDecimalFormatter.setMinimumFractionDigits(2);}
	
	/*
	 * Members of this class
	 */
	private static final double VERY_SMALL = 1.5E-06;
	
	private static final double EPSILON = 1E-12;
	
	final double[][] m_afData;
	public final int m_iRows;
	public final int m_iCols;
	
	/**
	 * Constructor 1. Creates a matrix from a two-dimension array.
	 * @param data a two-dimension array of double
	 */
	public Matrix(double data[][]) {
		this(data.length, data[0].length);
		for (int i = 0; i < m_iRows; i++)
			for (int j = 0; j < m_iCols; j++)
				setValueAt(i, j, data[i][j]);
	}
	
	/**
	 * Constructor 2. Creates a column vector from an array of double
	 * @param data an array of double instances.
	 */
	public Matrix(double data[]) {
		this(data.length, 1);
		for (int i = 0; i < m_iRows; i++)
			setValueAt(i, 0, data[i]);
	}
	
	/**
	 * Constructor 3. Creates a column vector with all the values found in the List instance.
	 * @param list a List of Number-derived instances
	 */
	public Matrix(List<? extends Number> list) {
		this(list.size(), 1);
		Number number;
		for (int i = 0; i < m_iRows; i++) {
			number = list.get(i);
			setValueAt(i, 0, number.doubleValue());
		}
	}
	
	/**
	 * Constructor 4. Creates a matrix with the elements starting from a given number with a particular increment.
	 * @param iRows number of rows
	 * @param iCols number of columns
	 * @param from first element of the matrix
	 * @param iIncrement increment for the next elements
	 */
	public Matrix(int iRows, int iCols, double from, double iIncrement) {
		this(iRows,iCols);
		double value = from;
		for (int i = 0; i < m_iRows; i++) {
			for (int j = 0; j < m_iCols; j++) {
				setValueAt(i, j, value);
				value += iIncrement;
			}
		}
	}

	/**
	 * Basic constructor. Creates a matrix with all elements set to 0.
	 * @param iRows number of rows
	 * @param iCols number of columns
	 */
	public Matrix(int iRows, int iCols) {
		this(iRows, iCols, true);
	}

	/**
	 * Protected constructor which allows for the former implementation. 
	 * @param iRows number of rows
	 * @param iCols number of columns
	 * @param newImplementation if true, the column vector are stored in transposed manner so that they take less memory. 
	 * By default, it is set to true.
	 */
	protected Matrix(int iRows, int iCols, boolean newImplementation) {
		if (iRows <= 0 || iCols <= 0) {
			throw new InvalidParameterException("The number of rows or columns must be equal to or greater than 1!");
		}
		if (iCols == 1 && newImplementation) {
//			m_afData = new double[iCols][iRows];		// the array is stored as a row vector for better memory management
			m_afData = contructInternalArray(iCols, iRows);		// the array is stored as a row vector for better memory management
		} else {
//	  		m_afData = new double[iRows][iCols];
			m_afData = contructInternalArray(iRows, iCols);
		}
		m_iRows = iRows;
		m_iCols = iCols;
	}

	
	protected double[][] contructInternalArray(int iRows, int iCols) {
		return new double[iRows][iCols];
	}
	
	/**
	 * Create an array from the Matrix instance. Note that this array is not the internal array. It is a copy and consequently repeated calls to 
	 * this method with large matrices might be computationally intensive.
	 * @return a 2-dimension array
	 */
	public final double[][] toArray() {
		double[][] arr = new double[m_iRows][m_iCols];
		for (int i = 0; i < m_iRows; i++) {
			for (int j = 0; j < m_iCols; j++) {
				arr[i][j] = getValueAt(i, j);
			}
		}
		return arr;
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
	public void setValueAt(int i, int j, double value) {
		if (isNewImplementationForColumnVector()) {	// the vector is actually transposed for a better memory management
			m_afData[j][i] = value;
		} else {
			m_afData[i][j] = value;
		}
	}

	/**
	 * Return the value at row i and column j.
	 * @param i the row index
	 * @param j the column index
	 * @return the entry
	 */
	public double getValueAt(int i, int j) {
		return isNewImplementationForColumnVector() ? m_afData[j][i] : m_afData[i][j];
	}
	
	/**
	 * This method add matrix m to the current matrix.
	 * @param m the matrix to be added
	 * @return the result in a new Matrix instance
	 */
	public Matrix add(Matrix m) {
		if (!isTheSameDimension(m)) {
			throw new UnsupportedOperationException("This instance and the Matrix m are not of the same dimension!");
		}
		Matrix mat = new Matrix(m_iRows, m_iCols);
		for (int i = 0; i < m_iRows; i++) {
			for (int j = 0; j < m_iCols; j++) {
				mat.setValueAt(i, j, getValueAt(i, j) + m.getValueAt(i, j));
			}
		}
		
		return mat;
	}

	/**
	 * This method tests whether if any element of the Matrix object is 
	 * different from parameter d. 
	 * @param d the value to be checked
	 * @return true if at least one element is different from d
	 */
	public final boolean anyElementDifferentFrom(double d) {
		boolean bool = false;
		for (int i = 0; i < m_iRows; i++) {
			for (int j = 0; j < m_iCols; j++) {
				if (Math.abs(getValueAt(i, j) - d) > EPSILON)
					bool = true;
			}
		}
		return bool;
	}

	/**
	 * This method tests whether if any element of the Matrix object is 
	 * larger than parameter d. 
	 * @param d the value to be checked
	 * @return true if at least one element is larger than d
	 */
	public final boolean anyElementLargerThan(double d) {
		for (int i = 0; i < this.m_iRows; i++) {
			for (int j = 0; j < this.m_iCols; j++) {
				if (getValueAt(i, j) > d) {
					return true;
				}
			}
		}
		return false;
	}

	
	/**
	 * This method tests whether if any element of the Matrix object is 
	 * smaller than or equal to parameter d. 
	 * @param d the value to be checked
	 * @return true if at least one element is larger than d
	 */
	public final boolean anyElementSmallerOrEqualTo(double d) {
		for (int i = 0; i < this.m_iRows; i++) {
			for (int j = 0; j < this.m_iCols; j++) {
				if (getValueAt(i, j) <= d) {
					return true;
				}
			}
		}
		return false;
	}

	/**
	 * Check if any element is a "Not a Number" (NaN)
	 * @return true if at least one element is NaN
	 */
	public boolean anyElementNaN() {
		for (int i = 0; i < this.m_iRows; i++) {
			for (int j = 0; j < this.m_iCols; j++) {
				if (Double.isNaN(this.getValueAt(i, j)))
					return true;
			}
		}
		return false;
	}

	
	/**
	 * This method return a vector that contains the diagonal element of this Matrix instance.
	 * A check is implemented to make sure this is a square matrix.
	 * @return the resulting matrix
	 * @throws UnsupportedOperationException if the matrix is not square
	 */
	public final Matrix diagonalVector() {
		if (!isSquare()) {
			throw new UnsupportedOperationException("Matrix.diagonalVector() : The input matrix is not square");
		} else  {
			Matrix oMat = new Matrix(m_iRows, 1);
			for (int i = 0; i < m_iRows; i++) {
				oMat.setValueAt(i, 0, getValueAt(i, i));
			}
			return oMat;
		}
	}

	/**
	 * This method returns true if the matrix contains at least one value NaN.
	 * @return a boolean
	 */
	public final boolean doesContainAnyNaN() {
		for (int i = 0; i < m_iRows; i++) {
			for (int j = 0; j < m_iCols; j++) {
				if (Double.isNaN(getValueAt(i, j))) {
					return true;
				}
			}
		}
		return false;
	}

	/**
	 * Compute the elementwise division of this by m.
	 * @param m a Matrix instance of the same dimensions
	 * @return the resulting matrix
	 */
	public Matrix elementWiseDivide(Matrix m) {
		if (isTheSameDimension(m)) {
			Matrix oMat = new Matrix(this.m_iRows,this.m_iCols);
			for (int i = 0; i < this.m_iRows; i++) {
				for (int j = 0; j < this.m_iCols; j++) {
					oMat.setValueAt(i, j, getValueAt(i, j) / m.getValueAt(i, j));
				}
			}
			return oMat;
		} else {
			throw new UnsupportedOperationException("The matrix m does not have the same dimensions than the current matrix!");
		}
	}
	
	/**
	 * This method compute the elementwise product of this x m
	 * @param m the matrix that contains the elements to be multiplied with.
	 * @return a Matrix instance
	 */
	public Matrix elementWiseMultiply(Matrix m) {
		if (isTheSameDimension(m)) {
			Matrix oMat = new Matrix(this.m_iRows,this.m_iCols);
			for (int i = 0; i < this.m_iRows; i++) {
				for (int j = 0; j < this.m_iCols; j++) {
					if (getValueAt(i, j) != 0d && m.getValueAt(i, j) != 0d) {
						oMat.setValueAt(i, j, getValueAt(i, j) * m.getValueAt(i, j));
					}
				}
			}
			return oMat;
		} else {
			throw new UnsupportedOperationException("The matrix m does not have the same dimensions than the current matrix!");
		}
	}
	
	/**
	 * Compute the exponential of the elements of this matrix.
	 * @return the results in a matrix
	 */
	public Matrix expMatrix() {
		Matrix matrix = new Matrix(m_iRows, m_iCols);
		for (int i = 0; i < matrix.m_iRows; i++) {
			for (int j = 0; j < matrix.m_iCols; j++) {
				matrix.setValueAt(i, j, Math.exp(getValueAt(i, j)));
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
	public final Matrix getSubMatrix(int startRow, int endRow, int startColumn, int endColumn) {
		int iRows = endRow - startRow + 1;
		int iCols = endColumn - startColumn + 1;
		Matrix mat = new Matrix(iRows, iCols);
		for (int i = 0; i < iRows; i++) {
			for (int j = 0; j < iCols; j++) {
				mat.setValueAt(i, j, getValueAt(startRow + i, startColumn + j));
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
	public final Matrix getSubMatrix(List<Integer> rowIndex, List<Integer> columnIndex, boolean sortIndices) { 
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
		
		Matrix outputMatrix = new Matrix(rowIndex.size(), columnIndex.size());
		for (int i = 0; i < rowIndex.size(); i++) {
			for (int j = 0; j < columnIndex.size(); j++) {
				outputMatrix.setValueAt(i, j, getValueAt(rowIndex.get(i), columnIndex.get(j)));
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
	 * @return a Matrix instance
	 */
	public final Matrix getSubMatrix(List<Integer> rowIndex, List<Integer> columnIndex) { // it should not be sorted here?????
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
		
		Matrix outputMatrix = new Matrix(rowIndex.size(), columnIndex.size());
		for (int i = 0; i < rowIndex.size(); i++) {
			for (int j = 0; j < columnIndex.size(); j++) {
				outputMatrix.setValueAt(i, j, getValueAt(rowIndex.get(i), columnIndex.get(j)));
			}
		}
	
		return outputMatrix;
	}

	/**
	 * This method checks if this is a column vector
	 * @return a boolean that is true if this is a column vector
	 */
	public final boolean isColumnVector() {return m_iCols == 1;}
	
	/**
	 * This method checks whether this matrix is a diagonal matrix, i.e. with all its off-diagonal 
	 * elements being equal to zero.
	 * @return a boolean
	 */
	public boolean isDiagonalMatrix() {
		if (isSquare()) {
			for (int i = 0; i < m_iRows; i++) {
				for (int j = 0; j < m_iCols; j++) {
					if (i != j) {
						if (Math.abs(getValueAt(i, j)) != 0d) {
							return false;
						}
					}
				}
			}
			return true;
		} else {
			return false;
		}
	}

	
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
	 * This method tests whether the matrix is symmetric. 
	 * @return true if the matrix is symmetric or false otherwise
	 */
	public boolean isSymmetric() {
		boolean valid = true;
		if (!isSquare()) {
			return false;
		}
		
		outerLoop:
		for (int i = 0; i < m_iRows; i++) {
			for (int j = i + 1; j < m_iCols; j++) {
				if (Math.abs(getValueAt(j, i)) < 1E-50) {		// equal to 0
					if (Math.abs(getValueAt(i, j)) > 1E-50) {	// not equal to 0
						valid = false;
						break outerLoop;
					}
				} else {
					double ratio = getValueAt(i, j) / getValueAt(j, i);
					if (Math.abs(ratio - 1) > VERY_SMALL) {
						valid = false;
						break outerLoop;
					}
				}
			}
		}
		return valid;
	}

	/**
	 * This method checks whether or not this and m have the same dimensions
	 * @param m a Matrix instance
	 * @return boolean
	 */
	public final boolean isTheSameDimension(Matrix m) {
		boolean output = false;
		if (m_iCols == m.m_iCols) {
			if (m_iRows == m.m_iRows) {
				output = true;
			}
		}
		return output;
	}

	/**
	 * Compute the logarithm of the elements of this matrix 
	 * @return the results in a Matrix instance
	 * @throws UnsupportedOperationException if one element of the matrix is smaller than or equal to 0
	 */
    public Matrix logMatrix() {
    	boolean valid = true;		// default is valid
    	Matrix matrix = new Matrix(m_iRows, m_iCols);
    	outerloop:
    		for (int i = 0; i < matrix.m_iRows; i++) {
    			for (int j = 0; j < matrix.m_iCols; j++) {
    				if (getValueAt(i, j) <= 0d) {
    					valid = false;
    					break outerloop;
    				}
    				matrix.setValueAt(i, j, Math.log(getValueAt(i, j)));
    			}
    		}
    	if (valid) {
    		return matrix;
    	} else {
    		throw new UnsupportedOperationException("Matrix.logMatrix() : At least one argument value for the log function is smaller or equal to 0");
    	}
    }

	/**
	 * Create a new matrix in which the current matrix represents the first diagonal block and matrix m represents the second
	 * diagonal block.
	 * @param m the matrix to be diagonally blocked
	 * @return the result in a new Matrix instance
	 */
	public Matrix matrixDiagBlock(Matrix m) {
		int m1Row = m_iRows;
		int m1Col = m_iCols;
		int m2Row = m.m_iRows;
		int m2Col = m.m_iCols;
		Matrix matrix = new Matrix(m1Row + m2Row, m1Col + m2Col);
		matrix.setSubMatrix(this, 0, 0);
		matrix.setSubMatrix(m, m1Row, m1Col);
		return matrix;
	}

	
	/**
	 * Compute a diagonal matrix from a row or a column vector.
	 * @return the resulting matrix
	 * @throws UnsupportedOperationException if this matrix is not a vector
	 */
	public final DiagonalMatrix matrixDiagonal() {
		if (!isRowVector() && !isColumnVector()) {
			throw new UnsupportedOperationException("Matrix.matrixDiagonal() : The input matrix is not a vector");
		} else {
			int dim;
			if (isColumnVector()) {
				dim = m_iRows;
			} else {
				dim = m_iCols;
			}
			DiagonalMatrix matrix = new DiagonalMatrix(dim); 
			for (int i = 0; i < dim; i++) {
				if (isColumnVector()) {
					matrix.setValueAt(i, i, getValueAt(i, 0));
				} else {
					matrix.setValueAt(i, i, getValueAt(0, i));
				}
			}
			return matrix;
		}
	}

	/**
	 * Create a new matrix which is the stack of this and matrix m. 
	 * @param m the matrix to stack.
	 * @param stackOver true if the stack is vertically or false if horizontally
	 * @return the stacked matrix
	 */
	public final Matrix matrixStack(Matrix m, boolean stackOver) {
		int m1Row = m_iRows;
		int m1Col = m_iCols;
		int m2Row = m.m_iRows;
		int m2Col = m.m_iCols;
		if (m1Col == m2Col || m1Row == m2Row) {
			if (stackOver) {
				Matrix matrix = new Matrix(m1Row + m2Row, m1Col);
				matrix.setSubMatrix(this, 0, 0);
				matrix.setSubMatrix(m, m1Row, 0);
				return matrix;
			} else {
				Matrix matrix = new Matrix(m1Row,m1Col+m2Col);
				matrix.setSubMatrix(this, 0, 0);
				matrix.setSubMatrix(m, 0, m1Col);
				return matrix;
			}
		} else {
			throw new UnsupportedOperationException("Matrix m cannot be stacked on the current matrix because their dimensions do not match");
		}
	}

	/**
	 * Compute the matrix product of this x m
	 * @param m a Matrix type object 
	 * @return a matrix type object that contains the result of the matrix multiplication
	 */
	public Matrix multiply(Matrix m) {
		if (m_iCols != m.m_iRows) {
			throw new UnsupportedOperationException("The matrix m cannot multiply the current matrix for the number of rows is incompatible!");
		} else {
			Matrix mat = new Matrix(m_iRows, m.m_iCols);
			for (int i_this = 0; i_this < m_iRows; i_this++) {
				for (int j_m = 0; j_m < m.m_iCols; j_m++ ) {
					for (int j_this = 0; j_this < m_iCols; j_this++) {
						int i_m = j_this;
						if (getValueAt(i_this, j_this) != 0d && m.getValueAt(i_m, j_m) != 0d) {
							double newValue = mat.getValueAt(i_this, j_m) + getValueAt(i_this, j_this) * m.getValueAt(i_m, j_m);
							mat.setValueAt(i_this, j_m, newValue);
						}
					}
				}
			}
			return mat;
		}
	}

	/**
	 * Reset all the elements of this Matrix instance to 0.
	 */
	public final void resetMatrix() {
		for (int i = 0; i < m_iRows; i++) {
			for (int j = 0; j < m_iCols; j++) {
				setValueAt(i, j, 0d);
			}
		}
	}
	
	/**
	 * Add the scalar d to all the elements of the current matrix.
	 * @param d the scalar to be added
	 * @return the result in a new Matrix instance
	 */
	public Matrix scalarAdd(double d) {
		Matrix mat = new Matrix(m_iRows, m_iCols);
		for (int i = 0; i < m_iRows; i++) {
			for (int j = 0; j < m_iCols; j++) {
				mat.setValueAt(i, j, getValueAt(i, j) + d);
			}
		}
		return mat;
	}
	
	/**
	 * Multiply the elements of the current matrix by the scalar d.
	 * @param d the multiplier
	 * @return the result in a new Matrix instance
	 */
	public Matrix scalarMultiply(double d) {
		Matrix mat = new Matrix(m_iRows, m_iCols);
		for (int i = 0; i < m_iRows; i++) {
			for (int j = 0; j < m_iCols; j++) {
				mat.setValueAt(i, j, getValueAt(i, j) * d);
			}
		}
		return mat;
	}


	/**
	 * Replace some elements of the matrix by those that are contained in matrix m.
	 * @param m a Matrix instance 
	 * @param i the row index of the first element to be changed
	 * @param j the column index of the first element to be changed
	 */
	public final void setSubMatrix(Matrix m, int i, int j) {
		if (this instanceof SymmetricMatrix) {
			throw new UnsupportedOperationException("This Matrix instance does not support the setSubMatrix method!");
		} else {
			for (int ii = 0; ii < m.m_iRows; ii++) {
				for (int jj = 0; jj < m.m_iCols; jj++) {
					setValueAt(i + ii, j + jj, m.getValueAt(ii, jj));
				}
			}
		}
	}

	
	/**
	 * Create a square symmetric matrix from a vector. <br>
	 * <br>
	 * Checks are implemented to make sure the vector has the appropriate
	 * number of elements.
	 * @return the resulting matrix
	 */
	public final SymmetricMatrix squareSym() {
		if (!isColumnVector()) {
			throw new UnsupportedOperationException("The current matrix is not a column vector!");
		} else {
			double numberElem = m_iRows;
			Double numberRow = (-1.0 + Math.sqrt(1.0+8*numberElem))*0.5;
			int nbRow = numberRow.intValue();
			if (Math.abs(numberRow.doubleValue() - nbRow) > Matrix.VERY_SMALL) {	// check if numberRow is an integer, if not it means the matrix is not square
				throw new UnsupportedOperationException("The number of elements contained in the imput column vector is not appropriate to transform the matrix into a square symmetric matrix!");
			} else {
				SymmetricMatrix matrix = new SymmetricMatrix(nbRow);
				int pointer = 0;
				for (int i = 0; i < nbRow; i++) {
					Matrix subMat = getSubMatrix(pointer, pointer + i, 0, 0);
					for (int j = 0; j < subMat.m_iRows; j++) {
						matrix.setValueAt(j, i, subMat.getValueAt(j, 0));
					}
//					matrix.setSubMatrix(getSubMatrix(pointer, pointer + i, 0, 0), 0, i);
//					matrix.setSubMatrix(getSubMatrix(pointer, pointer + i, 0, 0).transpose(), i, 0);
					pointer += i + 1; 
				}
				return matrix;
			}
		}
	}
	

	/**
	 * Subtract matrix m from this matrix.
	 * @param m the matrix to be subtracted
	 * @return the result in a new Matrix instance
	 */
	public Matrix subtract(Matrix m) {
		if (!isTheSameDimension(m)) {
			throw new UnsupportedOperationException("This instance and the Matrix m are not of the same dimension!");
		}
		Matrix mat = new Matrix(m_iRows, m_iCols);
		for (int i = 0; i < m_iRows; i++) {
			for (int j = 0; j < m_iCols; j++) {
				mat.setValueAt(i, j, getValueAt(i, j) - m.getValueAt(i, j));
			}
		}
		return mat;
	}

	/**
	 * Compute the trace of the matrix.
	 * @return a double
	 */
	public final double getTrace() {
		if (!isSquare()) {
			throw new UnsupportedOperationException("The trace operation requires the matrix to be square!");
		}
		double sum = 0;
		for (int i = 0; i < m_iRows; i++) {
			sum += getValueAt(i, i);
		}
		return sum;
	}

	/**
	 * Create a transposed matrix.
	 * @return the transposed matrix in a new Matrix instance
	 */
	public Matrix transpose() {
		Matrix matrix = new Matrix(m_iCols, m_iRows);
		for (int i = 0; i < m_iRows; i++) {
			for (int j = 0; j < m_iCols; j++) {
				matrix.setValueAt(j, i, getValueAt(i, j));
			}
		}
		return matrix;
	}
	
	/**
	 * Compute the power of the seed by the elements of the matrix. For example, if the first element of this
	 * matrix is 2, the first element of the resulting matrix will be seed ^ 2.
	 * @param seed a double
	 * @return a Matrix instance
	 */
	public Matrix powMatrix(double seed) {
		Matrix matrix = new Matrix(m_iRows, m_iCols);
		for (int i = 0; i < matrix.m_iRows; i++) {
			for (int j = 0; j < matrix.m_iCols; j++) {
				matrix.setValueAt(i, j, Math.pow(seed, getValueAt(i, j)));
			}
		}
		return matrix;
	}
	
	/**
	 * Compute the elements of the matrix at a given power.
	 * @param power a double
	 * @return a Matrix instance
	 */
	public Matrix elementWisePower(double power) {
		Matrix matrix = new Matrix(m_iRows, m_iCols);
		for (int i = 0; i < matrix.m_iRows; i++) {
			for (int j = 0; j < matrix.m_iCols; j++) {
				matrix.setValueAt(i, j, Math.pow(getValueAt(i, j), power));
			}
		}
		return matrix;
	}
	

	/**
	 * Repeat this matrix a given number of times in each dimension.
	 * @param nrow the number of times to repeat in row-wise direction
	 * @param ncol the number of times to repeat in column-wise direction
	 * @return the resulting matrix
	 */
	public final Matrix repeat(int nrow, int ncol) {
		Matrix resultingMatrix = new Matrix(m_iRows * nrow, m_iCols * ncol);
		for (int i = 0; i < nrow; i++) {
			for (int j = 0; j < ncol; j++) {
				resultingMatrix.setSubMatrix(this, i * m_iRows, j * m_iCols);
			}
		}
		return resultingMatrix;
	}
	
	/**
	 * Remove some elements in a particular matrix and create a row vector. 
	 * @param index is the index of the elements to be removed
	 * @return a row vector
	 */
	public final Matrix removeElements(List<Integer> index) {
		Matrix oMat = new Matrix(1, m_iRows * m_iCols - index.size());
		int pointer = 0;
		for (int i=0; i < m_iRows; i++) {
			for (int j=0; j < m_iCols; j++) {
				if (!index.contains(i * m_iCols + j)) {
					oMat.setValueAt(0, pointer, getValueAt(i, j));
					pointer++;
				}
			}
		}
		return oMat;
	}
	
	/**
	 * Return the elements defined by the List indices in a row vector.
	 * @param indices a List of indices
	 * @return a row vector
	 */
	public final Matrix getElements(List<Integer> indices) {
		Matrix oMat = new Matrix(1, indices.size());
		int pointer = 0;
		for (int i = 0; i < m_iRows; i++) {
			for (int j = 0; j < m_iCols; j++) {
				if (indices.contains(i * m_iCols + j)) {
					oMat.setValueAt(0, pointer, getValueAt(i, j));
					pointer++;
				}
			}
		}
		return oMat;
	}
	
	/**
	 * Replace the elements of the matrix designated through the indices by the values
	 * in the row vector m.
	 * @param indices a List of Integer representing the indices
	 * @param m a Matrix instance
	 */
	public final void setElements(List<Integer> indices, Matrix m) {
		if (this instanceof SymmetricMatrix) {
			throw new UnsupportedOperationException("This Matrix instance does not support the setElement method!");
		}
		if (!m.isColumnVector()) {
			throw new InvalidParameterException("Parameter m must be a row vector!");
		}
		for (int i = 0; i < m_iRows; i++) {
			for (int j = 0; j < m_iCols; j++) {
				if (indices.contains(i * m_iCols + j)) {
					setValueAt(i, j, m.getValueAt(indices.indexOf(i * m_iCols + j), 0));
				}
			}
		}
	}

	/**
	 * Add the elements of the parameter matrix to those designated through the indices.
	 * @param indices a List of Integer representing the indices
	 * @param m a Matrix instance
	 */
	public final void addElementsAt(List<Integer> indices, Matrix m) {
		if (this instanceof SymmetricMatrix) {
			throw new UnsupportedOperationException("This Matrix instance does not support the addElementsAt method!");
		}
		if (!m.isColumnVector()) {
			throw new InvalidParameterException("Parameter m must be a row vector!");
		}
		for (int i = 0; i < m_iRows; i++) {
			for (int j = 0; j < m_iCols; j++) {
				if (indices.contains(i * m_iCols + j)) {
					double newValue = getValueAt(i, j) + m.getValueAt(indices.indexOf(i * m_iCols + j), 0);
					setValueAt(i, j, newValue);
				}
			}
		}
	}

	/**
	 * Return a List of integers, which represent the index of the elements
	 * that are equal to parameter d. The index is calculated as i * m_iCols + j.
	 * @param d the value that is checked for
	 * @return a List of integers
	 */
	public final List<Integer> getLocationIndex(double d) {
		List<Integer> list = new ArrayList<Integer>();
		for (int i=0; i < m_iRows; i++) {
			for (int j=0; j < m_iCols; j++) {
				if (Math.abs(getValueAt(i, j) - d) < VERY_SMALL) {
					list.add(i * m_iCols + j);
				}
			}
		}
		return list;
	}

	/**
	 * Compute the sum of all the elements in the Matrix instance.
	 * @return a double
	 */
	public final double getSumOfElements() {
		double sum = 0d;
		for (int i = 0; i < m_iRows; i++) {
			for (int j = 0; j < m_iCols; j++) {
				sum += getValueAt(i, j);
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
	public final double getSumOfElements(int startRow, int endRow, int startColumn, int endColumn) {
		if (endRow >= this.m_iRows || endColumn >= this.m_iCols) {
			throw new InvalidParameterException("The specified end row or end column exceeds the capacity of the matrix!");
		} else if (startRow < 0 || startRow > endRow) {
			throw new InvalidParameterException("The specified start row is either negative or larger than the end row!");
		} else if (startColumn < 0 || startColumn > endColumn) {
			throw new InvalidParameterException("The specified start column is either negative or larger than the end column!");
		}
		double sum = 0d;
		for (int i = startRow; i <= endRow; i++) {
			for (int j = startColumn; j <= endColumn; j++) {
				sum += getValueAt(i, j);
			}
		}
		return sum;
	}
	
	/**
	 * Create a Matrix instance that contains the absolute values
	 * of the original matrix.
	 * @return a Matrix instance that contains the absolute values
	 */
	public Matrix getAbsoluteValue() {
		Matrix oMat = new Matrix(m_iRows, m_iCols);
		for (int i = 0; i < m_iRows; i++) {
			for (int j = 0; j < m_iCols; j++) {
				oMat.setValueAt(i, j, Math.abs(getValueAt(i, j)));
			}
		}
		return oMat;
	}

	/**
	 * Return the number of elements in a Matrix object.
	 * @return the number of elements (integer)
	 */
	public final int getNumberOfElements() {
		return m_iRows * m_iCols;
	}
	
	/**
	 * Calculate the Kronecker product of this by the m Matrix object.
	 * @param m a Matrix instance
	 * @return the resulting product (a matrix instance)
	 */
	public Matrix getKroneckerProduct(Matrix m) {
		Matrix result = new Matrix(m_iRows * m.m_iRows, m_iCols * m.m_iCols);
		for (int i1 = 0; i1 < m_iRows; i1++) {
			for (int j1 = 0; j1 < m_iCols; j1++) {
				for (int i2 = 0; i2 < m.m_iRows; i2++) {
					for (int j2 = 0; j2 < m.m_iCols; j2++) {
						result.setValueAt(i1 * m.m_iRows + i2, j1 * m.m_iCols + j2, getValueAt(i1, j1) * m.getValueAt(i2, j2));
					}
				}
			}
		}
		return result;
	}
	
	Matrix getIsserlisMatrixOnlyForTestPurpose() {
		if (!isSymmetric()) {
			throw new UnsupportedOperationException("Matrix.getIsserlisMatrix: this matrix is not symmetric!");
		} else {
			Matrix output = new Matrix(m_iRows * m_iRows, m_iCols * m_iCols);
			double covariance;
			int indexRow;
			int indexCol;
			for (int i = 0; i < m_iRows; i++) {
				for (int j = i; j < m_iCols; j++) {
					for (int iPrime = 0; iPrime < m_iRows; iPrime++) {
						for (int jPrime = 0; jPrime < m_iCols; jPrime++) {
							covariance = getValueAt(i, j) * getValueAt(iPrime, jPrime) + 
									getValueAt(i, iPrime) * getValueAt(j, jPrime) +
									getValueAt(i, jPrime) * getValueAt(j, iPrime);
							indexRow = i * m_iRows + iPrime;
							indexCol = j * m_iCols + jPrime;
							output.setValueAt(indexRow, indexCol, covariance);
							if (indexRow != indexCol) {
								output.setValueAt(indexCol, indexRow, covariance);
							}
						}
					}
				}
			}
			
			return output;
		}
	}


	@Override
	public Matrix getDeepClone() {
		Matrix oMat = new Matrix(m_iRows, m_iCols);
		for (int i = 0; i < m_iRows; i++) {
			for (int j = 0; j < m_iCols; j++) {
				oMat.setValueAt(i, j, getValueAt(i, j));
			}
		}
		return oMat;
	}
	
	/**
	 * Returns a representation of the matrix content.
	 */
	@Override
	public final String toString() {
		String outputString = "{";
		for (int i = 0; i < m_iRows; i ++) {
			outputString += convertArrayToString(i);
			if (i == m_iRows - 1) {
				outputString += "}";
			} else {
				if (isColumnVector()) {
					outputString += ", ";
				} else {
					outputString += ", \n";
				}
			}
			if (outputString.length() > 5000) {
				outputString += "...";
				break;
			}
		}
		return outputString;
	}

	private String convertArrayToString(int i) {
		String outputString = "";
		for (int j = 0; j < m_iCols; j++) {
			if (j > 0) {
				outputString = outputString.concat(" ");
			}
			double absValue = Math.abs(getValueAt(i, j));
			if (absValue > 0.1 && absValue < 1E3) {
				outputString = outputString.concat("[" + SimpleDecimalFormatter.format(getValueAt(i, j)) + "]");
			} else {
				outputString = outputString.concat("[" + ScientificFormatter.format(getValueAt(i, j)) + "]");
			}
		}
		return outputString;
	}
	
	/**
	 * Return the LU decomposition of this Matrix instance. The diagonal in the lower triangle is 1.
	 * @return an array of two matrices, the first and the second being the lower and the upper triangle, respectively
	 * @throws UnsupportedOperationException if the matrix is not square
	 */
	public final Matrix[] getLUDecomposition() {
		Matrix[] outputMatrices = new Matrix[2];
		if (!isSquare()) {
			throw new UnsupportedOperationException("Matrix.getLUDecomposition(): The matrix is not square!");
		} else {
			Matrix l = new Matrix(m_iRows, m_iCols);
			Matrix u = new Matrix(m_iRows, m_iCols);
		    for (int i = 0; i < m_iRows; i++) {
		        l.setValueAt(i, i, 1d);
		        for (int j = i; j < m_iRows; j++) {
		            double sum = 0;
		            for (int s = 0; s <= i - 1; s++) {
		                sum += l.getValueAt(i, s) * u.getValueAt(s, j);
		            }
		            u.setValueAt(i, j, getValueAt(i, j) - sum);
		        }
		        for (int iii = i + 1; iii < m_iRows; iii++) {
		            double sum = 0;
		            for(int s = 0; s <= i - 1; s++) {
		                sum += l.getValueAt(iii, s) * u.getValueAt(s, i);
		            }
		            if (u.getValueAt(i, i) == 0d) {
		            	throw new UnsupportedOperationException("The determinant cannot be calculated because of a division by 0!");
		            }
		            l.setValueAt(iii, i, (getValueAt(iii, i) - sum) / u.getValueAt(i, i));
		        }
		    }
		    outputMatrices[0] = l;
		    outputMatrices[1] = u;
		    return outputMatrices;
		}
	}

	private void swap(int i, int j, boolean columnWise) {
		if (columnWise) {
			if (i >= m_iCols || j >= m_iCols) {
				throw new UnsupportedOperationException("Columns cannot be swapped as their indices are out of bound!");
			} else {
				double d;
				for (int k = 0; k < m_iRows; k++) {
					d = getValueAt(k, i);
					setValueAt(k, i, getValueAt(k, j));
					setValueAt(k, j, d);
				}
			}
		} else {
			if (i >= m_iRows || j >= m_iRows) {
				throw new UnsupportedOperationException("Columns cannot be swapped as their indices are out of bound!");
			} else {
				double d;
				for (int k = 0; k < m_iCols; k++) {
					d = getValueAt(i, k);
					setValueAt(i, k, getValueAt(j, k));
					setValueAt(j, k, d);
				}
			}
		}
	}
	
	protected void swapAlongTheDiagonal(int i, int j) {
		if (!isSquare()) {
			throw new UnsupportedOperationException("The matrix is not square!");
		}
		if (i >= m_iCols || j >= m_iCols) {
			throw new UnsupportedOperationException("The index is out of bound!");
		} else {
			swap(i,j,true);
			swap(j,i,false);
		}
	}
	

	/**
	 * Return a list of indices that define the blocks in the matrix.
	 * @return a List of List of Integers
	 */
	protected List<List<Integer>> getBlockConfiguration() {
		if (!isSquare()) {
			throw new UnsupportedOperationException("The matrix is not square!");
		}
		
		List<Integer> remainingIndex = new ArrayList<Integer>();
		for (int i = 0; i < m_iCols; i++) {
			remainingIndex.add(i);
		}

		List<List<Integer>> blocks = new ArrayList<List<Integer>>();
		
		if (!isSymmetric()) {
			blocks.add(remainingIndex);
			return blocks;
		} else {
			List<Integer> potentialBlock;
			while (!remainingIndex.isEmpty()) {
				int i = 0;
				potentialBlock = new ArrayList<Integer>();
				potentialBlock.add(remainingIndex.get(i));
				while (i < potentialBlock.size()) {
					int indexI = potentialBlock.get(i);
					for (int j = remainingIndex.indexOf(indexI) + 1; j < remainingIndex.size(); j++) {
						int indexJ = remainingIndex.get(j);
						if (getValueAt(indexI, indexJ) != 0) {
							if (!potentialBlock.contains(indexJ)) {
								potentialBlock.add(indexJ);
							}
							for (int k = i + 1; k < j; k++) {
								int indexK = remainingIndex.get(k);
								if (getValueAt(indexK, indexJ) != 0) {
									if (!potentialBlock.contains(indexK)) {
										potentialBlock.add(indexK);
									}
								}							
							}
						}
					}
					i++;
				}
				blocks.add(potentialBlock);
				remainingIndex.removeAll(potentialBlock);
			}
			return blocks;
		}
	}
	
	
	/**
	 * Compute the minor of this matrix, i.e. the determinant of the Matrix that 
	 * contains all the elements of the original matrix except those in row i and column j.
	 * @param i the index of the row to be omitted
	 * @param j the index of the column to be omitted
	 * @return the minor
	 */
	public final double getMinor(int i, int j) {
		if (i >= m_iRows) {
			throw new InvalidParameterException("The index i is not within the bound of this matrix!");
		} else if (j >= m_iCols) {
			throw new InvalidParameterException("The index j is not within the bound of this matrix!");
		} else if (getNumberOfElements() == 1) {
			throw new UnsupportedOperationException("The matrix only has one element!");
		} else {
			Matrix m = new Matrix(m_iRows - 1, m_iCols - 1);
			int index_i = 0;
			for (int ii = 0; ii < m_iRows; ii++) {		// iterate in the current matrix
				if (ii != i) {
					int index_j = 0;
					for (int jj = 0; jj < m_iCols; jj++) {
						if (jj != j) {
							m.setValueAt(index_i, index_j, getValueAt(ii, jj));
							index_j++;
						}
					}
					index_i++;
				}
			}
			return m.getDeterminant();
		}
	}
	
	/**
	 * Compute the cofactor of this matrix with respect to the element i,j.
	 * @param i the row index of the element
	 * @param j the column index of the element
	 * @return the cofactor matrix
	 */
	public final double getCofactor(int i, int j) {
		double minor = getMinor(i,j);
		double multiplicator = 1d;
		if ((i + j) % 2 != 0) {
			multiplicator = -1d;
		}
		return minor * multiplicator;
	}
	
	
	/**
	 * Compute the determinant of this matrix using Laplace's method for small matrices and LU decomposition
	 * for larger matrices.
	 * @return a double
	 * @throws UnsupportedOperationException if the matrix is not square
	 */ 
	public final double getDeterminant() {
		double determinant = 0;
		if (!isSquare()) {
			throw new UnsupportedOperationException("The matrix is not square!");
		} else if (m_iRows == 1) {
			return getValueAt(0, 0);
		} else if (m_iRows == 2) {
			return getValueAt(0, 0) * getValueAt(1, 1) - getValueAt(0, 1) * getValueAt(1, 0);
		} else if (m_iRows <= SizeBeforeSwitchingToLUDecompositionInDeterminantCalculation) {
			for (int j = 0; j < m_iRows; j++) {
				if (getValueAt(0, j) != 0d) {
					determinant += getValueAt(0, j) * getCofactor(0, j);
				}
			}
		} else {
			Matrix triangle = getLUDecomposition()[1];
			determinant = 1d;
			for (int i = 0; i < triangle.m_iRows; i++) {
				determinant *= triangle.getValueAt(i, i);
			}
		}
		return determinant;
	}

	/**
	 * Compute the adjugate matrix of this matrix.
	 * @return a Matrix instance
	 */
	public final Matrix getAdjugateMatrix() {
		Matrix adjugate = new Matrix(m_iRows, m_iCols);
		for (int i = 0; i < m_iRows; i++) {
			for (int j = 0; j < m_iCols; j++) {
				adjugate.setValueAt(j, i, getCofactor(i, j));		// i and j are inversed for adjugate to ensure the transposition
			}
		}
		return adjugate;
	}
	
	protected Matrix getInternalInverseMatrix() {
		if (m_iRows == 1) {
			Matrix output = new Matrix(1,1);
			output.setValueAt(0, 0, 1d / getValueAt(0, 0));
			return output;
		} else if (m_iRows > SizeBeforeSwitchingToLUDecompositionInDeterminantCalculation) {
			int index = m_iCols / 2;
			Matrix output = new Matrix(m_iRows, m_iCols);
			Matrix a = getSubMatrix(0, index - 1, 0, index - 1);
			Matrix b = getSubMatrix(0, index - 1, index, m_iCols - 1);
			Matrix c = getSubMatrix(index, m_iRows - 1, 0, index - 1);
			Matrix d = getSubMatrix(index, m_iRows - 1, index, m_iCols - 1);
			Matrix invD = d.getInternalInverseMatrix();
			Matrix tmp = b.multiply(invD).multiply(c);
			Matrix invComplement = a.subtract(tmp).getInternalInverseMatrix();
			output.setSubMatrix(invComplement, 0, 0);
			output.setSubMatrix(invComplement.multiply(b).multiply(invD).scalarMultiply(-1d), 0, index);
			output.setSubMatrix(invD.multiply(c).multiply(invComplement).scalarMultiply(-1d), index, 0);
			output.setSubMatrix(invD.multiply(c).multiply(invComplement).multiply(b).multiply(invD).add(invD), index, index);
			return output;
		} else {
			double determinant = getDeterminant();
			if (determinant == 0) {
				throw new UnsupportedOperationException("The matrix cannot be inverted as its determinant is equal to 0!");
			} else {
				return getAdjugateMatrix().scalarMultiply(1d / determinant);
			}
		}
	}
	
	
	/**
	 * Ensure that no elements in the matrix are lower than value by setting them to value if needed. <p>
	 * This method changes the matrix in place 
	 * @param value the threshold. Lower values are set to this value.
	 */
	public final void clampIfLowerThan(double value) {
		if (isNewImplementationForColumnVector()) {
			for (int i = 0; i < m_iRows; i++) {
				for (int j = 0; j < m_iCols; j++) {
					this.m_afData[j][i] = this.m_afData[j][i] < value ? value : this.m_afData[j][i];   
				}
			}
		}
		else {
			for (int i = 0; i < m_iRows; i++) {
				for (int j = 0; j < m_iCols; j++) {
					this.m_afData[i][j] = this.m_afData[i][j] < value ? value : this.m_afData[i][j];   
				}
			}
		}
	}
	
	/**
	 * Ensure that no elements in the matrix are higher than value by setting them to value if needed.<p>
	 * This method changes the matrix in place 
	 * @param value the threshold. Higher values are set to this value.
	 */
	public final void clampIfHigherThan(double value) {
		if (isNewImplementationForColumnVector()) {
			for (int i = 0; i < m_iRows; i++) {
				for (int j = 0; j < m_iCols; j++) {
					this.m_afData[j][i] = this.m_afData[j][i] > value ? value : this.m_afData[j][i];   
				}
			}
		}
		else {
			for (int i = 0; i < m_iRows; i++) {
				for (int j = 0; j < m_iCols; j++) {
					this.m_afData[i][j] = this.m_afData[i][j] > value ? value : this.m_afData[i][j];   
				}
			}
		}
	}
	
	/**
	 * Compute the inverse of this matrix.
	 * @return the inverse matrix
	 */
	public Matrix getInverseMatrix() {
		if (isDiagonalMatrix()) {		// procedure for diagonal matrices
			return DiagonalMatrix.forceConversionToDiagonalMatrix(this).getInverseMatrix();
		} 
		List<List<Integer>> indices = getBlockConfiguration();
		if (indices.size() == 1) {
			return getInternalInverseMatrix();
		} else {
			Matrix inverseMatrix = new Matrix(m_iRows, m_iCols);
			for (List<Integer> blockIndex : indices) {
				Matrix invSubMatrix = getSubMatrix(blockIndex, blockIndex).getInternalInverseMatrix();
				for (int i = 0; i < blockIndex.size(); i++) {
					for (int j = 0; j < blockIndex.size(); j++) {
						inverseMatrix.setValueAt(blockIndex.get(i), blockIndex.get(j), invSubMatrix.getValueAt(i, j));
					}
				}
			}
			return inverseMatrix;
		}
		
	}

	/**
	 * Create an identity matrix of dimension i.
	 * @param dim the dimension of the matrix
	 * @return a Matrix instance
	 */
	public static DiagonalMatrix getIdentityMatrix(int dim) {
		if (dim <= 0) {
			throw new InvalidParameterException("The dim argument must be larger than 0!");
		}
		DiagonalMatrix mat = new DiagonalMatrix(dim);
		for (int i = 0; i < dim; i++) {
			mat.setValueAt(i, i, 1d);
		}
		return mat;
	}
	
	@Override
	public final boolean equals(Object obj) {
		if (obj instanceof Matrix) {
			Matrix mat = (Matrix) obj;
			if (mat.m_iCols != m_iCols || mat.m_iRows != m_iRows) {
				return false;
			} else {
				double jLength = -1;
				for (int i = 0; i < m_afData.length; i++) {
					if (jLength == -1) {
						jLength = m_afData[i].length;
					}
					for (int j = 0; j < jLength; j++) {
						if (getValueAt(i, j) != mat.getValueAt(i, j)) {
							return false;
						}
					}
				}
				return true;
//				return !subtract(mat).getAbsoluteValue().anyElementLargerThan(1E-20);
			}
		} else {
			return false;
		}
	}
	
	public static void main(String[] args) {
		Matrix matrix = new Matrix(5,4, 1, 100);
		System.out.println(matrix);
		matrix = new Matrix(4,3, -0.05, 5E-2);
		System.out.println(matrix);
		matrix = new Matrix(5,1, 0, 2);
		System.out.println(matrix);
	}
	
}
