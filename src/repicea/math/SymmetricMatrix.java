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

import java.util.List;

@SuppressWarnings("serial")
public class SymmetricMatrix extends Matrix implements UnmodifiableMatrix {

	public SymmetricMatrix(int size) {
		super(size, size, true);
	}
	
	@Override
	protected double[][] contructInternalArray(int iRows, int iCols) {
		double[][] mainArray = new double[iRows][];
		for (int i = 0; i < mainArray.length; i++) {
			mainArray[i] = new double[iCols - i]; 
		}
		return mainArray;
	}

	@Override
	public void setValueAt(int i, int j, double value) {
		if (j >= i) {
			m_afData[i][j - i] = value;
		} else {
			m_afData[j][i - j] = value;
		}
	}
	
	/**
	 * Return the value at row i and column j.
	 * @param i
	 * @param j
	 * @return a double
	 */
	@Override
	public double getValueAt(int i, int j) {
		return j >= i ? m_afData[i][j - i] : m_afData[j][i - j];
	}

	@Override
	public final boolean isSymmetric() {
		return true;
	}

	@Override
	public Matrix elementWiseDivide(Matrix m) {
		if (isTheSameDimension(m) && m.isSymmetric()) {
			SymmetricMatrix oMat = new SymmetricMatrix(m_iRows);
			for (int i = 0; i < this.m_iRows; i++) {
				for (int j = i; j < this.m_iCols; j++) {
					oMat.setValueAt(i, j, getValueAt(i, j) / m.getValueAt(i, j));
				}
			}
			return oMat;
		} else {
			return super.elementWiseDivide(m);
		}
	}

	
	@Override
	public boolean anyElementNaN() {
		for (int i = 0; i < this.m_iRows; i++) {
			for (int j = i; j < this.m_iCols; j++) {
				if (Double.isNaN(this.getValueAt(i, j)))
					return true;
			}
		}
		return false;
	}

	
	@Override
	public Matrix elementWiseMultiply(Matrix m) {
		if (isTheSameDimension(m) && m.isSymmetric()) {
			SymmetricMatrix oMat = new SymmetricMatrix(m_iRows);
			for (int i = 0; i < this.m_iRows; i++) {
				for (int j = i; j < this.m_iCols; j++) {
					if (getValueAt(i, j) != 0d && m.getValueAt(i, j) != 0d) {
						oMat.setValueAt(i, j, getValueAt(i, j) * m.getValueAt(i, j));
					}
				}
			}
			return oMat;
		} else {
			return super.elementWiseMultiply(m);
		}
	}

	@Override
	public final SymmetricMatrix expMatrix() {
		SymmetricMatrix matrix = new SymmetricMatrix(m_iRows);
		for (int i = 0; i < matrix.m_iRows; i++) {
			for (int j = i; j < matrix.m_iCols; j++) {
				matrix.setValueAt(i, j, Math.exp(getValueAt(i, j)));
			}
		}
		return matrix;
	}
	
	@Override
    public SymmetricMatrix logMatrix() {
    	boolean valid = true;		// default is valid
    	SymmetricMatrix matrix = new SymmetricMatrix(m_iRows);
    	outerloop:
    		for (int i = 0; i < matrix.m_iRows; i++) {
    			for (int j = i; j < matrix.m_iCols; j++) {
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


	

	@Override
	public Matrix multiply(Matrix m) {
		if (m_iCols != m.m_iRows) {
			throw new UnsupportedOperationException("The matrix m cannot multiply the current matrix for the number of rows is incompatible!");
		} else {
			if (m.equals(this)) {	// multiplied by itself yields a SymmetricMatrix instance
				SymmetricMatrix mat = new SymmetricMatrix(m_iRows);
				for (int i_this = 0; i_this < m_iRows; i_this++) {
					for (int j_m = i_this; j_m < m.m_iCols; j_m++ ) {
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
			} else {
				return super.multiply(m);
			}
		}
	}

	
	@Override
	public final SymmetricMatrix scalarAdd(double d) {
		SymmetricMatrix mat = new SymmetricMatrix(m_iRows);
		for (int i = 0; i < m_iRows; i++) {
			for (int j = i; j < m_iCols; j++) {
				mat.setValueAt(i, j, getValueAt(i, j) + d);
			}
		}
		return mat;
	}

	@Override
	public SymmetricMatrix scalarMultiply(double d) {
		SymmetricMatrix mat = new SymmetricMatrix(m_iRows);
		for (int i = 0; i < m_iRows; i++) {
			for (int j = i; j < m_iCols; j++) {
				mat.setValueAt(i, j, getValueAt(i, j) * d);
			}
		}
		return mat;
	}

	
	@Override
	public final Matrix add(Matrix m) {
		if (!isTheSameDimension(m)) {
			throw new UnsupportedOperationException("This instance and the Matrix m are not of the same dimension!");
		}
		if (m.isSymmetric()) {
			SymmetricMatrix mat = new SymmetricMatrix(m_iRows);
			for (int i = 0; i < m_iRows; i++) {
				for (int j = i; j < m_iCols; j++) {
					mat.setValueAt(i, j, getValueAt(i, j) + m.getValueAt(i, j));
				}
			}
			return mat;
		} else {
			return super.add(m);
		}
	}

	
	@Override
	public final Matrix subtract(Matrix m) {
		if (!isTheSameDimension(m)) {
			throw new UnsupportedOperationException("This instance and the Matrix m are not of the same dimension!");
		}
		if (m.isSymmetric()) {
			SymmetricMatrix mat = new SymmetricMatrix(m_iRows);
			for (int i = 0; i < m_iRows; i++) {
				for (int j = i; j < m_iCols; j++) {
					mat.setValueAt(i, j, getValueAt(i, j) - m.getValueAt(i, j));
				}
			}
			return mat;
		} else {
			return super.subtract(m);
		}
	}

	/**
	 * Create a vector of the values corresponding to a symmetric matrix.
	 * @return a nx1 Matrix 
	 */
	public final Matrix symSquare() {
		int numberOfElements = (m_iCols + 1) * m_iCols / 2;
		Matrix outputMatrix = new Matrix(numberOfElements, 1);
		int pointer = 0;
		Matrix tmp;
		for (int i = 0; i < m_iCols; i++) {
			tmp = getSubMatrix(i, i, 0, i).transpose();		// transpose required to get a column vector
			outputMatrix.setSubMatrix(tmp, pointer, 0);
			pointer += tmp.getNumberOfElements();
		}
		return outputMatrix; 
	}


	@Override
	public final Matrix transpose() {
		return this;
	}

	@Override
	public final SymmetricMatrix powMatrix(double seed) {
		SymmetricMatrix matrix = new SymmetricMatrix(m_iRows);
		for (int i = 0; i < matrix.m_iRows; i++) {
			for (int j = i; j < matrix.m_iCols; j++) {
				matrix.setValueAt(i, j, Math.pow(seed, getValueAt(i, j)));
			}
		}
		return matrix;
	}
	
	@Override
	public SymmetricMatrix elementWisePower(double power) {
		SymmetricMatrix matrix = new SymmetricMatrix(m_iRows);
		for (int i = 0; i < matrix.m_iRows; i++) {
			for (int j = i; j < matrix.m_iCols; j++) {
				matrix.setValueAt(i, j, Math.pow(getValueAt(i, j), power));
			}
		}
		return matrix;
	}
	

	@Override
	public SymmetricMatrix getAbsoluteValue() {
		SymmetricMatrix oMat = new SymmetricMatrix(m_iRows);
		for (int i = 0; i < m_iRows; i++) {
			for (int j = i; j < m_iCols; j++) {
				oMat.setValueAt(i, j, Math.abs(getValueAt(i, j)));
			}
		}
		return oMat;
	}

	
	/**
	 * Calculate the Kronecker product of this by the m Matrix object.
	 * @param m a Matrix instance
	 * @return the resulting product (a matrix instance)
	 */
	public final Matrix getKroneckerProduct(Matrix m) {
		if (m.isSymmetric()) {
			SymmetricMatrix result = new SymmetricMatrix(m_iRows * m.m_iRows);
			for (int i1 = 0; i1 < m_iRows; i1++) {
				for (int j1 = i1; j1 < m_iCols; j1++) {
					for (int i2 = 0; i2 < m.m_iRows; i2++) {
						for (int j2 = 0; j2 < m.m_iCols; j2++) {
							result.setValueAt(i1 * m.m_iRows + i2, j1 * m.m_iCols + j2, getValueAt(i1, j1) * m.getValueAt(i2, j2));
						}
					}
				}
			}
			return result;
		} else {
			return super.getKroneckerProduct(m);
		}
	}
	
	/**
	 * Create a Matrix that corresponds to the Isserlis theorem given that matrix this is
	 * a variance-covariance matrix.
	 * @return a SymmetricMatrix instance
	 */
	public final SymmetricMatrix getIsserlisMatrix() {
		SymmetricMatrix output = new SymmetricMatrix(m_iRows * m_iRows);
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
//						if (indexRow != indexCol) {
//							output.setValueAt(indexCol, indexRow, covariance);
//						}
					}
				}
			}
		}

		return output;
	}

	
	/**
	 * This method compute the lower triangle of the Cholesky decomposition.
	 * Checks are implemented to make sure that this is square and symmetric.
	 * @return the resulting matrix 
	 * @throws UnsupportedOperationException if the Cholesky factorisation cannot be completed
	 */
    public Matrix getLowerCholTriangle() {
    	int m1Row = m_iRows;
    	Matrix matrix = new Matrix(m1Row,m1Row);
    	double dTmp;
    	for (int i = 0; i < m1Row; i++) {
    		for (int j = 0; j <= i; j++) {
    			if (j == i) {
    				dTmp = 0;
    				for (int k = 0; k <= i - 1; k++) {
    					dTmp += matrix.getValueAt(i, k) * matrix.getValueAt(i, k);
    				}
    				matrix.setValueAt(i, j, Math.sqrt(getValueAt(i, j) - dTmp));
    			} else {
    				dTmp = 0;
    				for (int k = 0; k <= j - 1; k++) {
    					dTmp += matrix.getValueAt(i, k) * matrix.getValueAt(j, k);
    				}
    				matrix.setValueAt(i, j, 1d / matrix.getValueAt(j, j) * (getValueAt(i, j) - dTmp));
    			}
    			if (Double.isNaN(matrix.getValueAt(i, j))) {
    				throw new UnsupportedOperationException("Matrix.lowerChol(): the lower triangle of the Cholesky decomposition cannot be calculated because NaN have been generated!");
    			}
    		}
    	}
    	return matrix;
    }

	
	/**
	 * Check if the matrix is positive definite. The check is based on the Cholesky factorization. If the factorization can
	 * be computed the method returns true.
	 * @return true if it is or false otherwise
	 */
	public final boolean isPositiveDefinite() {
		try {
			getLowerCholTriangle();
			return true;
		} catch (UnsupportedOperationException e) {
			return false;
		}
	}

	@Override
	public SymmetricMatrix getDeepClone() {
		SymmetricMatrix oMat = new SymmetricMatrix(m_iRows);
		for (int i = 0; i < m_iRows; i++) {
			for (int j = i; j < m_iCols; j++) {
				oMat.setValueAt(i, j, getValueAt(i, j));
			}
		}
		return oMat;
	}
	
	@Override
	protected SymmetricMatrix getInternalInverseMatrix() {
		Matrix m = super.getInternalInverseMatrix();
		return forceConversionToSymmetricMatrix(m);
	}

	/**
	 * Try to convert a Matrix instance to a SymmetricMatrix instance.
	 * @param m a Matrix instance
	 * @return a SymmetricMatrix instance
	 */
	public static SymmetricMatrix convertToSymmetricIfPossible(Matrix m) {
		if (m instanceof SymmetricMatrix) {
			return (SymmetricMatrix) m;
		} else if (!m.isSymmetric())
			throw new UnsupportedOperationException("The Matrix instance m is not symmetric!");
		else {
			return forceConversionToSymmetricMatrix(m);
		}
	}
	
	static SymmetricMatrix forceConversionToSymmetricMatrix(Matrix m) {
		if (m instanceof SymmetricMatrix) {
			return (SymmetricMatrix) m;
		} else {
			if (!m.isSquare()) {
				throw new UnsupportedOperationException("Matrix m must be a square matrix!");
			} else {
				SymmetricMatrix sm = new SymmetricMatrix(m.m_iRows);
				for (int i = 0; i < m.m_iRows; i++) {
					for (int j = i; j < m.m_iRows; j++) {
						sm.setValueAt(i, j, m.getValueAt(i, j));
					}
				}
				return sm;
			}
		}
		
	}
	
	
	@Override
	public Matrix matrixDiagBlock(Matrix m) {
		if (m.isSymmetric()) {
			int m1Row = m_iRows;
			int m2Row = m.m_iRows;
			SymmetricMatrix matrix = new SymmetricMatrix(m1Row + m2Row);
			for (int i = 0; i < m1Row; i++) {
				for (int j = i; j < m1Row; j++) {
					matrix.setValueAt(i, j, this.getValueAt(i, j));
				}
			}
			for (int i = 0; i < m2Row; i++) {
				for (int j = i; j < m2Row; j++) {
					matrix.setValueAt(i + m1Row, j + m1Row, m.getValueAt(i, j));
				}
			}
			return matrix;
		} else {
			return super.matrixDiagBlock(m);
		}
	}

	
	@Override
	public SymmetricMatrix getInverseMatrix() {
		if (isDiagonalMatrix()) {		// procedure for diagonal matrices
			return DiagonalMatrix.forceConversionToDiagonalMatrix(this).getInverseMatrix();
		} 
		List<List<Integer>> indices = getBlockConfiguration();
		if (indices.size() == 1) {
			return getInternalInverseMatrix();
		} else {
			SymmetricMatrix inverseMatrix = new SymmetricMatrix(m_iRows);
			for (List<Integer> blockIndex : indices) {
				Matrix invSubMatrix = getSubMatrix(blockIndex, blockIndex).getInternalInverseMatrix();
				for (int i = 0; i < blockIndex.size(); i++) {
					for (int j = i; j < blockIndex.size(); j++) {
						inverseMatrix.setValueAt(blockIndex.get(i), blockIndex.get(j), invSubMatrix.getValueAt(i, j));
					}
				}
			}
			return inverseMatrix;
		}
		
	}
	
}
