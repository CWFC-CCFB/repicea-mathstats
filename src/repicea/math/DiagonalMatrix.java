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

public class DiagonalMatrix extends SymmetricMatrix {

	public DiagonalMatrix(int size) {
		super(size);
	}

	@Override
	protected double[][] contructInternalArray(int iRows, int iCols) {
		double[][] mainArray = new double[1][iRows];
		return mainArray;
	}

	@Override
	public boolean anyElementNaN() {
		for (int i = 0; i < this.m_iRows; i++) {
			if (Double.isNaN(this.getValueAt(i, i)))
				return true;
		}
		return false;
	}

	
	@Override
	public final void setValueAt(int i, int j, double value) {
		if (j != i) {
			throw new UnsupportedOperationException("The DiagonalMatrix instance only allows for setting the values on the diagonal!");
		} else {
			m_afData[0][i] = value;
		}
	}
	
	/**
	 * Return the value at row i and column j.
	 * @param i
	 * @param j
	 * @return a double
	 */
	@Override
	public final double getValueAt(int i, int j) {
		if (i < 0 || i >= this.m_iRows) {
			throw new UnsupportedOperationException("The index i exceeds the dimension of the Matrix instance!");
		}
		if (j < 0 || j >= this.m_iCols) {
			throw new UnsupportedOperationException("The index i exceeds the dimension of the Matrix instance!");
		}
		return j == i ? m_afData[0][i] : 0d;
	}

	@Override
	public final boolean isDiagonalMatrix() {
		return true;
	}
	
	@Override
	public final DiagonalMatrix elementWiseDivide(Matrix m) {
		if (isTheSameDimension(m)) {
			DiagonalMatrix oMat = new DiagonalMatrix(m_iRows);
			for (int i = 0; i < this.m_iRows; i++) {
				oMat.setValueAt(i, i, getValueAt(i, i) / m.getValueAt(i, i));
			}
			return oMat;
		} else {
			throw new UnsupportedOperationException("The matrix m does not have the same dimensions than the current matrix!");
		}
	}

	@Override
	public final DiagonalMatrix elementWiseMultiply(Matrix m) {
		if (isTheSameDimension(m)) {
			DiagonalMatrix oMat = new DiagonalMatrix(m_iRows);
			for (int i = 0; i < this.m_iRows; i++) {
				oMat.setValueAt(i, i, getValueAt(i, i) * m.getValueAt(i, i));
			}
			return oMat;
		} else {
			throw new UnsupportedOperationException("The matrix m does not have the same dimensions than the current matrix!");
		}
	}

	@Override
    public final SymmetricMatrix logMatrix() {
		throw new UnsupportedOperationException("The DiagonalMatrix class does not support the logMatrix method!");
    }


	@Override
	public final DiagonalMatrix scalarMultiply(double d) {
		DiagonalMatrix mat = new DiagonalMatrix(m_iRows);
		for (int i = 0; i < m_iRows; i++) {
			mat.setValueAt(i, i, getValueAt(i, i) * d);
		}
		return mat;
	}
	
	@Override
	public final DiagonalMatrix elementWisePower(double power) {
		DiagonalMatrix matrix = new DiagonalMatrix(m_iRows);
		for (int i = 0; i < matrix.m_iRows; i++) {
			matrix.setValueAt(i, i, Math.pow(getValueAt(i, i), power));
		}
		return matrix;
	}
	

	@Override
	public final DiagonalMatrix getAbsoluteValue() {
		DiagonalMatrix oMat = new DiagonalMatrix(m_iRows);
		for (int i = 0; i < m_iRows; i++) {
			oMat.setValueAt(i, i, Math.abs(getValueAt(i, i)));
		}
		return oMat;
	}

	

	/**
	 * This method compute the lower triangle of the Cholesky decomposition.
	 * Checks are implemented to make sure that this is square and symmetric.
	 * @return the resulting matrix 
	 * @throws UnsupportedOperationException if the Cholesky factorisation cannot be completed
	 */
    public final DiagonalMatrix getLowerCholTriangle() {
    	DiagonalMatrix matrix = new DiagonalMatrix(m_iRows);
    	for (int i = 0; i < m_iRows; i++) {
    		matrix.setValueAt(i, i, Math.sqrt(getValueAt(i, i)));
    	}
    	return matrix;
    }

	@Override
	public final DiagonalMatrix getDeepClone() {
		DiagonalMatrix oMat = new DiagonalMatrix(m_iRows);
		for (int i = 0; i < m_iRows; i++) {
			oMat.setValueAt(i, i, getValueAt(i, i));
		}
		return oMat;
	}
	
	@Override
	public final Matrix multiply(Matrix m) {
		if (m_iCols != m.m_iRows) {
			throw new UnsupportedOperationException("The matrix m cannot multiply the current matrix for the number of rows is incompatible!");
		} else {
			if (m.equals(this)) {	// multiplied by itself yields a SymmetricMatrix instance
				DiagonalMatrix mat = new DiagonalMatrix(m_iRows);
				for (int i = 0; i < m_iRows; i++) {
					double originalValue = getValueAt(i,i);
					mat.setValueAt(i, i, originalValue * originalValue);
				}
				return mat;
			} else {
				return super.multiply(m);
			}
		}
	}

	
	@Override
	protected final DiagonalMatrix getInternalInverseMatrix() {
		DiagonalMatrix m = this.elementWisePower(-1);
		return m;
	}

	
	@Override
	public final Matrix matrixDiagBlock(Matrix m) {
		if (m.isDiagonalMatrix()) {
			int m1Row = m_iRows;
			int m2Row = m.m_iRows;
			DiagonalMatrix matrix = new DiagonalMatrix(m1Row + m2Row);
			for (int i = 0; i < m1Row; i++) {
				matrix.setValueAt(i, i, this.getValueAt(i, i));
			}
			for (int i = 0; i < m2Row; i++) {
				matrix.setValueAt(i + m1Row, i + m1Row, m.getValueAt(i, i));
			}
			return matrix;
		} else {
			return super.matrixDiagBlock(m);
		}
	}

	
	@Override
	public final DiagonalMatrix getInverseMatrix() {
		return getInternalInverseMatrix();
	}

	static DiagonalMatrix forceConversionToDiagonalMatrix(Matrix m) {
		if (m instanceof DiagonalMatrix) {
			return (DiagonalMatrix) m;
		} else {
			if (!m.isDiagonalMatrix()) {
				throw new UnsupportedOperationException("Some off diagonal elements are different from 0!");
			} else {
				DiagonalMatrix sm = new DiagonalMatrix(m.m_iRows);
				for (int i = 0; i < m.m_iRows; i++) {
					sm.setValueAt(i, i, m.getValueAt(i, i));
				}
				return sm;
			}
		}
		
	}

	
}
