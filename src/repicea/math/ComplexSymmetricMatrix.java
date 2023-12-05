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

/**
 * A Hermitian matrix implementation.<p>
 * 
 * The implementation is based on a full array. It overrides the method
 * setValueAt to make sure the element j,i is equal to element i,j.
 * @author Mathieu Fortin - November 2023
 */
public class ComplexSymmetricMatrix extends ComplexMatrix {

	/**
	 * Constructor. <p>
	 * 
	 * The ComplexMatrix argument is checked for symmetry. If the
	 * check passes the ComplexSymmetricMatrix object is constructed.
	 * 
	 * @param mat a ComplexMatrix instance
	 * 
	 * @throws UnsupportedOperationException if the ComplexMatrix argument is not symmetric
	 */
	public ComplexSymmetricMatrix(ComplexMatrix mat) {
		this(mat.m_iRows);
		if (!mat.isSymmetric()) {
			throw new UnsupportedOperationException("The ComplexMatrix argument must be a symmetric matrix!");
		}
		for (int i = 0; i < mat.m_iRows; i++) {
			for (int j = i; j < mat.m_iCols; j++) {
				setValueAt(i, j, mat.getValueAt(i,j));
			}
		}
	}

	private ComplexSymmetricMatrix(int size) {
		super(size, size);
	}

	@Override
	protected ComplexNumber[][] contructInternalArray(int iRows, int iCols) {
		ComplexNumber[][] mainArray = new ComplexNumber[iRows][];
		for (int i = 0; i < mainArray.length; i++) {
			mainArray[i] = new ComplexNumber[iCols - i]; 
		}
		return mainArray;
	}

	@Override
	public void setValueAt(int i, int j, ComplexNumber value) {
		if (j >= i) {
			m_afData[i][j - i] = value;
		} else {
			m_afData[j][i - j] = value;
		}
	}

	/**
	 * Return the value at row i and column j.
	 * @param i the row index
	 * @param j the column index
	 * @return the entry
	 */
	@Override
	public ComplexNumber getValueAt(int i, int j) {
		return j >= i ? m_afData[i][j - i] : m_afData[j][i - j];
	}

	@Override
	public final boolean isSymmetric() {
		return true;
	}

}
