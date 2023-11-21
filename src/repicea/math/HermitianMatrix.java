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
 * setValueAt to make sure the element j,i is the conjugate of element i,j.
 * @author Mathieu Fortin - November 2023
 */
public class HermitianMatrix extends ComplexMatrix {

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
	public HermitianMatrix(ComplexMatrix mat) {
		this(mat.m_iRows);
		if (!mat.isSquare()) {
			throw new UnsupportedOperationException("The ComplexMatrix argument must be a square matrix!");
		}
		for (int i = 0; i < mat.m_iRows; i++) {
			for (int j = i; j < mat.m_iCols; j++) {
				ComplexNumber valueUp = mat.getValueAt(i,j);
				ComplexNumber valueBelow = mat.getValueAt(j,i);
				if (valueUp.equals(valueBelow.getComplexConjugate())) {
					setValueAt(i,j, valueUp);
				} else {
					throw new UnsupportedOperationException("The ComplexMatrix argument is not a Hermitian matrix!");
				}
			}
		}
	}

	private HermitianMatrix(int size) {
		super(size, size);
	}

	@Override
	protected ComplexNumber[][] contructInternalArray(int iRows, int iCols) {
		ComplexNumber[][] mainArray = new ComplexNumber[iRows][iCols];
		return mainArray;
	}

	@Override
	public void setValueAt(int i, int j, ComplexNumber value) {
		m_afData[i][j] = value;
		if (i != j) {
			m_afData[j][i] = value.getComplexConjugate();
		}
	}

}
