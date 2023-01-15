/*
 * This file is part of the repicea library.
 *
 * Copyright (C) 2009-2015 Mathieu Fortin for Rouge-Epicea
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
package repicea.stats.integral;

import java.security.InvalidParameterException;
import java.util.ArrayList;
import java.util.List;

import repicea.math.LogFunctionWrapper;
import repicea.math.MathematicalFunction;
import repicea.math.Matrix;
import repicea.math.SymmetricMatrix;

/**
 * The InternalLogWrapperFunction class is used by the LaplaceApproximation and 
 * AdaptativeGaussHermiteQuadrature classes.
 * @author Mathieu Fortin - December 2015
 */
@SuppressWarnings("serial")
class InternalLogWrapperFunction extends LogFunctionWrapper {

	private final double mPart;
	private final Matrix invG;
	private final Matrix originalParameterValues;
	private final List<Integer> parameterIndices;

	/**
	 * Constructor
	 * @param originalFunction the nested function
	 * @param gMatrix the variance-covariance matrix of the deviate
	 */
	InternalLogWrapperFunction(MathematicalFunction originalFunction, List<Integer> parameterIndices, Matrix gMatrix) {
		super(originalFunction);
		if (!gMatrix.isSymmetric()) {
			throw new InvalidParameterException("Matrix G is supposed to be symmetric!");
		}
		//				this.gMatrix = gMatrix;
		this.invG = gMatrix.getInverseMatrix();
		double n = gMatrix.m_iRows;
		double gDeterminant = gMatrix.getDeterminant();
		mPart = - 0.5 * n * Math.log(2d * Math.PI) - 0.5 * Math.log(gDeterminant);
		this.parameterIndices = parameterIndices;
		originalParameterValues = getParametersFromNestedFunction();
	}

	@Override
	public Double getValue() {
		Matrix u = getParametersFromNestedFunction().subtract(originalParameterValues);
		return super.getValue() + mPart - 0.5 * u.transpose().multiply(invG).multiply(u).getValueAt(0, 0); 
	}

	@Override
	public Matrix getGradient() {
		List<Integer> columnIndex = new ArrayList<Integer>();
		columnIndex.add(0);
		Matrix u = getParametersFromNestedFunction().subtract(originalParameterValues);
		return super.getGradient().getSubMatrix(parameterIndices, columnIndex).subtract(invG.multiply(u));
	}


	@Override
	public SymmetricMatrix getHessian() {
		return SymmetricMatrix.convertToSymmetricIfPossible(super.getHessian().getSubMatrix(parameterIndices, parameterIndices).subtract(invG));
	}

	private Matrix getParametersFromNestedFunction() {
		Matrix parameterValues = new Matrix(parameterIndices.size(), 1);
		for (int i = 0; i < parameterIndices.size(); i++) {
			parameterValues.setValueAt(i, 0, getOriginalFunction().getParameterValue(parameterIndices.get(i)));
		}
		return parameterValues;
	}

}
