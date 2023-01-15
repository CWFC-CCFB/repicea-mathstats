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
package repicea.stats;

import java.io.Serializable;

import repicea.math.AbstractMathematicalFunction;
import repicea.math.Matrix;
import repicea.math.SymmetricMatrix;

/**
 * This class specifies a linear term in the form x * beta.
 * @author Mathieu Fortin - October 2011
 */
@SuppressWarnings("serial")
public final class LinearStatisticalExpression extends AbstractMathematicalFunction implements Serializable {

	protected Matrix gradient;
	protected SymmetricMatrix hessian;
	
	@Override
	public Double getValue() {
		if (getNumberOfParameters() != getNumberOfVariables()) {
			throw new IllegalArgumentException("Incompatible vectors");
		} 
		double productResult = 0;
		for (int i = 0; i < getNumberOfParameters(); i++) {
			productResult += getVariableValue(i) * getParameterValue(i);
		}
		return productResult;
	}

	@Override
	public Matrix getGradient() {
		if (gradient == null || gradient.m_iRows != getNumberOfVariables()) {			// create a gradient matrix only once or only if the number of variables in x changes
			gradient = new Matrix(getNumberOfVariables(),1);								
		}
		
		for (int i = 0; i < getNumberOfVariables(); i++) {
			gradient.setValueAt(i, 0, getVariableValue(i));							// update the value in the gradient matrix
		}
		
		return gradient;
	}

	@Override
	public SymmetricMatrix getHessian() {
		if (hessian == null || hessian.m_iCols != getNumberOfVariables()) {				// create a hessian matrix only once or only if the number of variables in x changes
			hessian = new SymmetricMatrix(getNumberOfVariables());
		}
		return hessian;
	}
	

	
}
