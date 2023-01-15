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
package repicea.stats.model;

import repicea.math.AbstractMathematicalFunctionWrapper;
import repicea.math.Matrix;
import repicea.math.SymmetricMatrix;
import repicea.math.utility.MatrixUtility;

/**
 * A simple composite log likelihood for distribution models. <br>
 * <br>
 * There is only a vector of the response variable and no explanatory variables.
 * @author Mathieu Fortin - July 2022
 */
@SuppressWarnings("serial")
public class SimpleCompositeLogLikelihood extends AbstractMathematicalFunctionWrapper implements CompositeLogLikelihood {

	private final Matrix yValues;
	
	public SimpleCompositeLogLikelihood(IndividualLogLikelihood innerLogLikelihoodFunction, Matrix yValues) {
		super(innerLogLikelihoodFunction);
		this.yValues = yValues;
	}
		
	@Override
	public IndividualLogLikelihood getOriginalFunction() {return (IndividualLogLikelihood) super.getOriginalFunction();}
	
	@Override
	public Double getValue() {
		double loglikelihood = 0;
		for (int i = 0; i < yValues.m_iRows; i++) {
			setValuesInLikelihoodFunction(i);
			loglikelihood += getOriginalFunction().getValue();
		}
		return loglikelihood;
	}

	@Override
	public Matrix getGradient() {
		Matrix resultingGradient = new Matrix(getOriginalFunction().getNumberOfParameters(), 1);
		for (int i = 0; i < yValues.m_iRows; i++) {
			setValuesInLikelihoodFunction(i);
//			MatrixUtility.add(resultingGradient, getOriginalFunction().getGradient());
			resultingGradient = resultingGradient.add(getOriginalFunction().getGradient());
		}
		return resultingGradient;
	}

	@Override
	public SymmetricMatrix getHessian() {
		SymmetricMatrix resultingHessian = new SymmetricMatrix(getOriginalFunction().getNumberOfParameters());
		for (int i = 0; i < yValues.m_iRows; i++) {
			if (i == 208) {
				int u = 0;
			}
			setValuesInLikelihoodFunction(i);
			Matrix hessianToAdd = getOriginalFunction().getHessian();
//			if (hessianToAdd.anyElementNaN()) {
//				throw new UnsupportedOperationException("The hessian contains NaN!");
//			}
			resultingHessian = (SymmetricMatrix) resultingHessian.add(hessianToAdd);
		}
		return resultingHessian;
	}

	protected void setValuesInLikelihoodFunction(int index) {
		getOriginalFunction().setYVector(yValues.getSubMatrix(index, index, 0, 0));
	}

		
	@Override
	public void setParameters(Matrix beta) {
		getOriginalFunction().setParameters(beta);
	}

	@Override
	public Matrix getParameters() {return getOriginalFunction().getParameters();}

	@Override
	public void reset() {}
	
	
}
