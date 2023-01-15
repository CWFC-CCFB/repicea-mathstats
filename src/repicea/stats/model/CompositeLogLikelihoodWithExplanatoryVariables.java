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

import repicea.math.Matrix;

/**
 * The log-likelihood of a model given a vector of the response variable and a matrix
 * of explanatory variables.
 * @author Mathieu Fortin - July 2022
 */
@SuppressWarnings("serial")
public class CompositeLogLikelihoodWithExplanatoryVariables extends SimpleCompositeLogLikelihood {

	protected Matrix xValues;
	
	public CompositeLogLikelihoodWithExplanatoryVariables(IndividualLogLikelihood innerLogLikelihoodFunction, Matrix xValues, Matrix yValues) {
		super(innerLogLikelihoodFunction, yValues);
		this.xValues = xValues;
	}
		
	protected void setValuesInLikelihoodFunction(int index) {
		super.setValuesInLikelihoodFunction(index);
		getOriginalFunction().setVariables(xValues.getSubMatrix(index, index, 0, xValues.m_iCols - 1));
	}

	/**
	 * This method returns all the predicted values.
	 * @return a Matrix instance
	 */
	public Matrix getPredictions() {
		Matrix predictedValues = new Matrix(xValues.m_iRows, 1);
		for (int i = 0; i < xValues.m_iRows; i++) {
			setValuesInLikelihoodFunction(i);
			predictedValues.setSubMatrix(getOriginalFunction().getPredictionVector(), i, 0);
		}
		return predictedValues;
	}
	
	/**
	 * Resets this composite likelihood to its initial values.
	 */
	public void reset() {}

	
}
