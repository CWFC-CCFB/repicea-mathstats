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
 * This interface ensures the model can provide its predicted values and its residuals
 * @author Mathieu Fortin - 2022
 */
public interface PredictableModel {

	
	/**
	 * Produce a column vector of predicted values.<p>
	 * These values are calculated from the sample.
	 * @return a Matrix instance
	 * @throws UnsupportedOperationException if the estimator has not converged
	 */
	public default Matrix getPredicted() throws UnsupportedOperationException {
		return getPredicted(null);
	}

	/**
	 * Produce a column vector of predicted values.<p>
	 * These values are calculated using the xMatrix argument.
	 * If this argument is null, then the predictions are calculated from 
	 * the sample.
	 * @param xMatrix a design matrix
	 * @return a Matrix instance
	 * @throws UnsupportedOperationException if the estimator has not converged
	 */
	public Matrix getPredicted(Matrix xMatrix) throws UnsupportedOperationException;
	
	
	/**
	 * Return a column vector of residuals, that is observed values minus predictions.
	 * @return a Matrix instance
	 */
	public Matrix getResiduals();

}
