/*
 * This file is part of the repicea library.
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
package repicea.stats.model;

import repicea.math.MathematicalFunction;
import repicea.math.Matrix;

/**
 * An interface that makes sure the instance can be interpreted as a likelihood
 * @author Mathieu Fortin - 2012
 */
public interface Likelihood extends MathematicalFunction {
	
	/**
	 * This method sets the vector of observed values.
	 * @param yVector a row vector (Matrix instance)
	 */
	public void setYVector(Matrix yVector);

	/**
	 * This method returns the vector of observed values.
	 * @return a Matrix instance
	 */
	public Matrix getYVector();

	/**
	 * This method returns the prediction associated with the observation.
	 * @return a Matrix instance
	 */
	public Matrix getPredictionVector();
	
}
