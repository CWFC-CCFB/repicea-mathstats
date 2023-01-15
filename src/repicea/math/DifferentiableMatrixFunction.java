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
package repicea.math;

/**
 * This interface returns the first and second derivatives of a matrix function.
 * @author Mathieu Fortin - November 2012
 * @param <ParameterID> an enum that defines the parameter index
 */
public interface DifferentiableMatrixFunction<ParameterID> {
	
	/**
	 * This method returns a vector that contains the first derivatives of the function with respect to its parameters.
	 * @param parameter a parameter with respect to the matrix function has to be derived
	 * @return a Matrix instance
	 */
	public Matrix getGradient(ParameterID parameter);
	
	
	/**
	 * This method returns a matrix that contains the second derivatives of the matrix function with respect 
	 * to its parameters.
	 * @param parameter1 the index of the first parameter 
	 * @param parameter2 the index of the second parameter
	 * @return a Matrix instance
	 */
	public Matrix getHessian(ParameterID parameter1, ParameterID parameter2);


}
