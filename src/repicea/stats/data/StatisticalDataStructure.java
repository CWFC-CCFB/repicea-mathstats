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
package repicea.stats.data;


import java.util.List;

import repicea.math.Matrix;

/**
 * This interface is the basic interface for any data structure.
 * @author Mathieu Fortin - October 2011
 */
public interface StatisticalDataStructure {

	
	/**
	 * Returns the number of observations in the data structure.
	 * @return an integer
	 */
	public int getNumberOfObservations();
	
	/**
	 * Set the model definition.
	 * @param modelDefinition the model definition (e.g. y ~ x1 + x2 + x3)
	 */
	public void setModelDefinition(String modelDefinition);
	
//	/**
//	 * This method builds the matrices that are required to fit the model. 
//	 * @throws StatisticalDataException
//	 */
//	public void constructMatrices() throws StatisticalDataException;
	
	/**
	 * Return the design matrix of the fixed effects.
	 * @return a Matrix instance
	 */
	public Matrix constructMatrixX();

	/**
	 * Return the vector of response variables.
	 * @return a Matrix instance
	 */
	public Matrix constructVectorY();
	
	/**
	 * This method returns true if the model has an intercept or false otherwise.
	 * @return a boolean instance
	 */
	public boolean isInterceptModel();
	
	/**
	 * This method determines whether the model has an intercept
	 * @param isInterceptModel true if the model has an intercept or false otherwise
	 */
	public void setInterceptModel(boolean isInterceptModel);
	
	
	
	/**
	 * This method returns the DataSet instance behind the StatisticalDataStructure.
	 * @return a DataSet instance
	 */
	public DataSet getDataSet();
	
	
	/**
	 * This method provides the possible value for a particular dummy variable.
	 * @param fieldName the name of the field containing the variable.
	 * @param refClass the name of the class that is the reference class (can be null).
	 * @return a List of values
	 */
	public List getPossibleValueForDummyVariable(String fieldName, String refClass);


	/**
	 * Returns the index of this field in the model definition. The "0" refers to the intercept.
	 * @param effect
	 * @return an integer (-1 if the effect is not found)
	 */
	public int indexOfThisEffect(String effect);

	/**
	 * Returns the list of effects and interactions in the model.
	 * @return
	 */
	public List<String> getEffectList();
}
