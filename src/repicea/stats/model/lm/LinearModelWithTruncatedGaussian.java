/*
 * This file is part of the repicea-mathstats library.
 *
 * Copyright (C) 2023 His Majesty the King in right of Canada
 * Author: Mathieu Fortin, Canadian Wood Fibre Centre, CFS
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
package repicea.stats.model.lm;

import java.util.ArrayList;
import java.util.List;

import repicea.math.Matrix;
import repicea.stats.data.DataSet;
import repicea.stats.data.GenericStatisticalDataStructure;
import repicea.stats.data.StatisticalDataException;
import repicea.stats.data.StatisticalDataStructure;
import repicea.stats.estimators.MaximumLikelihoodEstimator;
import repicea.stats.estimators.MaximumLikelihoodEstimator.MaximumLikelihoodCompatibleModel;
import repicea.stats.model.AbstractStatisticalModel;
import repicea.stats.model.CompositeLogLikelihood;

// TODO check if the PredictableModel interface is really needed here 
public class LinearModelWithTruncatedGaussian extends AbstractStatisticalModel implements MaximumLikelihoodCompatibleModel {

	private final StatisticalDataStructure dataStruct;
	private Matrix matrixX;
	private Matrix vectorY;

	/**
	 * Constructor.
	 * @param dataSet a DataSet instance
	 * @param modelDefinition a model definition
	 */
	public LinearModelWithTruncatedGaussian(DataSet dataSet, String modelDefinition) {
		super();
		dataStruct = new GenericStatisticalDataStructure(dataSet);

		try {
			setModelDefinition(modelDefinition);
		} catch (StatisticalDataException e) {
			System.out.println("Unable to define this model : " + modelDefinition);
			e.printStackTrace();
		}

	}
	

	@Override
	public double getConvergenceCriterion() {
		// TODO Auto-generated method stub
		return 0;
	}

	@Override
	public CompositeLogLikelihood getCompleteLogLikelihood() {
		// TODO Auto-generated method stub
		return null;
	}

	@Override
	protected MaximumLikelihoodEstimator instantiateDefaultEstimator() {
		return new MaximumLikelihoodEstimator(this);
	}
	
	@Override
	public int getNumberOfObservations() {
		return dataStruct.getNumberOfObservations();
	}

	@Override
	public boolean isInterceptModel() {return dataStruct.isInterceptModel();}

	@Override
	public List<String> getEffectList() {return dataStruct.getEffectList();}

	@Override
	public List<String> getOtherParameterNames() {return new ArrayList<String>();}


}
