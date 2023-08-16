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
package repicea.stats.model.lm;

import java.util.ArrayList;
import java.util.List;

import repicea.math.Matrix;
import repicea.stats.data.DataSet;
import repicea.stats.data.GenericStatisticalDataStructure;
import repicea.stats.data.StatisticalDataException;
import repicea.stats.data.StatisticalDataStructure;
import repicea.stats.estimators.Estimator;
import repicea.stats.estimators.OLSEstimator;
import repicea.stats.estimators.OLSEstimator.OLSCompatibleModel;
import repicea.stats.model.AbstractStatisticalModel;
import repicea.stats.model.PredictableModel;

/**
 * The LinearModel is a traditional model fitted with an Ordinary Least Squares estimator.
 * @author Mathieu Fortin - November 2012
 */
public class LinearModel extends AbstractStatisticalModel implements PredictableModel, OLSCompatibleModel {

	private final StatisticalDataStructure dataStruct;
	private Matrix matrixX;
	private Matrix vectorY;
	
	/**
	 * Constructor.
	 * @param dataSet a DataSet instance
	 * @param modelDefinition a model definition
	 */
	public LinearModel(DataSet dataSet, String modelDefinition) {
		super();
		dataStruct = new GenericStatisticalDataStructure(dataSet);

		try {
			setModelDefinition(modelDefinition);
		} catch (StatisticalDataException e) {
			System.out.println("Unable to define this model : " + modelDefinition);
			e.printStackTrace();
		}
	}

	
	/*
	 * Useless (non-Javadoc)
	 * @see repicea.stats.model.StatisticalModel#setParameters(repicea.math.Matrix)
	 */
	@Override
	public void setParameters(Matrix beta) {}

//	@Override
//	public Matrix getParameters() {
//		return getEstimator().getParameterEstimates().getMean();
//	}
	
	/**
	 * This method returns the residual variance only if the optimizer is an instance
	 * of OLSOptimizer.
	 * @return a Matrix with a single element
	 */
	public double getResidualVariance() {
		if (getEstimator() instanceof OLSEstimator) {
			return ((OLSEstimator) getEstimator()).getResidualVariance().getMean().getValueAt(0, 0);
		} else {
			return -1d;
		}
	}

	@Override
	public Matrix getPredicted() {
		return getMatrixX().multiply(getParameters());
	}

	@Override
	public Matrix getResiduals() {
		return getVectorY().subtract(getPredicted());
	}

	@Override
	protected Estimator instantiateDefaultEstimator() {
		return new OLSEstimator(this);
	}

	protected void setModelDefinition(String modelDefinition) throws StatisticalDataException {
		super.setModelDefinition(modelDefinition);
		dataStruct.setModelDefinition(modelDefinition);
		vectorY = dataStruct.constructVectorY();
		matrixX = dataStruct.constructMatrixX();
	}


	@Override
	public Matrix getMatrixX() {return matrixX;}


	@Override
	public Matrix getVectorY() {return vectorY;}


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
