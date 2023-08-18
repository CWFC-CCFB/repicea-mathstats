/*
 * This file is part of the repicea-mathstats library.
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

import java.security.InvalidParameterException;
import java.util.ArrayList;
import java.util.List;

import repicea.math.Matrix;
import repicea.stats.data.DataSet;
import repicea.stats.data.GenericStatisticalDataStructure;
import repicea.stats.data.StatisticalDataException;
import repicea.stats.data.StatisticalDataStructure;
import repicea.stats.estimators.Estimator;
import repicea.stats.estimators.MaximumLikelihoodEstimator;
import repicea.stats.estimators.MaximumLikelihoodEstimator.MaximumLikelihoodCompatibleModel;
import repicea.stats.estimators.OLSEstimator;
import repicea.stats.estimators.OLSEstimator.OLSCompatibleModel;
import repicea.stats.model.AbstractStatisticalModel;
import repicea.stats.model.CompositeLogLikelihood;
import repicea.stats.model.CompositeLogLikelihoodWithExplanatoryVariables;
import repicea.stats.model.PredictableModel;

/**
 * The LinearModel is a traditional model fitted with an Ordinary Least Squares estimator.
 * @author Mathieu Fortin - November 2012
 */
public class LinearModel extends AbstractStatisticalModel implements PredictableModel, OLSCompatibleModel, MaximumLikelihoodCompatibleModel {

	private final StatisticalDataStructure dataStruct;
	private Matrix matrixX;
	private Matrix vectorY;
	private double convergenceCriterion = 1E-8;
	protected final Matrix startingValues;

	
	/**
	 * Constructor for OLS estimation.
	 * @param dataSet a DataSet instance
	 * @param modelDefinition a model definition
	 */
	public LinearModel(DataSet dataSet, String modelDefinition) {
		super();
		checkNonNullValues(dataSet, modelDefinition);
		dataStruct = new GenericStatisticalDataStructure(dataSet);
		startingValues = null;
		try {
			setModelDefinition(modelDefinition);
		} catch (StatisticalDataException e) {
			System.out.println("Unable to define this model : " + modelDefinition);
			e.printStackTrace();
		}
	}

	/**
	 * Constructor for maximum likelihood estimation.
	 * @param dataSet a DataSet instance
	 * @param modelDefinition a model definition
	 * @param startingValues a Matrix of starting values
	 */
	public LinearModel(DataSet dataSet, String modelDefinition, Matrix startingValues) {
		super();
		checkNonNullValues(dataSet, modelDefinition);
		dataStruct = new GenericStatisticalDataStructure(dataSet);
		if (startingValues == null) {
			throw new InvalidParameterException("The parameter startingValues cannot be non null!");
		}
		if (!startingValues.isColumnVector()) {
			throw new InvalidParameterException("The parameter startingValues must be a column vector!");
		}
		this.startingValues = startingValues.getDeepClone();
		try {
			setModelDefinition(modelDefinition);
			if (startingValues.m_iRows != matrixX.m_iCols + 1) {
				throw new InvalidParameterException("Considering the model definition, the parameter startingValues should be a column vector of " + (matrixX.m_iCols + 1) + "elements. The last of them should be the residual variance!");
			}
		} catch (StatisticalDataException e) {
			System.out.println("Unable to define this model : " + modelDefinition);
			e.printStackTrace();
		}
	}

	private static void checkNonNullValues(DataSet dataSet, String modelDefinition) {
		if (dataSet == null) {
			throw new InvalidParameterException("The parameter dataSet cannot be null!");
		}
		if (modelDefinition == null) {
			throw new InvalidParameterException("The parameter modelDefinition cannot be null!");
		}
	}
	/*
	 * Useless (non-Javadoc)
	 * @see repicea.stats.model.StatisticalModel#setParameters(repicea.math.Matrix)
	 */
	@Override
	public void setParameters(Matrix beta) {}


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
	public String toString() {
		return this.getEstimator() instanceof OLSEstimator ? 
				"Linear model fitted with OLS estimator" : 
					"Linear model fitted with maximum likelihood estimator";
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
		return startingValues == null ? 
				new OLSEstimator(this) :
					new MaximumLikelihoodEstimator(this);
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
	public List<String> getOtherParameterNames() {
		if (this.getEstimator() instanceof MaximumLikelihoodEstimator) {
			List<String> otherParms = new ArrayList<String>();
			otherParms.add("sigma2");
			return otherParms;
		} else {
			return new ArrayList<String>();
		}
	}


	@Override
	public double getConvergenceCriterion() {return convergenceCriterion;}

	@Override
	public CompositeLogLikelihood getCompleteLogLikelihood() {
		return new CompositeLogLikelihoodWithExplanatoryVariables(new GaussianIndividualLogLikelihood(getMatrixX().m_iCols, startingValues), getMatrixX(), getVectorY());
	}
	
}
